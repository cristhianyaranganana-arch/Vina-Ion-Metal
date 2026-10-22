#!/usr/bin/env python3
"""Vina-Ion-Metal receptor preparation (audits #5 and #11).

VinaMetal 1.0.0 shipped a thin passthrough to Meeko's mk_prepare_receptor.py.
That left three practical failure modes:

  1. Meeko cannot type Cu2+/Ni2+ (no templates, not in its metal list), so the
     Cu/Ni parameters of the scoring engine were unreachable in practice.
  2. Meeko is strict about alternate locations, missing atoms and unknown
     HETATM records, and the wrapper offered no preprocessing help.
  3. No pre-flight coordination report: occupancy problems were only visible
     after docking, in the runtime CN warnings.

This script is a self-contained pipeline:

  * PDB preprocessing: multi-model handling, alternate-location resolution,
    optional water/HETATM removal, numbering-gap and missing-backbone checks.
  * Metal detection by element (column 77-78) with name/residue fallbacks that
    disambiguate cases such as CA (calcium ion) vs CA (alpha carbon).
  * PARAMETRIZED ions (Zn, Mg, Ca, Mn, Fe, Cu, Ni) are appended to the final
    PDBQT with their specific AutoDock types - including Cu and Ni, which Meeko
    cannot emit (audit #5).
  * UNPARAMETRIZED metals (Na, K, Co, Cd, Hg, ...) are preserved with standard
    Vina treatment and a warning is printed.
  * Coordination pre-flight report: geometric N/O/S donor count around each
    metal, estimated ligand CN target, and warnings for empty coordination
    spheres (mirrors the runtime diagnostics of audit #10).

Two conversion backends:

  --backend meeko (default when available)
      Metals are stripped, the protein is converted by Meeko's
      mk_prepare_receptor.py (recommended: explicit hydrogens, real charges),
      then the typed metal records are appended to the Meeko output.

  --backend simple (no external dependencies)
      Built-in united-atom converter with standard amino-acid templates.
      Valid for the vina/vinardo scoring functions (charges are ignored
      there); NOT suitable for AD4 map generation. If the input PDB carries
      hydrogens, polar hydrogens are kept so receptor H-bond donors are typed
      correctly; without hydrogens the receptor is acceptor-only (warned).

Examples:

    scripts/prepare_receptor.py --pdb receptor.pdb --out receptor.pdbqt
    scripts/prepare_receptor.py --pdb rec.pdb --out rec.pdbqt --backend simple
    scripts/prepare_receptor.py --pdb rec.pdb --out rec.pdbqt --metal-charge "Cu=2,Fe=3"

Requires: Python 3.6+ (standard library only) for --backend simple;
Meeko (https://github.com/forlilab/Meeko) for the default backend.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

VERSION = "1.0.1"

# ---------------------------------------------------------------------------
# Metal knowledge base (mirrors src/lib/metal_parameters.cpp)
# ---------------------------------------------------------------------------

# element -> (ad_type, ideal_cn, geometry, cn_cutoff, default_formal_charge)
PARAMETRIZED_METALS = {
    "Zn": ("Zn", 4, "tetrahedral",   3.0, 2),
    "Mg": ("Mg", 6, "octahedral",    3.0, 2),
    "Ca": ("Ca", 6, "octahedral",    3.3, 2),
    "Mn": ("Mn", 6, "octahedral",    3.1, 2),
    "Fe": ("Fe", 6, "octahedral",    3.0, 2),
    "Cu": ("Cu", 4, "square-planar", 2.9, 2),
    "Ni": ("Ni", 6, "octahedral",    3.0, 2),
}

# Metals recognized by Vina's generic metal typing (standard Vina behavior,
# no coordination scoring)
UNPARAMETRIZED_METALS = {"Na", "K", "Co", "Cd", "Hg", "U", "Ni4", "Pt", "Se"}
UNPARAMETRIZED_METALS -= set(PARAMETRIZED_METALS)

ALL_METAL_ELEMENTS = set(PARAMETRIZED_METALS) | {"Na", "K", "Co", "Cd", "Hg", "U", "Pt", "Gd", "Sm", "W", "Mo", "V", "Cr", "Ti", "Al", "Li", "Cs", "Rb", "Sr", "Ba"}

DONOR_ELEMENTS = {"N", "O", "S"}

STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE",  # selenomethionine treated as MET-like
    "HSD", "HSE", "HSP", "HIE", "HID", "HIP",  # protonation-specific His variants
}

BACKBONE_ATOMS = {"N", "CA", "C", "O"}

# United-atom AutoDock typing templates for the standard amino acids.
# Maps side-chain atom name -> AD type; backbone handled separately.
# Polar hydrogen names are listed per residue so that hydrogens present in the
# input PDB can be typed HD (donors) while non-polar hydrogens are dropped.
SIDECHAIN_TYPES = {
    "ALA": {}, "GLY": {}, "PRO": {}, "VAL": {}, "LEU": {}, "ILE": {},
    "PHE": {},  # all carbons
    "MET": {"SD": "SA"},
    "MSE": {"SE": "SA"},
    "CYS": {"SG": "SA"},
    "SER": {"OG": "OA"},
    "THR": {"OG1": "OA"},
    "TRP": {"NE1": "N"},  # donor amide N (polar H upgrades to donor)
    "TYR": {"OH": "OA"},
    "ASN": {"OD1": "OA", "ND2": "N"},
    "GLN": {"OE1": "OA", "NE2": "N"},
    "ASP": {"OD1": "OA", "OD2": "OA"},
    "GLU": {"OE1": "OA", "OE2": "OA"},
    "LYS": {"NZ": "N"},  # protonated donor, non-acceptor
    "ARG": {"NE": "N", "NH1": "N", "NH2": "N"},  # guanidinium donors
    "HIS": {"ND1": "NA", "NE2": "NA"},  # acceptors; H presence upgrades to donor+acceptor
    "HSD": {"ND1": "NA", "NE2": "NA"},
    "HSE": {"ND1": "NA", "NE2": "NA"},
    "HSP": {"ND1": "NA", "NE2": "NA"},
    "HIE": {"ND1": "NA", "NE2": "NA"},
    "HID": {"ND1": "NA", "NE2": "NA"},
    "HIP": {"ND1": "NA", "NE2": "NA"},
}

POLAR_HYDROGENS = {
    "ALA": set(), "GLY": set(), "VAL": set(), "LEU": set(), "ILE": set(),
    "PHE": set(), "MET": set(), "PRO": set(),
    "CYS": {"HG"},
    "SER": {"HG"},
    "THR": {"HG1"},
    "TRP": {"HE1"},
    "TYR": {"HH"},
    "ASN": {"HD21", "HD22"},
    "GLN": {"HE21", "HE22"},
    "ASP": set(),
    "GLU": set(),
    "LYS": {"HZ1", "HZ2", "HZ3"},
    "ARG": {"HE", "HH11", "HH12", "HH21", "HH22"},
    "HIS": {"HD1", "HE2"},
    "HSD": {"HD1"}, "HSE": {"HE2"}, "HSP": {"HD1", "HE2"},
    "HID": {"HD1"}, "HIE": {"HE2"}, "HIP": {"HD1", "HE2"},
    "MSE": set(),
}


def warn(msg):
    print(f"WARNING: {msg}", file=sys.stderr)


def info(msg):
    print(f"[prepare_receptor] {msg}", file=sys.stderr)


def normalize_element(el):
    """PDB writes element symbols upper-case in columns 77-78 ('ZN', 'CL');
    normalize to chemical convention ('Zn', 'Cl') for the lookup tables."""
    if not el:
        return ""
    if len(el) == 1:
        return el.upper()
    return el[0].upper() + el[1:].lower()


# ---------------------------------------------------------------------------
# PDB parsing
# ---------------------------------------------------------------------------

class Atom:
    __slots__ = ("record", "serial", "name", "altloc", "resname", "chain",
                 "resseq", "icode", "x", "y", "z", "occ", "bfac", "element",
                 "charge", "ad_type")

    def __init__(self, line):
        self.record = line[0:6].strip()
        self.serial = line[6:11].strip()
        self.name = line[12:16].strip()
        self.altloc = line[16:17].strip()
        self.resname = line[17:20].strip()
        self.chain = line[21:22].strip()
        self.resseq = line[22:26].strip()
        self.icode = line[26:27].strip()
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.occ = float(line[54:60]) if line[54:60].strip() else 1.0
        self.bfac = float(line[60:66]) if line[60:66].strip() else 0.0
        self.element = normalize_element(line[76:78].strip())
        self.charge = 0.0
        self.ad_type = None

    @property
    def res_id(self):
        return (self.chain, self.resseq, self.icode)

    def distance_to(self, other):
        return ((self.x - other.x) ** 2 + (self.y - other.y) ** 2
                + (self.z - other.z) ** 2) ** 0.5


def parse_pdb(path, altloc_mode="first", keep_waters=False, remove_het=False):
    """Parse a PDB file into atom lists.

    Returns (protein_atoms, het_atoms, waters, metals, problems) where metals
    are HETATM/ATOM records identified as metal ions by element.
    """
    protein, het, waters, metals, other_metal_candidates = [], [], [], [], []
    problems = []
    seen_models = False
    in_first_model = True
    altloc_seen = {}

    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            rec = line[0:6].strip() if len(line) >= 6 else ""
            if rec == "MODEL":
                if seen_models:
                    in_first_model = False  # only process the first model
                seen_models = True
                continue
            if rec == "ENDMDL":
                in_first_model = seen_models and in_first_model
                continue
            if not in_first_model:
                continue
            if rec not in ("ATOM", "HETATM"):
                continue
            if len(line) < 54:
                problems.append(f"short {rec} record skipped: {line!r}")
                continue
            try:
                atom = Atom(line)
            except ValueError as exc:
                problems.append(f"unparseable {rec} record skipped ({exc}): {line!r}")
                continue

            # alternate location resolution
            if atom.altloc and atom.altloc != " ":
                key = (atom.name, atom.res_id)
                if key not in altloc_seen:
                    altloc_seen[key] = atom.altloc
                    if altloc_mode not in ("first", atom.altloc):
                        info(f"altLoc '{atom.altloc}' kept for {atom.name} "
                             f"{atom.resname}{atom.resseq} (first-seen policy)")
                elif altloc_mode == "first":
                    continue  # drop other altlocs

            if atom.element == "" or len(atom.element) > 2:
                atom.element = infer_element(atom)

            # Metal-ion cross-check by residue identity: the element column is
            # often blank or misplaced for ion records (audit #11 robustness).
            if atom.element not in ALL_METAL_ELEMENTS:
                guessed = guess_metal_from_resname(atom)
                if guessed:
                    atom.element = guessed

            if atom.element in ALL_METAL_ELEMENTS and is_metal_record(atom):
                metals.append(atom)
            elif atom.resname in STANDARD_AMINO_ACIDS:
                protein.append(atom)
            elif atom.resname in ("HOH", "WAT", "DOD", "TIP3"):
                if keep_waters:
                    waters.append(atom)
            elif atom.record == "HETATM":
                if remove_het:
                    continue
                if is_metal_record(atom) and atom.element in ALL_METAL_ELEMENTS:
                    metals.append(atom)
                else:
                    het.append(atom)
            else:
                protein.append(atom)

    return protein, het, waters, metals, problems


def infer_element(atom):
    """Best-effort element inference when columns 77-78 are blank."""
    name = atom.name
    res = atom.resname
    if res in STANDARD_AMINO_ACIDS:
        if name in BACKBONE_ATOMS:
            return "N" if name == "N" else ("O" if name == "O" else "C")
        t = SIDECHAIN_TYPES.get(res, {}).get(name)
        if t is not None:
            # AD type first letter is the element for C/N/O/S templates
            return t[0]
        if name.startswith(("H", "D")):
            return "H"
        return "C"
    # residue name that IS an element symbol (e.g. '   CA' ion vs ALA's CA atom)
    if res.capitalize() in ALL_METAL_ELEMENTS:
        return res.capitalize()
    if name and name[0].isalpha():
        two = name[0:2].capitalize()
        if two in ALL_METAL_ELEMENTS:
            return two
        return name[0].upper()
    return ""


def is_metal_record(atom):
    """A record is a metal ion when element is a metal AND (residue looks like
    an ion OR the atom name matches the element). Disambiguates e.g. CA ion vs
    alpha carbon (which has element C anyway)."""
    el = atom.element
    if el not in ALL_METAL_ELEMENTS:
        return False
    res = atom.resname.strip().upper()
    name = atom.name.strip().upper()
    if res in STANDARD_AMINO_ACIDS:
        # e.g. a "ZN" atom misnamed inside an amino-acid residue: only a metal
        # if the atom name equals the element symbol
        return name.startswith(el.upper())
    if res in ("HOH", "WAT", "DOD", "TIP3"):
        return False
    return True


def guess_metal_from_resname(atom):
    """Recover metal identity when the element column is missing or misplaced.

    Ion records in the wild frequently look like
        HETATM 8723 NA    NA A 401 ...      (element column blank/wrong)
    A residue whose name is a bare metal symbol (optionally with a valence
    suffix like ZN2/FE3) and whose atom name matches it is an ion. This also
    disambiguates 'CA' the ion from 'CA' the alpha carbon (alanine's resname
    is ALA, not CA)."""
    res = atom.resname.strip().upper()
    name = atom.name.strip().upper()
    res_base = res.rstrip("0123456789+-")
    name_base = name.rstrip("0123456789+-")
    if not res_base or not name_base:
        return None
    element = res_base[0] + res_base[1:].lower()  # 'NA' -> 'Na'
    if element not in ALL_METAL_ELEMENTS:
        return None
    if res_base == name_base or name_base.startswith(res_base):
        return element
    return None


def check_chain_integrity(protein):
    """Report numbering gaps and missing backbone atoms (audit #11)."""
    problems = []
    by_chain = {}
    for a in protein:
        by_chain.setdefault(a.chain, {}).setdefault(
            (a.resseq, a.icode, a.resname), []).append(a)
    for chain, residues in by_chain.items():
        keys = sorted(residues.keys(),
                      key=lambda k: (int(k[0]) if k[0].isdigit() else 0, k[1]))
        prev_num = None
        for key in keys:
            num, icode, resname = key
            try:
                n = int(num)
            except ValueError:
                continue
            if prev_num is not None and n > prev_num + 1:
                problems.append(
                    f"chain {chain or '-'}: numbering gap {prev_num} -> {n} "
                    f"({n - prev_num - 1} residue(s) missing or non-consecutive)")
            names = {a.name for a in residues[key]}
            for bb in ("N", "CA", "C"):
                if bb not in names:
                    problems.append(
                        f"chain {chain or '-'} {resname}{num}: missing backbone atom {bb}")
                    break
            prev_num = n
    return problems


# ---------------------------------------------------------------------------
# Coordination pre-flight report (mirrors the runtime diagnostics, audit #10)
# ---------------------------------------------------------------------------

def coordination_report(metals, protein, het, waters, metal_charges):
    """Geometric N/O/S donor count around each metal and CN target estimate."""
    print("\nCoordination pre-flight report (geometric N/O/S donors per metal):")
    print("-" * 72)
    lines = []
    for m in metals:
        el = m.element
        cutoff = PARAMETRIZED_METALS.get(el, (None, None, None, 3.2, 2))[3]
        ideal_cn = PARAMETRIZED_METALS.get(el, (None, None, None, 3.2, None))[1]
        donors = 0
        donor_atoms = []
        neighbors = list(protein) + list(het) + list(waters)
        for a in neighbors:
            if a.element not in DONOR_ELEMENTS:
                continue
            if m.distance_to(a) <= cutoff:
                donors += 1
                donor_atoms.append(f"{a.name}-{a.resname}{a.resseq}")
        if el in PARAMETRIZED_METALS:
            target = max(0, (ideal_cn or 0) - donors)
            status = "parametrized -> coordination scoring ACTIVE" if donors > 0 \
                else "parametrized -> coordination scoring ACTIVE (empty sphere!)"
            lines.append(
                f"  {el:2s} ({m.resname}{m.resseq} chain {m.chain or '-'}): "
                f"{donors} receptor donor(s) within {cutoff:.1f} A "
                f"-> ligand CN target ~{target} (ideal CN {ideal_cn})")
            if donors == 0:
                lines.append(
                    f"      WARNING: no N/O/S donors within {cutoff:.1f} A. The metal may be"
                    " an artifact, part of a\n"
                    "      partially-modeled site, or residues may be missing. Check the "
                    "structure before docking.")
            elif donor_atoms:
                shown = ", ".join(donor_atoms[:6]) + (" ..." if len(donor_atoms) > 6 else "")
                lines.append(f"      donors: {shown}")
        else:
            lines.append(
                f"  {el:2s} ({m.resname}{m.resseq} chain {m.chain or '-'}): "
                f"UNPARAMETRIZED -> standard Vina interactions (no coordination scoring)")
    for ln in lines:
        print(ln)
    print("-" * 72)


# ---------------------------------------------------------------------------
# PDBQT output
# ---------------------------------------------------------------------------

def pdbqt_line(atom, ad_type, charge, serial):
    """Format one receptor atom in PDBQT (columns per parse_pdbqt_atom_string:
    serial 7-11, x/y/z 31-38/39-46/47-54, charge 69-76, type 78-79)."""
    record = "ATOM  "
    name = atom.name
    if len(name) < 4:
        name = " " + name  # 1-2 letter element symbols: name starts at col 14
    altloc = " "  # altlocs were resolved during parsing
    chain = atom.chain if atom.chain else " "
    try:
        resseq = int(atom.resseq)
    except ValueError:
        resseq = 0
    icode = atom.icode if atom.icode else " "
    # columns: 1-6 record, 7-11 serial, 12 blank, 13-16 name, 17 altloc,
    # 18-20 resname, 21 blank, 22 chain, 23-26 resseq, 27 icode, 28-30 blank,
    # 31-54 xyz, 55-60 occ, 61-66 bfac, 67-76 charge, 77 blank, 78-79 type
    line = (
        f"{record}{serial:5d} {name:<4s}{altloc:1s}{atom.resname:>3s} {chain:1s}"
        f"{resseq:4d}{icode:1s}   "
        f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
        f"{atom.occ:6.2f}{atom.bfac:6.2f}"
        f"{charge:10.3f} {ad_type:>2s}"
    )
    assert len(line) == 79, f"PDBQT line width {len(line)} != 79: {line!r}"
    return line


def simple_convert(protein, het, waters, metals, metal_charges):
    """Built-in united-atom converter (no external dependencies).

    Heavy atoms are typed from standard amino-acid templates; polar hydrogens
    present in the input are kept (typed HD) so Vina can type receptor H-bond
    donors; non-polar hydrogens are dropped (united-atom convention). Charges
    are zero: the vina/vinardo scoring functions ignore receptor charges
    (AD4 map generation is NOT supported by this backend).
    """
    lines = []
    serial = 1
    n_h_kept = 0
    for atom in protein + het + waters:
        res = atom.resname if atom.resname in SIDECHAIN_TYPES else None
        el = atom.element
        if el == "H":
            if res and atom.name in POLAR_HYDROGENS.get(res, set()):
                ad_type = "HD"
                n_h_kept += 1
            elif not res and atom.name in ("H",):  # backbone amide H of unknown residue
                ad_type = "HD"
                n_h_kept += 1
            else:
                continue  # non-polar hydrogen: united-atom merge
        elif res:
            if atom.name in BACKBONE_ATOMS:
                ad_type = {"N": "N", "CA": "C", "C": "C", "O": "OA"}[atom.name]
                if res == "PRO" and atom.name == "N":
                    ad_type = "N"
            else:
                t = SIDECHAIN_TYPES.get(res, {}).get(atom.name)
                if t is None:
                    t = {"C": "C", "N": "N", "O": "OA", "S": "SA", "P": "P",
                         "F": "F", "CL": "Cl", "BR": "Br", "I": "I"}.get(el.upper(), "C")
                ad_type = t
        else:
            ad_type = {"C": "C", "N": "N", "O": "OA", "S": "SA", "P": "P",
                       "F": "F", "CL": "Cl", "BR": "Br", "I": "I"}.get(el.upper(), "C")
        lines.append(pdbqt_line(atom, ad_type, 0.0, serial))
        serial += 1

    lines.extend(metal_pdbqt_lines(metals, metal_charges, serial))
    if n_h_kept == 0:
        warn("input PDB has no (kept) hydrogens: the receptor will be typed "
             "acceptor-only (no H-bond donors). Metal coordination scoring is "
             "unaffected; for full H-bond typing use --backend meeko or a "
             "protonated PDB.")
    return "\n".join(lines) + "\n", serial


def metal_pdbqt_lines(metals, metal_charges, start_serial):
    """Emit metal ion records with the correct AutoDock types (audit #5).

    Parametrized ions (Zn/Mg/Ca/Mn/Fe/Cu/Ni) get their specific AD type so the
    coordination scoring engine picks them up - this is the path Meeko cannot
    produce for Cu/Ni. Unparametrized metals keep their element symbol and the
    standard Vina generic-metal treatment.
    """
    lines = []
    serial = start_serial
    for m in metals:
        el = m.element
        if el in PARAMETRIZED_METALS:
            ad_type = PARAMETRIZED_METALS[el][0]
            charge = float(metal_charges.get(el, PARAMETRIZED_METALS[el][4]))
        else:
            ad_type = el  # e.g. "Na", "K": Vina's non-AD metal typing
            charge = float(metal_charges.get(el, 1.0))
            warn(f"unparametrized metal {el} kept with standard Vina treatment "
                 f"({m.resname}{m.resseq}); supported ions: "
                 + ", ".join(sorted(PARAMETRIZED_METALS)))
        lines.append(pdbqt_line(m, ad_type, charge, serial))
        serial += 1
    return lines


# ---------------------------------------------------------------------------
# Meeko backend
# ---------------------------------------------------------------------------

MEEKO_TOOL = "mk_prepare_receptor.py"


def meeko_convert(protein_pdb_path, out_pdbqt, extra_args):
    """Run Meeko's mk_prepare_receptor on the protein-only PDB."""
    exe = shutil.which(MEEKO_TOOL)
    if exe is None:
        raise RuntimeError(
            f"'{MEEKO_TOOL}' not found on PATH. Install Meeko "
            "(https://github.com/forlilab/Meeko) or use --backend simple.")
    cmd = [exe, "--pdb", protein_pdb_path, "--out", out_pdbqt]
    if extra_args:
        cmd.extend(extra_args)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            "Meeko failed (exit {}):\n{}\nHint: common causes are unknown "
            "HETATM records (try --remove-het), missing atoms or unresolved "
            "altlocs; see --help for the preprocessing flags, or use "
            "--backend simple for a dependency-free conversion.".format(
                proc.returncode, (proc.stdout + proc.stderr)[-1500:]))
    if not os.path.exists(out_pdbqt):
        # some Meeko versions write '<out>.pdbqt' or into a directory
        alt = out_pdbqt + ".pdbqt"
        if os.path.exists(alt):
            return alt
        raise RuntimeError("Meeko did not produce the expected output file.")
    return out_pdbqt


def write_protein_pdb(path, atoms):
    """Write a minimal but column-correct PDB for the Meeko backend."""
    with open(path, "w") as fh:
        serial = 1
        for a in atoms:
            name = a.name if len(a.name) >= 4 else " " + a.name
            chain = a.chain if a.chain else " "
            try:
                resseq = int(a.resseq)
            except ValueError:
                resseq = 0
            icode = a.icode if a.icode else " "
            record = "HETATM" if a.record == "HETATM" else "ATOM  "
            fh.write(
                f"{record}{serial:5d} {name:<4s}{' ':1s}{a.resname:>3s} {chain:1s}"
                f"{resseq:4d}{icode:1s}   "
                f"{a.x:8.3f}{a.y:8.3f}{a.z:8.3f}{a.occ:6.2f}{a.bfac:6.2f}"
                f"{'':>10s}{a.element:>2s}\n")
            serial += 1
        fh.write("END\n")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def parse_metal_charges(spec):
    charges = {}
    if not spec:
        return charges
    for pair in spec.split(","):
        pair = pair.strip()
        if not pair:
            continue
        if "=" not in pair:
            raise SystemExit(f"ERROR: --metal-charge expects El=charge pairs, got '{pair}'")
        el, val = pair.split("=", 1)
        el = el.strip().capitalize()
        try:
            charges[el] = float(val)
        except ValueError:
            raise SystemExit(f"ERROR: --metal-charge: '{val}' is not a number (pair '{pair}')")
    return charges


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Vina-Ion-Metal receptor preparation (audits #5/#11): "
                    "robust PDB preprocessing, metal typing incl. Cu/Ni, and a "
                    "coordination pre-flight report.")
    ap.add_argument("--pdb", required=True, help="input receptor PDB")
    ap.add_argument("--out", required=True, help="output PDBQT path")
    ap.add_argument("--backend", choices=("auto", "meeko", "simple"), default="auto",
                    help="conversion backend: auto (meeko if available, else simple), "
                         "meeko (explicit hydrogens + charges; recommended), or "
                         "simple (no dependencies; vina/vinardo scoring only)")
    ap.add_argument("--altloc", default="first",
                    help="alternate location policy: 'first' (default) or a specific "
                         "altloc character to keep")
    ap.add_argument("--keep-waters", action="store_true",
                    help="keep water molecules (off by default)")
    ap.add_argument("--remove-het", action="store_true",
                    help="remove HETATM records other than metals and waters")
    ap.add_argument("--metal-charge", default="",
                    help="formal charges for metal records, e.g. 'Cu=2,Fe=3' "
                         "(defaults: parametrized ions +2)")
    ap.add_argument("--meeko-args", default="",
                    help="extra arguments forwarded to mk_prepare_receptor.py "
                         "(meeko backend only)")
    ap.add_argument("--version", action="version",
                    version=f"Vina-Ion-Metal prepare_receptor {VERSION}")
    args = ap.parse_args(argv)

    if not os.path.exists(args.pdb):
        raise SystemExit(f"ERROR: input PDB not found: {args.pdb}")

    metal_charges = parse_metal_charges(args.metal_charge)

    protein, het, waters, metals, problems = parse_pdb(
        args.pdb, altloc_mode=args.altloc, keep_waters=args.keep_waters,
        remove_het=args.remove_het)

    if not protein and not het:
        raise SystemExit("ERROR: no protein/HETATM atoms parsed from the input PDB.")

    for p in problems:
        warn(p)

    integrity = check_chain_integrity(protein)
    for p in integrity:
        warn(p)

    info(f"parsed {len(protein)} protein atoms, {len(het)} HETATM atoms, "
         f"{len(waters)} waters, {len(metals)} metal ions")

    # --- coordination pre-flight (before conversion: catches prep problems) ---
    if metals:
        coordination_report(metals, protein, het, waters, metal_charges)

    # --- conversion ---
    backend = args.backend
    if backend == "auto":
        backend = "meeko" if shutil.which(MEEKO_TOOL) else "simple"
        info(f"backend auto -> using '{backend}'")

    if backend == "meeko":
        with tempfile.TemporaryDirectory() as tmp:
            protein_pdb = os.path.join(tmp, "protein_only.pdb")
            # Meeko handles the protein (and kept HETATM/waters); metals are
            # stripped because Meeko cannot type Cu/Ni (audit #5) and are
            # re-appended below with the correct AD types.
            write_protein_pdb(protein_pdb, protein + het + waters)
            meeko_out = os.path.join(tmp, "meeko_out.pdbqt")
            try:
                produced = meeko_convert(protein_pdb, meeko_out,
                                         args.meeko_args.split() if args.meeko_args else [])
            except RuntimeError as exc:
                raise SystemExit(f"ERROR: {exc}")
            with open(produced) as fh:
                body = fh.read()
        extra = []
        if metals:
            # continue atom numbering after the Meeko output
            n_atoms = sum(1 for ln in body.splitlines()
                          if ln.startswith(("ATOM", "HETATM")))
            extra = metal_pdbqt_lines(metals, metal_charges, n_atoms + 1)
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w") as fh:
            fh.write(body)
            if not body.endswith("\n"):
                fh.write("\n")
            for ln in extra:
                fh.write(ln + "\n")
    else:
        body, _ = simple_convert(protein, het, waters, metals, metal_charges)
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w") as fh:
            fh.write(body)

    summary = []
    if metals:
        param = sorted({m.element for m in metals if m.element in PARAMETRIZED_METALS})
        unparam = sorted({m.element for m in metals if m.element not in PARAMETRIZED_METALS})
        if param:
            summary.append("coordination-scoring ions: " + ", ".join(param))
        if unparam:
            summary.append("unparametrized metals (standard Vina): " + ", ".join(unparam))
    info(f"wrote {args.out} (backend: {backend}"
         + ("; " + "; ".join(summary) if summary else "") + ")")
    if backend == "simple":
        info("NOTE: simple backend charges are zero -> use only with the "
             "vina/vinardo scoring functions (not AD4 maps).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
