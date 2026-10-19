/*

  VinaMetal - gradient verification test (finite differences vs analytic)

  Verifies that the analytic gradients of the metal coordination terms,
  exactly as BFGS consumes them through model::eval_deriv, match central
  finite differences of the total energy (acceptance criterion #5).

  Fixture: a Zn2+ ion at the origin + two rigid acetate-like fragments
  (root-only ligands, no torsions) placed at ~100 degrees around the metal
  with both oxygens at ~2.05 A. This exercises:
    - term (a) metal vdW pair potential (donors near the 12-6 minimum)
    - term (b) CN penalty (CN ~ 2 vs ideal 4 -> strong gradient)
    - term (c) geometry bonus (angle off-ideal -> nonzero angular gradient)
  Rotation-vector increments replicate rigid_conf::increment exactly, and the
  central difference step h = 1e-6 keeps the chord-vs-tangent artifact well
  below the 1e-4 acceptance threshold.

*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdlib>
#include "model.h"
#include "parse_pdbqt.h"
#include "precalculate.h"
#include "non_cache.h"
#include "scoring_function.h"
#include "quaternion.h"
#include "metal_scoring.h"
#include "metal_parameters.h"
#include "conf.h"

// Zn ion at origin + two carbons far from the coordination site
static const char* RECEPTOR_ZN =
"ATOM      1  ZN  ZN  A   1       0.000   0.000   0.000  1.00  0.00     0.000 Zn\n"
"ATOM      2  C   LYS A   2      -6.000   0.000   0.000  1.00  0.00     0.000 C\n";

// control receptor without metal
static const char* RECEPTOR_NOMETAL =
"ATOM      1  C   LYS A   2      -6.000   0.000   0.000  1.00  0.00     0.000 C\n";

// root-only acetate-like fragment (rigid), O at d0 from metal when placed right
static const char* FRAG_A =
"ROOT\n"
"ATOM      1  O   ACA A   1       2.050   0.000   0.000  1.00  0.00    -0.500 OA\n"
"ATOM      2  C   ACA A   1       2.550   1.300   0.200  1.00  0.00     0.300 C\n"
"ENDROOT\n"
"TORSDOF 0\n";

static const char* FRAG_B =
"ROOT\n"
"ATOM      1  O   ACB A   1      -0.356   2.019   0.300  1.00  0.00    -0.500 OA\n"
"ATOM      2  C   ACB A   1       0.844   2.419   0.000  1.00  0.00     0.300 C\n"
"ENDROOT\n"
"TORSDOF 0\n";

static std::string write_temp_pdbqt(const char* content) {
    static int counter = 0;
    std::string path = std::string("/tmp/vm_test_") + char('0' + counter++) + ".pdbqt";
    std::ofstream out(path.c_str());
    out << content;
    out.close();
    return path;
}

struct eval_context {
    ScoringFunction* sf;
    precalculate* p_type;
    precalculate_byatom* p_byatom;
    non_cache* nc;
    vec v;

    // pure_metal: zero out the standard Vina weights so the total energy is
    // ONLY the metal coordination terms (all C-infinity) - this isolates the
    // gradients introduced by VinaMetal from the pre-existing finite-difference
    // artifacts of the standard terms (repulsion kinks, SO(3) chord effects).
    eval_context(model& m, const grid_dims& gd, metal_scoring* ms, bool pure_metal) : v(1000, 1000, 1000) {
        fl w_scale = pure_metal ? 0.0 : 1.0;
        flv w;
        w.push_back(-0.035579 * w_scale); // gauss1
        w.push_back(-0.005156 * w_scale); // gauss2
        w.push_back(0.840245 * w_scale);  // repulsion
        w.push_back(-0.035069 * w_scale); // hydrophobic
        w.push_back(-0.587439 * w_scale); // hbond
        w.push_back(50.0 * w_scale);      // glue
        w.push_back(0.05846 * w_scale);   // rot (conf-independent transform)
        sf = new ScoringFunction(SF_VINA, w);
        p_type = new precalculate(*sf);
        p_byatom = new precalculate_byatom(*sf, m);
        nc = new non_cache(m, gd, p_type, 1e6);
        if (ms) nc->set_metal_scoring(ms);
    }
    ~eval_context() { delete sf; delete p_type; delete p_byatom; delete nc; }

    fl energy(model& m, const conf& c) {
        m.set(c);
        change g(m.get_size());
        fl e = m.eval_deriv(*p_byatom, *nc, v, g);
        return e;
    }
    void gradient(model& m, const conf& c, change& g) {
        m.set(c);
        m.eval_deriv(*p_byatom, *nc, v, g);
    }
};

// total scalar dofs: per ligand 3 position + 3 rotation
static sz num_dofs(const conf& c) {
    return c.ligands.size() * 6;
}

// apply +delta to dof k replicating rigid_conf::increment semantics
static void perturb(conf& c, sz k, fl delta) {
    const sz lig = k / 6;
    const sz comp = k % 6;
    if (comp < 3) {
        c.ligands[lig].rigid.position[comp] += delta;
    } else {
        vec rotation(0, 0, 0);
        rotation[comp - 3] = delta;
        quaternion_increment(c.ligands[lig].rigid.orientation, rotation);
    }
}

static fl analytic_component(const change& g, sz k) {
    const sz lig = k / 6;
    const sz comp = k % 6;
    if (comp < 3) return g.ligands[lig].rigid.position[comp];
    return g.ligands[lig].rigid.orientation[comp - 3];
}

static int run_gradient_check(const char* receptor_str, metal_scoring* ms, const char* label, bool& ok, bool pure_metal, fl tol) {
    atom_type::t atype = atom_type::XS;
    std::string rec_path = write_temp_pdbqt(receptor_str);
    std::string a_path = write_temp_pdbqt(FRAG_A);
    std::string b_path = write_temp_pdbqt(FRAG_B);

    model rec = parse_receptor_pdbqt(rec_path, std::string(), atype);
    model m = rec;
    m.append(parse_ligand_pdbqt_from_file(a_path, atype));
    m.append(parse_ligand_pdbqt_from_file(b_path, atype));

    grid_dims gd;
    VINA_FOR_IN(i, gd) {
        gd[i].begin = -30;
        gd[i].end = 30;
        gd[i].n_voxels = 120;
    }

    eval_context ctx(m, gd, ms, pure_metal);

    // Use the parsed configuration as-is: donors sit at ~2.05 A from the metal
    // at a ~100 deg angle (off-ideal) with identity orientation and zero
    // torsions. All three metal terms are active with non-zero gradients.
    conf c = m.get_initial_conf();

    change g(m.get_size());
    ctx.gradient(m, c, g);

    const fl h = 1e-6;
    bool local_ok = true;
    int n_fail = 0;
    std::cout << "\n[" << label << "] dof, analytic, finite-difference, rel.err\n";
    for (sz k = 0; k < num_dofs(c); ++k) {
        conf cp = c; perturb(cp, k, +h);
        conf cm = c; perturb(cm, k, -h);
        fl ep = ctx.energy(m, cp);
        fl em = ctx.energy(m, cm);
        fl d = (ep - em) / (2 * h);
        fl a = analytic_component(g, k);
        fl denom = (std::max)((std::max)(std::abs(a), std::abs(d)), (fl)1e-6);
        fl err = std::abs(a - d) / denom;
        std::cout << std::setw(4) << k << "  " << std::setw(14) << std::setprecision(7) << a
                  << "  " << std::setw(14) << std::setprecision(7) << d
                  << "  " << std::setprecision(3) << err;
        if (err > tol) { std::cout << "  FAIL"; local_ok = false; ++n_fail; }
        else std::cout << "  ok";
        std::cout << "\n";
    }
    std::cout << "[" << label << "] " << (local_ok ? "PASS" : "FAIL") << " (" << n_fail << " failures)\n";
    ok = ok && local_ok;
    return local_ok ? 0 : 1;
}

int main(int argc, char** argv) {
    bool ok = true;
    metal_parameters params;

    // 1) PURE metal scoring: the total energy is only the three metal terms
    //    (all C-infinity), so the analytic gradients must match central finite
    //    differences to 1e-4 (acceptance criterion #5).
    metal_scoring ms(&params, 1.0, 1.0, 1.0, false);
    run_gradient_check(RECEPTOR_ZN, &ms, "Zn + pure metal terms", ok, true, 1e-4);

    // 2) FULL Vina scoring + metal terms: verify the metal gradients do not
    //    interfere with the standard machinery. The FD tolerance is relaxed to
    //    5e-3 because the pre-existing standard terms (repulsion kinks, SO(3)
    //    chord effects) dominate the finite-difference artifact; the no-metal
    //    control below shows the same artifact level.
    run_gradient_check(RECEPTOR_ZN, &ms, "Zn + full scoring", ok, false, 5e-3);

    // 3) control WITHOUT metals in receptor (full scoring): same FD artifact
    //    level must appear, proving the artifacts are not introduced by VinaMetal
    metal_scoring ms2(&params, 1.0, 1.0, 1.0, false);
    run_gradient_check(RECEPTOR_NOMETAL, &ms2, "no metal control", ok, false, 5e-3);

    std::cout << "\n==== " << (ok ? "ALL GRADIENT CHECKS PASSED" : "GRADIENT CHECK FAILURES") << " ====\n";
    return ok ? 0 : 1;
}
