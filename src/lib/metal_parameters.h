/*

  VinaMetal - metal_parameters module

  Copyright (c) 2024, VinaMetal developers
  Based on AutoDock Vina (Copyright (c) 2006-2010, The Scripps Research Institute)
  Licensed under the Apache License, Version 2.0.

  Metal ion parameter table for coordination-aware docking.
  Physics ported from AutoDock 4 (AD4.1_bound.dat) and Harding (Acta Cryst. D).

*/

#ifndef VINA_METAL_PARAMETERS_H
#define VINA_METAL_PARAMETERS_H

#include <string>
#include <vector>
#include <map>
#include "common.h"
#include "atom_constants.h"

// Preferred coordination geometries (O3c)
enum metal_geometry {
    METAL_GEOM_TETRAHEDRAL   = 0,
    METAL_GEOM_OCTAHEDRAL    = 1,
    METAL_GEOM_SQUARE_PLANAR = 2
};

// Ideal inter-ligand angles (degrees) per geometry
inline void metal_geometry_ideal_angles(metal_geometry g, flv& angles) {
    angles.clear();
    switch (g) {
        case METAL_GEOM_TETRAHEDRAL:
            angles.push_back(109.47);
            break;
        case METAL_GEOM_OCTAHEDRAL:
        case METAL_GEOM_SQUARE_PLANAR:
            angles.push_back(90.0);
            angles.push_back(180.0);
            break;
    }
}

inline std::string metal_geometry_name(metal_geometry g) {
    switch (g) {
        case METAL_GEOM_TETRAHEDRAL:   return "tetrahedral";
        case METAL_GEOM_OCTAHEDRAL:    return "octahedral";
        case METAL_GEOM_SQUARE_PLANAR: return "square-planar";
    }
    return "unknown";
}

struct metal_ion {
    std::string name;          // e.g. "Zn2+"
    std::string element;       // e.g. "Zn"
    sz ad_type;                // AD_TYPE_* of the receptor PDBQT atom (AD_TYPE_SIZE if unused)
    // --- AD4 atomic parameters (ported verbatim from AD4.1_bound.dat when available) ---
    fl Rii;                    // sum of vdW radii of two like atoms (Angstrom)
    fl epsii;                  // vdW well depth (kcal/mol)
    fl vol;                    // atomic solvation volume (Angstrom^3)
    fl solpar;                 // atomic solvation parameter
    // --- coordination model (Harding, Acta Cryst. D 2001/2006 + user-tunable) ---
    int formal_charge;
    int ideal_cn;              // ideal coordination number
    metal_geometry geometry;
    fl d_N;                    // mean metal-donor distance to N (Angstrom)
    fl d_O;                    // mean metal-donor distance to O (Angstrom)
    fl d_S;                    // mean metal-donor distance to S (Angstrom)
    fl cn_cutoff;              // coordination shell radius (Angstrom)
    fl geom_bonus_max;         // max geometric bonus per donor pair (kcal/mol), bounded
    // TODO(AD4Zn): port directional Zn parameters from AD4Zn.dat (nbp_r_eps) for Zn2+.

    fl radius() const { return Rii * 0.5; } // AD4 convention: Rii is the like-pair sum
};

struct metal_parameters {
    metal_parameters(); // loads the built-in table (AD4.1_bound.dat + Harding)

    void load_json(const std::string& path); // per-ion override, --metal_params

    sz size() const { return ions.size(); }
    const metal_ion& ion(sz i) const { return ions[i]; }

    // Oxidation-state selection for elements with multiple entries (Fe2+/Fe3+).
    // PDBQT does not encode oxidation states, so the user must declare them;
    // the default is the +2 state.
    void set_oxidation_state(const std::string& element, int state) { m_oxidation[element] = state; }
    int get_oxidation_state(const std::string& element) const {
        std::map<std::string, int>::const_iterator it = m_oxidation.find(element);
        return (it == m_oxidation.end()) ? 0 : it->second;
    }

    // Receptor atom (el == EL_TYPE_Met) -> ion index, or ions.size() if not parametrized.
    // ad_type is authoritative for Mg/Mn/Zn/Ca/Fe/Cu/Ni; the oxidation state map
    // selects among entries of the same element (Fe=2/3).
    sz index_for_ad_type(sz ad_type) const;
    bool has_oxidation_choice(const std::string& element) const; // elements with >1 entry

    std::vector<metal_ion> ions;

private:
    void add_defaults();
    std::map<std::string, int> m_oxidation;
};

#endif
