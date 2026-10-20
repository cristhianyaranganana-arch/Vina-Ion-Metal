/*

  VinaMetal - metal_parameters module

  Copyright (c) 2024, VinaMetal developers
  Based on AutoDock Vina (Copyright (c) 2006-2010, The Scripps Research Institute)
  Licensed under the Apache License, Version 2.0.

  Sources of the built-in values (DO NOT INVENT - rule of the porting):
   - Rii, epsii, vol, solpar : AutoDock4 repository, file AD4.1_bound.dat
       Mg  line 98 | Ca line 105 | Mn line 107 | Fe line 109 | Zn line 111
       donors: N/NA/NS lines 91-95 (Rii 3.50, eps 0.160), O/OA/OS lines 96-97
       (Rii 3.20, eps 0.200), S/SA lines 102-103 (Rii 4.00, eps 0.200)
   - ideal CN, geometry, metal-donor distances : Harding,
       "Geometry of metal-ligand interactions in proteins", Acta Cryst. D 57
       (2001) 401-411 and Acta Cryst. D 62 (2006) 678-682.
   - Cu2+/Ni2+ are NOT calibrated in AD4 (see Scripps "how to add new atom
       types"): values below are PLACEHOLDER-UNVALIDATED estimates and must be
       refined or overridden through --metal_params.

*/

#include "metal_parameters.h"
#include <boost/optional.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <cmath>
#include <fstream>

metal_parameters::metal_parameters() {
    add_defaults();
}

void metal_parameters::add_defaults() {
    ions.clear();
    metal_ion m;

    // ---------------- Zn2+ (AD4.1_bound.dat line 111) ----------------
    m = metal_ion();
    m.name = "Zn2+"; m.element = "Zn"; m.ad_type = AD_TYPE_Zn;
    m.Rii = 1.48; m.epsii = 0.550; m.vol = 1.70; m.solpar = -0.00110;
    m.formal_charge = 2;
    m.ideal_cn = 4; m.geometry = METAL_GEOM_TETRAHEDRAL;       // Harding 2001: tetrahedral dominant
    m.d_N = 2.05; m.d_O = 2.05; m.d_S = 2.32;                  // Harding 2001 Table 1
    m.cn_cutoff = 3.0;
    m.geom_bonus_max = 1.0;                                    // PLACEHOLDER: calibrated in benchmark phase
    ions.push_back(m);

    // ---------------- Mg2+ (AD4.1_bound.dat line 98) -----------------
    m = metal_ion();
    m.name = "Mg2+"; m.element = "Mg"; m.ad_type = AD_TYPE_Mg;
    m.Rii = 1.30; m.epsii = 0.875; m.vol = 1.56; m.solpar = -0.00110;
    m.formal_charge = 2;
    m.ideal_cn = 6; m.geometry = METAL_GEOM_OCTAHEDRAL;        // Harding 2001: octahedral dominant
    m.d_N = 2.10; m.d_O = 2.08; m.d_S = 2.50;
    m.cn_cutoff = 3.0;
    m.geom_bonus_max = 1.0;
    ions.push_back(m);

    // ---------------- Ca2+ (AD4.1_bound.dat line 105) ----------------
    m = metal_ion();
    m.name = "Ca2+"; m.element = "Ca"; m.ad_type = AD_TYPE_Ca;
    m.Rii = 1.98; m.epsii = 0.550; m.vol = 2.77; m.solpar = -0.00110;
    m.formal_charge = 2;
    m.ideal_cn = 6; m.geometry = METAL_GEOM_OCTAHEDRAL;        // Harding 2001: CN 6-8, we use 6
    m.d_N = 2.45; m.d_O = 2.40; m.d_S = 2.80;
    m.cn_cutoff = 3.3;
    m.geom_bonus_max = 1.0;
    ions.push_back(m);

    // ---------------- Mn2+ (AD4.1_bound.dat line 107) ----------------
    m = metal_ion();
    m.name = "Mn2+"; m.element = "Mn"; m.ad_type = AD_TYPE_Mn;
    m.Rii = 1.30; m.epsii = 0.875; m.vol = 2.14; m.solpar = -0.00110;
    m.formal_charge = 2;
    m.ideal_cn = 6; m.geometry = METAL_GEOM_OCTAHEDRAL;
    m.d_N = 2.20; m.d_O = 2.15; m.d_S = 2.45;
    m.cn_cutoff = 3.1;
    m.geom_bonus_max = 1.0;
    ions.push_back(m);

    // ---------------- Fe2+ (AD4.1_bound.dat line 109, generic Fe) ----
    // NOTE: AD4 only defines generic "Fe". Oxidation-state-specific refinement
    // is a TODO; both entries currently use the AD4 generic values.
    m = metal_ion();
    m.name = "Fe2+"; m.element = "Fe"; m.ad_type = AD_TYPE_Fe;
    m.Rii = 1.30; m.epsii = 0.010; m.vol = 1.84; m.solpar = -0.00110;
    m.formal_charge = 2;
    m.ideal_cn = 6; m.geometry = METAL_GEOM_OCTAHEDRAL;
    m.d_N = 2.00; m.d_O = 2.00; m.d_S = 2.25;
    m.cn_cutoff = 3.0;
    m.geom_bonus_max = 1.0;
    ions.push_back(m);

    // ---------------- Fe3+ (AD4.1_bound.dat line 109, generic Fe) ----
    // TODO: oxidation-state-specific parameters pending literature validation.
    m = metal_ion();
    m.name = "Fe3+"; m.element = "Fe"; m.ad_type = AD_TYPE_Fe;
    m.Rii = 1.30; m.epsii = 0.010; m.vol = 1.84; m.solpar = -0.00110;
    m.formal_charge = 3;
    m.ideal_cn = 6; m.geometry = METAL_GEOM_OCTAHEDRAL;
    m.d_N = 2.00; m.d_O = 2.00; m.d_S = 2.25;
    m.cn_cutoff = 3.0;
    m.geom_bonus_max = 1.0;
    ions.push_back(m);

    // ---------------- Cu2+ (NOT in AD4.1_bound.dat) ------------------
    // PLACEHOLDER-UNVALIDATED: estimate consistent with AD4 metal scales.
    // Must be validated or overridden via --metal_params (porting rule).
    m = metal_ion();
    m.name = "Cu2+"; m.element = "Cu"; m.ad_type = AD_TYPE_Cu;
    m.Rii = 1.40; m.epsii = 0.550; m.vol = 1.70; m.solpar = -0.00110;
    m.formal_charge = 2;
    m.ideal_cn = 4; m.geometry = METAL_GEOM_SQUARE_PLANAR;     // Jahn-Teller distorted Cu2+ (also see Vina issue #368 discussion)
    m.d_N = 2.00; m.d_O = 1.98; m.d_S = 2.22;
    m.cn_cutoff = 2.9;
    m.geom_bonus_max = 1.0;
    ions.push_back(m);

    // ---------------- Ni2+ (NOT in AD4.1_bound.dat) ------------------
    // PLACEHOLDER-UNVALIDATED: estimate consistent with AD4 metal scales.
    m = metal_ion();
    m.name = "Ni2+"; m.element = "Ni"; m.ad_type = AD_TYPE_Ni;
    m.Rii = 1.36; m.epsii = 0.550; m.vol = 1.70; m.solpar = -0.00110;
    m.formal_charge = 2;
    m.ideal_cn = 6; m.geometry = METAL_GEOM_OCTAHEDRAL;        // Harding 2006: six-coordinate dominant
    m.d_N = 2.05; m.d_O = 2.05; m.d_S = 2.25;
    m.cn_cutoff = 3.0;
    m.geom_bonus_max = 1.0;
    ions.push_back(m);
}

sz metal_parameters::index_for_ad_type(sz ad_type) const {
    // Pick the entry matching ad_type, honoring the user-selected oxidation
    // state for elements with more than one entry (default: +2).
    sz best = ions.size();
    VINA_FOR(i, ions.size()) {
        const metal_ion& m = ions[i];
        if (m.ad_type != ad_type) continue;
        if (has_oxidation_choice(m.element)) {
            int state = get_oxidation_state(m.element);
            if (state == 0) state = 2; // default: +2
            if (m.formal_charge == state) return i;
        } else {
            return i;
        }
        if (best == ions.size()) best = i; // fallback: first match
    }
    return best;
}

bool metal_parameters::has_oxidation_choice(const std::string& element) const {
    sz count = 0;
    VINA_FOR(i, ions.size())
        if (ions[i].element == element) ++count;
    return count > 1;
}

void metal_parameters::load_json(const std::string& path) {
    try {
        boost::property_tree::ptree pt;
        std::ifstream in(path.c_str());
        if (!in) {
            std::cerr << "ERROR: cannot open metal parameters file '" << path << "'\n";
            exit(EXIT_FAILURE);
        }
        boost::property_tree::read_json(in, pt);
        VINA_FOR_IN(i, ions) {
            metal_ion& m = ions[i];
            boost::optional<boost::property_tree::ptree&> node_opt =
                pt.get_child_optional(m.name);
            if (!node_opt) continue;
            const boost::property_tree::ptree& node = *node_opt;
            m.Rii     = node.get<fl>("Rii", m.Rii);
            m.epsii   = node.get<fl>("epsii", m.epsii);
            m.vol     = node.get<fl>("vol", m.vol);
            m.solpar  = node.get<fl>("solpar", m.solpar);
            m.ideal_cn = node.get<int>("ideal_cn", m.ideal_cn);
            std::string g = node.get<std::string>("geometry",
                metal_geometry_name(m.geometry));
            if (g == "tetrahedral")        m.geometry = METAL_GEOM_TETRAHEDRAL;
            else if (g == "octahedral")    m.geometry = METAL_GEOM_OCTAHEDRAL;
            else if (g == "square-planar") m.geometry = METAL_GEOM_SQUARE_PLANAR;
            m.d_N = node.get<fl>("d_N", m.d_N);
            m.d_O = node.get<fl>("d_O", m.d_O);
            m.d_S = node.get<fl>("d_S", m.d_S);
            m.cn_cutoff    = node.get<fl>("cn_cutoff", m.cn_cutoff);
            m.geom_bonus_max = node.get<fl>("geom_bonus_max", m.geom_bonus_max);
        }
    } catch (std::exception& e) {
        std::cerr << "ERROR: failed to parse metal parameters file '" << path
                  << "': " << e.what() << "\n";
        exit(EXIT_FAILURE);
    }
}
