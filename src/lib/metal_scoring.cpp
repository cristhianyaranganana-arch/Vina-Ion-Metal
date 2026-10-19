/*

  VinaMetal - metal_scoring module

  Copyright (c) 2024, VinaMetal developers
  Based on AutoDock Vina (Copyright (c) 2006-2010, The Scripps Research Institute)
  Licensed under the Apache License, Version 2.0.

*/

#include "metal_scoring.h"
#include "int_pow.h"
#include <cmath>
#include <algorithm>
#include <ostream>

const fl metal_vdw_cutoff = 8.0;    // Angstrom; matches Vina pair cutoff
const fl metal_pair_cap   = 100000; // AD4-style cap (potentials.h ad4_vdw)
const fl metal_cn_sigma   = 0.15;   // Angstrom; smoothness of CN sigmoid
const fl metal_ang_sigma  = 0.12;   // cosine units; smoothness of angular bonus

// Coordination donor set (XS typing, native Vina/Vinardo scoring).
// AD4 reference set is NA/NS/OA/OS/SA/S/N - see MAPA_DE_PORTERIA.md section 4.
// HD is excluded pending the full AD4Zn directional port (TODO).
bool metal_scoring::is_coordination_donor_type(sz xs) {
    switch (xs) {
        case XS_TYPE_N_A:
        case XS_TYPE_N_DA:
        case XS_TYPE_N_D:
        case XS_TYPE_O_A:
        case XS_TYPE_O_DA:
        case XS_TYPE_O_D:
        case XS_TYPE_S_P:
            return true;
        default:
            return false;
    }
}

bool metal_scoring::is_coordination_donor(const atom& a) {
    return is_coordination_donor_type(a.xs);
}

bool metal_scoring::is_metal_ion(const atom& a) const {
    if (a.el != EL_TYPE_Met) return false;
    return (m_params->index_for_ad_type(a.ad) != m_params->size());
}

bool metal_scoring::owns_pair(const atom& movable, const atom& rigid) const {
    if (!m_enabled) return false;
    if (rigid.el != EL_TYPE_Met) return false;
    if (!is_coordination_donor(movable)) return false;
    return is_metal_ion(rigid);
}

fl metal_scoring::donor_epsii(const atom& donor) const {
    if (donor.ad < AD_TYPE_SIZE) return ad_type_property(donor.ad).depth;
    return 0; // unknown type: no metal interaction
}

fl metal_scoring::donor_target_distance(const metal_ion& ion, const atom& donor) const {
    // Harding mean metal-donor distances by donor element
    switch (donor.el) {
        case EL_TYPE_N: return ion.d_N;
        case EL_TYPE_O: return ion.d_O;
        case EL_TYPE_S: return ion.d_S;
        default:        return ion.d_O; // fallback (should not happen for the donor set)
    }
}

// Smooth coordination count: sigmoid of the metal-donor distance.
// s(r) -> 1 for r << r_mid, s(r) -> 0 for r >> r_mid. C-infinity everywhere.
inline fl cn_sigmoid(fl r, fl r_mid) {
    return 1.0 / (1.0 + std::exp((r - r_mid) / metal_cn_sigma));
}
inline fl cn_sigmoid_derivative(fl r, fl r_mid, fl s) {
    return -s * (1.0 - s) / metal_cn_sigma;
}

// AD4 12-6 pair potential with the minimum placed at the coordination
// distance d0 (Harding) and depth epsij = sqrt(eps_M * eps_L) (AD4 rules).
// E(r) = epsij * ( (d0/r)^12 - 2 (d0/r)^6 )
inline fl metal_pair_energy(fl r, fl d0, fl epsij) {
    fl x6 = int_pow<6>(d0 / r);
    fl x12 = x6 * x6;
    fl e = epsij * (x12 - 2.0 * x6);
    return (std::min)(e, metal_pair_cap);
}
inline fl metal_pair_energy_deriv(fl r, fl d0, fl epsij) {
    // dE/dr = (12 * epsij / r) * ( (d0/r)^6 - (d0/r)^12 )
    fl x6 = int_pow<6>(d0 / r);
    fl x12 = x6 * x6;
    return (12.0 * epsij / r) * (x6 - x12);
}

fl metal_scoring::eval(const model& m, vecv* forces) const {
    if (!m_enabled) return 0;

    // Structural caches (metal/donor indices) depend only on atom topology.
    // They are normally prepared eagerly via prepare(); if a model with the
    // same topology (atom counts) is passed - e.g. a per-thread copy made by
    // parallel_mc - the cache stays valid and read-only (thread-safe). A model
    // with different topology triggers a rebuild (single-threaded paths only).
    if (!m_cache.prepared
        || m_cache.n_grid_atoms != m.grid_atoms.size()
        || m_cache.n_movable != m.num_movable_atoms()) {
        build_cache(m);
    }
    if (m_cache.metal_atoms.empty() || m_cache.donors.empty()) return 0;

    const szv& metal_atoms = m_cache.metal_atoms;
    const szv& metal_ions = m_cache.metal_ions;
    const szv& donors = m_cache.donors;

    fl e_total = 0;

    VINA_FOR_IN(k, metal_atoms) {
        const atom& m_atom = m.grid_atoms[metal_atoms[k]];
        const metal_ion& ion = m_params->ion(metal_ions[k]);
        const vec& mc = m_atom.coords;

        // -------- term (a): metal-specific vdW pair potential --------
        VINA_FOR_IN(d, donors) {
            const sz i = donors[d];
            const vec& pc = m.coords[i];
            vec diff; diff = pc - mc;
            fl r2 = sqr(diff);
            if (r2 >= sqr(metal_vdw_cutoff)) continue;
            const fl r = std::sqrt(r2);
            if (r < epsilon_fl) continue;
            const fl d0 = donor_target_distance(ion, m.atoms[i]);
            const fl epsij = std::sqrt(ion.epsii * donor_epsii(m.atoms[i]));
            const fl e = metal_pair_energy(r, d0, epsij);
            e_total += m_w_vdw * e;
            if (forces) {
                const fl dor = metal_pair_energy_deriv(r, d0, epsij);
                (*forces)[i] += (m_w_vdw * dor / r) * diff;
            }
        }

        // -------- term (b): coordination-number penalty --------
        // donors within the coordination shell
        szv shell; // positions within 'donors'
        flv shell_r;
        fl cn = 0;
        fl r_mid = ion.cn_cutoff - 1.5 * metal_cn_sigma;
        VINA_FOR_IN(d, donors) {
            const sz i = donors[d];
            vec diff; diff = m.coords[i] - mc;
            fl r2 = sqr(diff);
            if (r2 >= sqr(ion.cn_cutoff)) continue;
            fl r = std::sqrt(r2);
            if (r < epsilon_fl) continue;
            shell.push_back(i);
            shell_r.push_back(r);
            cn += cn_sigmoid(r, r_mid);
        }
        // The rigid receptor may already occupy part of the metal's
        // coordination sphere. Only the remaining sites are available to the
        // ligand, so the CN target is reduced accordingly (computed once in
        // build_cache; the receptor geometry is fixed).
        const fl cn_target = m_cache.target_cn[k];
        fl cn_dev = cn - cn_target;
        e_total += m_w_cn * sqr(cn_dev);
        if (forces) {
            VINA_FOR_IN(s, shell) {
                const sz i = shell[s];
                const fl r = shell_r[s];
                const fl s_val = cn_sigmoid(r, r_mid);
                const fl ds_dr = cn_sigmoid_derivative(r, r_mid, s_val);
                const fl dE_dr = 2.0 * m_w_cn * cn_dev * ds_dr;
                vec diff; diff = m.coords[i] - mc;
                (*forces)[i] += (dE_dr / r) * diff;
            }
        }

        // -------- term (c): coordination geometry bonus --------
        // bounded bonus for donor pairs matching the ideal inter-donor angles
        flv ideal_cos;
        {
            flv ideal_angles;
            metal_geometry_ideal_angles(ion.geometry, ideal_angles);
            VINA_FOR_IN(a_i, ideal_angles)
                ideal_cos.push_back(std::cos(ideal_angles[a_i] * pi / 180.0));
        }
        if (m_w_geom > 0 && shell.size() > 1) {
            VINA_FOR(s1, shell.size()) {
                VINA_RANGE(s2, s1 + 1, shell.size()) {
                    const sz i = shell[s1];
                    const sz j = shell[s2];
                    vec u; u = m.coords[i] - mc;
                    vec v; v = m.coords[j] - mc;
                    const fl ru = u.norm();
                    const fl rv = v.norm();
                    if (ru < epsilon_fl || rv < epsilon_fl) continue;
                    const fl cos_theta = (u * v) / (ru * rv);
                    // product-form aggregation of ideal angles: C-infinity,
                    // -> 1 when ANY ideal angle is matched, -> 0 otherwise
                    fl prod = 1.0;
                    VINA_FOR_IN(c_i, ideal_cos) {
                        fl dc = (cos_theta - ideal_cos[c_i]) / metal_ang_sigma;
                        prod *= (1.0 - std::exp(-sqr(dc)));
                    }
                    const fl f = 1.0 - prod;
                    // distance weighting: only coordinated donors count
                    const fl s_i = cn_sigmoid(shell_r[s1], r_mid);
                    const fl s_j = cn_sigmoid(shell_r[s2], r_mid);
                    const fl bonus = m_w_geom * ion.geom_bonus_max * s_i * s_j * f;
                    e_total -= bonus;
                    if (forces) {
                        // d(bonus)/d(cos_theta):
                        //   f = 1 - prod_k (1 - g_k),  g_k = exp(-((cos-c_k)/sigma)^2)
                        fl df_dcos = 0;
                        VINA_FOR_IN(c_i, ideal_cos) {
                            fl dc = (cos_theta - ideal_cos[c_i]) / metal_ang_sigma;
                            fl g = std::exp(-sqr(dc));
                            fl dg_dcos = -2.0 * dc * g / metal_ang_sigma;
                            fl prod_except = 1.0;
                            VINA_FOR_IN(c_j, ideal_cos) {
                                if (c_j == c_i) continue;
                                fl dc2 = (cos_theta - ideal_cos[c_j]) / metal_ang_sigma;
                                prod_except *= (1.0 - std::exp(-sqr(dc2)));
                            }
                            df_dcos += dg_dcos * prod_except;
                        }
                        // d(cos)/d(u) = v/(ru*rv) - cos * u/ru^2
                        vec dcos_du; dcos_du = (1.0 / (ru * rv)) * v - (cos_theta / (ru * ru)) * u;
                        vec dcos_dv; dcos_dv = (1.0 / (ru * rv)) * u - (cos_theta / (rv * rv)) * v;
                        const fl dE_dcos = m_w_geom * ion.geom_bonus_max * s_i * s_j * df_dcos; // d(-bonus)/d(cos)
                        // chain rule through the sigmoid distance weights:
                        // d s_i / d r_i < 0, d r_i / d coords_i = u/ru
                        const fl ds_i_dr = cn_sigmoid_derivative(shell_r[s1], r_mid, s_i);
                        const fl ds_j_dr = cn_sigmoid_derivative(shell_r[s2], r_mid, s_j);
                        // d(-bonus)/d(coords_i) = -W * [ s_j*f*ds_i_dr*u/ru + s_i*s_j*df_dcos*dcos_du ]
                        vec grad_i; grad_i = -(m_w_geom * ion.geom_bonus_max) *
                            ( s_j * f * ds_i_dr / ru * u + s_i * s_j * df_dcos * dcos_du );
                        vec grad_j; grad_j = -(m_w_geom * ion.geom_bonus_max) *
                            ( s_i * f * ds_j_dr / rv * v + s_i * s_j * df_dcos * dcos_dv );
                        // accumulate into minus_forces using the Vina convention
                        // minus_forces[x] += dE/d(coords_x)
                        (*forces)[i] += grad_i;
                        (*forces)[j] += grad_j;
                    }
                }
            }
        }
    }
    return e_total;
}

void metal_scoring::prepare(const model& m) {
    build_cache(m);
}

void metal_scoring::build_cache(const model& m) const {
    m_cache.prepared = true;
    m_cache.key = (const void*)&m;
    m_cache.n_grid_atoms = m.grid_atoms.size();
    m_cache.n_movable = m.num_movable_atoms();
    m_cache.metal_atoms.clear();
    m_cache.metal_ions.clear();
    m_cache.donors.clear();
    m_cache.occupied_cn.clear();
    m_cache.target_cn.clear();

    VINA_FOR(i, m.grid_atoms.size()) {
        const atom& a = m.grid_atoms[i];
        if (a.el != EL_TYPE_Met) continue;
        sz idx = m_params->index_for_ad_type(a.ad);
        if (idx == m_params->size()) continue; // unparametrized metal: standard Vina behavior
        m_cache.metal_atoms.push_back(i);
        m_cache.metal_ions.push_back(idx);

        // Count rigid-receptor donors already inside the coordination shell
        // (protein side chains, cofactors, retained waters...). Their sites are
        // not available to the ligand, so the effective ligand CN target is
        // reduced by this occupancy.
        const metal_ion& ion = m_params->ion(idx);
        int occupied = 0;
        VINA_FOR(j, m.grid_atoms.size()) {
            if (j == i) continue;
            const atom& donor = m.grid_atoms[j];
            if (!is_coordination_donor(donor)) continue;
            vec diff; diff = donor.coords - a.coords;
            if (sqr(diff) < sqr(ion.cn_cutoff)) ++occupied;
        }
        int target = ion.ideal_cn - occupied;
        if (target < 0) target = 0;
        m_cache.occupied_cn.push_back(occupied);
        m_cache.target_cn.push_back(target);
    }

    VINA_FOR(i, m.num_movable_atoms()) {
        const atom& a = m.atoms[i];
        if (is_coordination_donor(a) && a.ad < AD_TYPE_SIZE)
            m_cache.donors.push_back(i);
    }
}

void metal_scoring::report_cn_targets(std::ostream& out) const {
    if (!m_enabled || !m_cache.prepared) return;
    VINA_FOR_IN(k, m_cache.metal_atoms) {
        const metal_ion& ion = m_params->ion(m_cache.metal_ions[k]);
        out << "VinaMetal: " << ion.name << " CN target = ideal " << ion.ideal_cn
            << " - receptor-bound " << m_cache.occupied_cn[k]
            << " = " << m_cache.target_cn[k] << "\n";
    }
}

void metal_scoring::eval_breakdown(const model& m, fl& e_vdw, fl& e_cn, fl& e_geom, int& n_contacts) const {
    // reports by evaluating with one weight active at a time (no forces)
    e_vdw = 0; e_cn = 0; e_geom = 0; n_contacts = 0;
    if (!m_enabled) return;
    metal_scoring vdw_only(m_params, 1.0, 0, 0, m_cn_strict);
    metal_scoring cn_only(m_params, 0, 1.0, 0, m_cn_strict);
    metal_scoring geom_only(m_params, 0, 0, 1.0, m_cn_strict);
    e_vdw = vdw_only.eval(m, NULL);
    e_cn = cn_only.eval(m, NULL);
    e_geom = geom_only.eval(m, NULL);
    // count donor-metal contacts within the coordination shells
    VINA_FOR(i, m.grid_atoms.size()) {
        const atom& a = m.grid_atoms[i];
        if (a.el != EL_TYPE_Met) continue;
        sz idx = m_params->index_for_ad_type(a.ad);
        if (idx == m_params->size()) continue;
        const metal_ion& ion = m_params->ion(idx);
        VINA_FOR(j, m.num_movable_atoms()) {
            if (!is_coordination_donor(m.atoms[j])) continue;
            vec diff; diff = m.coords[j] - a.coords;
            if (sqr(diff) < sqr(ion.cn_cutoff)) ++n_contacts;
        }
    }
}
