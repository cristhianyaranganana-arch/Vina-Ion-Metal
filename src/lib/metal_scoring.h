/*

  VinaMetal - metal_scoring module

  Copyright (c) 2024, VinaMetal developers
  Based on AutoDock Vina (Copyright (c) 2006-2010, The Scripps Research Institute)
  Licensed under the Apache License, Version 2.0.

  Coordination-aware scoring terms, conditional on the presence of parametrized
  metal ions in the rigid receptor. Physics ported from AutoDock 4 with the
  coordination geometry model of Harding (Acta Cryst. D 2001/2006).

  Three terms (O3):
   (a) metal-specific vdW: AD4 12-6 pair potential with AD4 combination rules
       (Rij = sum of radii convention absorbed into a minimum placed at the
       Harding metal-donor distance; epsij = geometric mean of epsii), applied
       ONLY to coordination donor pairs; those pairs are excluded from the
       standard Vina potentials to avoid double counting.
   (b) coordination-number penalty: smooth differentiable count of donors in
       the coordination shell, penalizing deviation from the ideal CN.
   (c) coordination geometry bonus: bounded bonus for donor pairs whose
       M-L-L' angle matches the preferred geometry (tetrahedral 109.47 deg,
       octahedral / square-planar 90 and 180 deg).

  All terms are C-infinity smooth (no kinks) so BFGS can consume analytic
  gradients. Gradients are accumulated into model::minus_forces with the same
  convention as eval_interacting_pairs_deriv:
      forces[i] += (dE/dr / r) * (coords[i] - partner_coords)

  NOTE (grid vs direct): terms (b) and (c) are many-body (angle/counting) and
  cannot be precomputed into Vina's affinity maps; term (a) is distance-only
  but per-metal-ion (Zn vs Mg share XS_TYPE_Met_D), so all three terms are
  evaluated directly over movable donor x receptor metal pairs, following the
  Vina-Carb pattern for non-pairwise potentials.

*/

#ifndef VINA_METAL_SCORING_H
#define VINA_METAL_SCORING_H

#include <iosfwd>
#include "common.h"
#include "atom_constants.h"
#include "model.h"
#include "metal_parameters.h"

struct metal_scoring {
    metal_scoring() : m_enabled(false), m_params(NULL), m_w_vdw(0), m_w_cn(0), m_w_geom(0), m_cn_strict(false) {}
    metal_scoring(const metal_parameters* params, fl w_vdw, fl w_cn, fl w_geom, bool cn_strict)
        : m_enabled(true), m_params(params), m_w_vdw(w_vdw), m_w_cn(w_cn), m_w_geom(w_geom), m_cn_strict(cn_strict) {}

    bool enabled() const { return m_enabled; }
    fl weight_vdw() const { return m_w_vdw; }
    fl weight_cn()  const { return m_w_cn; }
    fl weight_geom()const { return m_w_geom; }

    // Total metal coordination energy of the current model configuration.
    // If forces != NULL, accumulates gradients into (*forces)[i] for movable
    // atoms (same convention as the rest of Vina).
    fl eval(const model& m, vecv* forces) const;

    // Precompute the structural caches (metal/donor indices) for the given
    // model. MUST be called (through Vina::refresh_metal_scoring paths) before
    // the parallel search: the cached indices are structural and remain valid
    // for per-thread copies of the same model, keeping eval() read-only and
    // thread-safe (O4 requirement).
    void prepare(const model& m);

    // Per-term breakdown for reporting (no forces).
    void eval_breakdown(const model& m, fl& e_vdw, fl& e_cn, fl& e_geom, int& n_contacts) const;

    // Human-readable summary of the receptor-occupancy adjustment applied to
    // each metal's CN target (requires prepare() to have been called).
    void report_cn_targets(std::ostream& out) const;

    // True if receptor atom a is a parametrized metal ion.
    bool is_metal_ion(const atom& a) const;
    // True if movable/receptor atom type can coordinate metals (donor set).
    static bool is_coordination_donor(const atom& a);
    static bool is_coordination_donor_type(sz xs); // XS-type variant (map population)
    // True if this movable donor x rigid atom pair must be handled by metal_scoring.
    bool owns_pair(const atom& movable, const atom& rigid) const;

private:
    bool m_enabled;
    const metal_parameters* m_params;
    fl m_w_vdw;
    fl m_w_cn;
    fl m_w_geom;
    bool m_cn_strict;

    // lazy structural caches; indices are topology-dependent only (atom types),
    // so they stay valid for copies of the same model used by worker threads
    struct model_cache {
        bool prepared;
        const void* key;       // model that built the cache (single-threaded paths)
        sz n_grid_atoms;
        sz n_movable;
        szv metal_atoms; // indices into m.grid_atoms
        szv metal_ions;  // indices into metal_parameters
        szv donors;      // indices into m.atoms (movable)
        std::vector<int> occupied_cn; // receptor donors already bound to each metal
        std::vector<int> target_cn;   // effective ligand CN target per metal
        model_cache() : prepared(false), key(NULL), n_grid_atoms(0), n_movable(0) {}
    };
    mutable model_cache m_cache;

    void build_cache(const model& m) const;

    fl donor_target_distance(const metal_ion& ion, const atom& donor) const; // d0
    fl donor_epsii(const atom& donor) const;
};

#endif
