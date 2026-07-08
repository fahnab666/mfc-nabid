#!/usr/bin/env python3
# 2D JWL products blast with Euler-Lagrange particles.
# Single-material JWL (PETN products) so the closure runs in its most robust
# configuration; a high-pressure circular region expands into lower-pressure
# products and sweeps a cloud of Lagrangian particles. Sized to run on 8 CPU
# ranks in well under 15 minutes on macOS.
import json

# PETN JWL products parameters (Q in J/kg; MFC derives jwl_E0 = rho0*Q).
jwl_A = 6.17e11
jwl_B = 2.11e9
jwl_R1 = 4.40
jwl_R2 = 1.20
jwl_omega = 0.25
jwl_rho0 = 1770.0
jwl_Cv = 900.0
jwl_Q = 5.96e6

# --- domain / grid (2D) ---
xmin, xmax = -0.05, 0.05
ymin, ymax = -0.05, 0.05
m = 399
n = 399

data = {
    "run_time_info": "T",
    # domain
    "x_domain%beg": xmin,
    "x_domain%end": xmax,
    "y_domain%beg": ymin,
    "y_domain%end": ymax,
    "m": m,
    "n": n,
    "p": 0,
    "cyl_coord": "F",
    # time: adaptive CFL, short physical window
    "cfl_adap_dt": "T",
    "cfl_target": 0.3,
    "cfl_const_dt": "F",
    "n_start": 0,
    "t_stop": 1.0e-4,
    "t_save": 1.0e-5,
    # numerics
    "model_eqns": 2,
    "num_fluids": 1,
    "num_patches": 2,
    "alt_soundspeed": "F",
    "mixture_err": "F",
    "mpp_lim": "F",
    "time_stepper": 3,
    "weno_order": 5,
    "mapped_weno": "T",
    "mp_weno": "T",
    "weno_eps": 1e-16,
    "avg_state": 2,
    "riemann_solver": 2,
    "wave_speeds": 1,
    "viscous": "F",
    # boundaries: extrapolation (non-reflecting-style; CBC/-6 is prohibited with eos_jwl)
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    # output
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "pres_wrt": "T",
    "parallel_io": "T",
    "lag_db_wrt": "T",
    # --- JWL products fluid (eos = 2) ---
    "fluid_pp(1)%eos": 2,
    "fluid_pp(1)%gamma": 1.0 / 0.4,
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%cv": jwl_Cv,
    "fluid_pp(1)%qv": 0.0,
    "fluid_pp(1)%qvp": 0.0,
    "fluid_pp(1)%G": 0.0,
    "fluid_pp(1)%jwl_A": jwl_A,
    "fluid_pp(1)%jwl_B": jwl_B,
    "fluid_pp(1)%jwl_R1": jwl_R1,
    "fluid_pp(1)%jwl_R2": jwl_R2,
    "fluid_pp(1)%jwl_omega": jwl_omega,
    "fluid_pp(1)%jwl_rho0": jwl_rho0,
    "fluid_pp(1)%jwl_Q": jwl_Q,
    "fluid_pp(1)%jwl_air_e0": 2.5575e5,
    "fluid_pp(1)%jwl_air_rho0": 1.225,
    # --- patch 1: background products (rectangle fills domain) ---
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.5 * (xmin + xmax),
    "patch_icpp(1)%y_centroid": 0.5 * (ymin + ymax),
    "patch_icpp(1)%length_x": xmax - xmin,
    "patch_icpp(1)%length_y": ymax - ymin,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": 10.0e9,
    "patch_icpp(1)%alpha_rho(1)": jwl_rho0,
    "patch_icpp(1)%alpha(1)": 1.0,
    # --- patch 2: high-pressure circular blast core ---
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%radius": 0.012,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": 30.0e9,
    "patch_icpp(2)%alpha_rho(1)": jwl_rho0,
    "patch_icpp(2)%alpha(1)": 1.0,
    # --- Euler-Lagrange particles ---
    "bubbles_lagrange": "F",
    "particles_lagrange": "T",
    "fd_order": 4,
    "lag_params%interpolation_order": 2,
    "lag_params%epsilonb": 1,
    "lag_params%charwidth": 5.0e-4,
    "lag_params%solver_approach": 2,
    "lag_params%nParticles_glb": 12000,
    "lag_params%valmaxvoid": 0.65,
    "lag_params%vel_model": 1,
    "lag_params%input_path": "input/lag_particles.dat",
    "lag_params%qs_drag_model": 2,
    "lag_params%stokes_drag": 0,
    "lag_params%added_mass_model": 1,
    "lag_params%pressure_force": "T",
    "lag_params%gravity_force": "F",
    "lag_params%collision_force": "T",
    "lag_params%qs_fluct_force": "F",
    "lag_params%write_bubbles": "T",
    "lag_params%mu_ref(1)": 1.7160e-5,
    "lag_params%suth(1)": 110.4,
    "particle_pp%rho0ref_particle": 1740.0,
    "particle_pp%cp_particle": 1000.0,
    "particle_pp%ksp_col": 10.0,
    "particle_pp%nu_col": 0.35,
    "particle_pp%E_col": 1.0e9,
    "particle_pp%cor_col": 0.7,
}

print(json.dumps(data, indent=4))
