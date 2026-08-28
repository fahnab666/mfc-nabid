#!/usr/bin/env python3
"""Daoud-Balachandar 2D moving-IBM curtain in a fixed three-level AMR corridor.

Mesh layout:
  L0: far upstream and downstream, extending to the domain boundaries.
  L1: fixed intermediate-resolution region closer to the curtain corridor.
  L2: fixed finest region containing the initial bed and expected travel.

The AMR hierarchy is stationary (amr_regrid_int=0), while the IBM particles
are force-driven (moving_ibm=2). Current MFC marks force-driven IBM with AMR
as unvalidated, and level-2 AMR with IBM must be run on one MPI rank.

Daoud inputs retained: Mach 1.66 shock, x=[-0.2,0.3] m, shock x=-0.005 m,
2 mm curtain, 115 micrometre particles, loading 0.21, particle density
2520 kg/m^3, restitution 0.7, and collision friction 0.4. Each 2D IBM circle
uses the physical mass of one 115 micrometre soda-lime sphere.
"""

import json
import math
import sys

GAMMA = 1.4
R_AIR = 287.058
X_BEG, X_END = -0.200, 0.300
Y_BEG, Y_END = -0.004, 0.004
SHOCK_X = -0.005

P1, T1, U1 = 82.7e3, 297.0, 0.0
RHO1 = P1 / (R_AIR * T1)
P2, RHO2, T2, U2 = 252.082e3, 2.073, 392.0, 304.146
MACH_SHOCK, SHOCK_SPEED = 1.66, 573.3138

PARTICLE_DIAMETER = 115.0e-6
PARTICLE_RADIUS = 0.5 * PARTICLE_DIAMETER
PARTICLE_DENSITY = 2520.0
PARTICLE_MASS = PARTICLE_DENSITY * math.pi * PARTICLE_DIAMETER**3 / 6.0

CURTAIN_X0, CURTAIN_X1 = 0.0, 0.002
CURTAIN_Y0, CURTAIN_Y1 = -0.00075, 0.00075
TARGET_LOADING = 0.21
curtain_area = (CURTAIN_X1 - CURTAIN_X0) * (CURTAIN_Y1 - CURTAIN_Y0)
particle_area = math.pi * PARTICLE_RADIUS**2
NUM_PARTICLES = round(TARGET_LOADING * curtain_area / particle_area)
ACHIEVED_LOADING = NUM_PARTICLES * particle_area / curtain_area

# L0/L1/L2 give approximately 5/10/20 cells per particle diameter.
L0_CELLS_PER_DIAMETER = 5
DX_TARGET = PARTICLE_DIAMETER / L0_CELLS_PER_DIAMETER
NX = round((X_END - X_BEG) / DX_TARGET)
NY = round((Y_END - Y_BEG) / DX_TARGET)
DX = (X_END - X_BEG) / NX
DY = (Y_END - Y_BEG) / NY

# Static L1 parent. MFC constructs L2 as its central half.
# Result: L1 ~[-7.5,18.5] mm and L2 ~[-1,12] mm in x.
L1_X0, L1_X1 = -0.0075, 0.0185
L1_Y0, L1_Y1 = -0.0020, 0.0020


def first_index(coord, origin, spacing):
    return max(0, math.ceil((coord - origin) / spacing - 0.5))


def last_index(coord, origin, spacing, count):
    return min(count - 1, math.floor((coord - origin) / spacing - 0.5))


BX0 = first_index(L1_X0, X_BEG, DX)
BX1 = last_index(L1_X1, X_BEG, DX, NX)
BY0 = first_index(L1_Y0, Y_BEG, DY)
BY1 = last_index(L1_Y1, Y_BEG, DY, NY)
AMR_STATIC_CAP = max(BX1 - BX0 + 1, BY1 - BY0 + 1)

# Reproduce the static L2 central-half construction for safety checks/reporting.
IX = max((BX1 - BX0 + 1) // 4, 10)
IY = max((BY1 - BY0 + 1) // 4, 10)
L2_BX0, L2_BX1 = BX0 + IX, BX1 - IX
L2_BY0, L2_BY1 = BY0 + IY, BY1 - IY
L2_X0, L2_X1 = X_BEG + L2_BX0 * DX, X_BEG + (L2_BX1 + 1) * DX
L2_Y0, L2_Y1 = Y_BEG + L2_BY0 * DY, Y_BEG + (L2_BY1 + 1) * DY

if not (L2_X0 < CURTAIN_X0 and L2_X1 >= 0.012):
    raise RuntimeError("L2 must cover the initial curtain and downstream travel")
if not (L2_Y0 < CURTAIN_Y0 and L2_Y1 > CURTAIN_Y1):
    raise RuntimeError("L2 must contain the complete transverse curtain strip")

MU_AIR = 1.835e-5
C2 = math.sqrt(GAMMA * P2 / RHO2)
MAX_SIGNAL_SPEED = max(math.sqrt(GAMMA * P1 / RHO1), U2 + C2)
DT_NOMINAL = 0.30 * min(DX, DY) / MAX_SIGNAL_SPEED
T_END = 1.0e-3
NUM_STEPS = math.ceil(T_END / DT_NOMINAL)
DT = T_END / NUM_STEPS
SAVE_EVERY = max(1, round(1.0e-5 / DT))
COLLISION_TIME = 20.0 * DT / 4.0

case = {
    "run_time_info": "T",
    "format": "silo",
    "precision": "double",
    "parallel_io": "T",
    "prim_vars_wrt": "T",
    "cons_vars_wrt": "T",
    "E_wrt": "T",
    "ib_state_wrt": "T",
    "rho_wrt": "T",
    "pres_wrt": "T",
    "vel_wrt(1)": "T",
    "vel_wrt(2)": "T",
    "mom_wrt(1)": "T",
    "mom_wrt(2)": "T",
    "alpha_wrt(1)": "T",
    "alpha_rho_wrt(1)": "T",
    "gamma_wrt": "T",
    "pi_inf_wrt": "T",
    "c_wrt": "T",
    "omega_wrt(3)": "T",
    "schlieren_wrt": "T",
    "schlieren_alpha(1)": 1.0,
    "qm_wrt": "T",
    "liutex_wrt": "T",
    "sim_data": "T",
    # cf_wrt requires surface tension and down_sample is 3D IGR-only.
    "cf_wrt": "F",
    "down_sample": "F",
    "x_domain%beg": X_BEG,
    "x_domain%end": X_END,
    "y_domain%beg": Y_BEG,
    "y_domain%end": Y_END,
    "m": NX - 1,
    "n": NY - 1,
    "p": 0,
    "cyl_coord": "F",
    "dt": DT,
    "t_step_start": 0,
    "t_step_stop": NUM_STEPS,
    "t_step_save": SAVE_EVERY,
    "num_patches": 2,
    "model_eqns": "5eq",
    "num_fluids": 1,
    "alt_soundspeed": "F",
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": "rk3",
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "weno_Re_flux": "T",
    "weno_avg": "T",
    "avg_state": "arithmetic",
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "T",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "viscous": "T",
    "fd_order": 4,
    # AMR rejects characteristic boundaries; use fixed left state and
    # extrapolation outflow. The transverse direction is periodic.
    "bc_x%beg": -17,
    "bc_x%end": -3,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    # Ambient pre-shock air.
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.5 * (X_BEG + X_END),
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": X_END - X_BEG,
    "patch_icpp(1)%length_y": Y_END - Y_BEG,
    "patch_icpp(1)%vel(1)": U1,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": P1,
    "patch_icpp(1)%alpha_rho(1)": RHO1,
    "patch_icpp(1)%alpha(1)": 1.0,
    # Post-shock gas to the left of x=-5 mm.
    "patch_icpp(2)%geometry": 3,
    "patch_icpp(2)%x_centroid": 0.5 * (X_BEG + SHOCK_X),
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%length_x": SHOCK_X - X_BEG,
    "patch_icpp(2)%length_y": Y_END - Y_BEG,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%vel(1)": U2,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": P2,
    "patch_icpp(2)%alpha_rho(1)": RHO2,
    "patch_icpp(2)%alpha(1)": 1.0,
    "fluid_pp(1)%gamma": 1.0 / (GAMMA - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%Re(1)": 1.0 / MU_AIR,
    "ib": "T",
    "num_ibs": 0,
    "num_particle_clouds": 1,
    "many_ib_patch_parallelism": "T",
    "ib_neighborhood_radius": 2,
    "particle_cloud(1)%x_centroid": 0.5 * (CURTAIN_X0 + CURTAIN_X1),
    "particle_cloud(1)%y_centroid": 0.5 * (CURTAIN_Y0 + CURTAIN_Y1),
    "particle_cloud(1)%z_centroid": 0.0,
    "particle_cloud(1)%length_x": CURTAIN_X1 - CURTAIN_X0,
    "particle_cloud(1)%length_y": CURTAIN_Y1 - CURTAIN_Y0,
    "particle_cloud(1)%length_z": 0.0,
    "particle_cloud(1)%num_particles": NUM_PARTICLES,
    "particle_cloud(1)%radius": PARTICLE_RADIUS,
    "particle_cloud(1)%mass": PARTICLE_MASS,
    "particle_cloud(1)%min_spacing": 0.0,
    "particle_cloud(1)%moving_ibm": 2,
    "particle_cloud(1)%seed": 8173,
    "particle_cloud(1)%cloud_geometry": 1,
    "particle_cloud(1)%packing_method": 1,
    "particle_cloud(1)%periodic": 0,
    "collision_model": 1,
    "coefficient_of_restitution": 0.7,
    "ib_coefficient_of_friction": 0.4,
    "collision_time": COLLISION_TIME,
    # Fixed AMR hierarchy: far field L0, corridor L1, bed+travel L2.
    "amr": "T",
    "amr_block_beg(1)": BX0,
    "amr_block_end(1)": BX1,
    "amr_block_beg(2)": BY0,
    "amr_block_end(2)": BY1,
    "amr_regrid_int": 0,
    "amr_ref_ratio": 2,
    "amr_max_level": 2,
    "amr_max_blocks": 2,
    # Prevent old builds from sizing GPU scratch to half the long x-domain.
    # New builds also cap static allocations automatically.
    "amr_max_grid_size": AMR_STATIC_CAP,
    "amr_subcycle": "T",
}

sys.stderr.write(
    f"Daoud 2D IBM-AMR: Np={NUM_PARTICLES}, phi_2D={ACHIEVED_LOADING:.5f}, "
    f"rho_p={PARTICLE_DENSITY:.1f} kg/m^3, m_p={PARTICLE_MASS:.6e} kg; "
    f"L0={NX}x{NY}, L1 x=[{L1_X0:.4f},{L1_X1:.4f}] m, "
    f"L2 x=[{L2_X0:.6f},{L2_X1:.6f}] m, "
    f"L2 y=[{L2_Y0:.6f},{L2_Y1:.6f}] m. Run with one MPI rank.\n"
)

print(json.dumps(case, indent=4))
