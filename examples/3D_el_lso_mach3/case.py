#!/usr/bin/env python3
import json
import math
import os

# Euler-Lagrange solid-particle bed in a Mach-3 ideal-gas stream.
# The generated input file contains 1000 stationary spheres at 20% volume fraction.
D = 0.1
radius = D / 2.0
num_particles = 1000
volume_fraction = 0.20
particle_volume = math.pi * D**3 / 6.0
lx = 25.0 * D
domain_volume = num_particles * particle_volume / volume_fraction
ly = math.sqrt(domain_volume / lx)
lz = ly

cells_per_diameter = 50
nx = round(lx / D * cells_per_diameter)
ny = round(ly / D * cells_per_diameter)
nz = round(lz / D * cells_per_diameter)

gamma = 1.4
rho_inf = 1.225
p_inf = 101325.0
c_inf = math.sqrt(gamma * p_inf / rho_inf)
u_inf = 3.0 * c_inf
dt = 0.25 * D / cells_per_diameter / (u_inf + c_inf)

case = {
    "run_time_info": "T",
    "x_domain%beg": -lx / 2.0,
    "x_domain%end": lx / 2.0,
    "y_domain%beg": -ly / 2.0,
    "y_domain%end": ly / 2.0,
    "z_domain%beg": -lz / 2.0,
    "z_domain%end": lz / 2.0,
    "m": nx - 1,
    "n": ny - 1,
    "p": nz - 1,
    "num_patches": 1,
    "model_eqns": "5eq",
    "num_fluids": 1,
    "time_stepper": "rk3",
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": 5000,
    "t_step_save": 250,
    "weno_order": 5,
    "fd_order": 4,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "mp_weno": "T",
    "weno_avg": "T",
    "weno_Re_flux": "T",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "avg_state": "arithmetic",
    "viscous": "T",
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "bc_z%beg": -3,
    "bc_z%end": -3,
    "particles_lagrange": "T",
    "particle_params%coupling": "one_way",
    "particle_params%nparticles_glb": num_particles,
    "particle_params%input_path": "input/particles.dat",
    "particle_params%vel_model": 0,
    "particle_params%qs_force": 1,
    "particle_params%interpolation_order": 2,
    "particle_params%pressure_gradient_force": "F",
    "particle_params%added_mass_force": 0,
    "particle_params%qs_fluct_force": "F",
    "particle_params%write_particles": "F",
    "particle_params%write_void_evol": "F",
    "particle_params%valmaxvoid": 0.9,
    "particle_pp%rho0ref_particle": 2500.0,
    "lso_filter": "T",
    "lso_filter_wrt": "T",
    "lso_stat_wrt": "F",
    "filter_sigma": D / 4.0,
    "lso_R_gas": 287.05,
    "lso_mu": 1.84e-5,
    "fluid_pp(1)%k_therm": 0.0262,
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "E_wrt": "T",
    "parallel_io": "T",
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_x": lx,
    "patch_icpp(1)%length_y": ly,
    "patch_icpp(1)%length_z": lz,
    "patch_icpp(1)%vel(1)": u_inf,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": p_inf,
    "patch_icpp(1)%alpha_rho(1)": rho_inf,
    "patch_icpp(1)%alpha(1)": 1.0,
    "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
    "fluid_pp(1)%eos": "ideal_gas",
    "fluid_pp(1)%cv": 287.05 / (gamma - 1.0),
    "fluid_pp(1)%Re(1)": rho_inf * u_inf * D / 1.84e-5,
}


def write_particle_file():
    input_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "input")
    os.makedirs(input_dir, exist_ok=True)
    with open(os.path.join(input_dir, "particles.dat"), "w", encoding="utf-8") as particle_file:
        for ix in range(10):
            x = -lx / 2.0 + radius + ix * (lx - D) / 9.0
            for iy in range(10):
                y = -ly / 2.0 + radius + iy * (ly - D) / 9.0
                for iz in range(10):
                    z = -lz / 2.0 + radius + iz * (lz - D) / 9.0
                    particle_file.write(f"{x:.8e} {y:.8e} {z:.8e} 0 0 0 {radius:.8e}\n")


write_particle_file()
print(json.dumps(case, indent=4))
