import dolfinx
from dolfinx.fem import petsc
from dolfinx.nls import petsc
import ufl
from basix.ufl import element, mixed_element

# import mpi4py
# mpi4py.rc.initialize = False
# mpi4py.rc.finalize = False
from mpi4py import MPI

import numpy as np

import os
from scipy.interpolate import interp1d

from Utilities import linstep, project


class Pretreatment:
    def __init__(self, verbose, show_plots):

        self._init_constants()
        self.verbose = verbose
        self.show_plots = show_plots

        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.num_procs = self.comm.Get_size()

        self.path_to_output = "output"
        self.path_to_data_files = os.path.join(self.path_to_output, "spatially_varying")

        os.makedirs(self.path_to_data_files, exist_ok=True)

    # ================================================================

    def _init_constants(self):

        # Physical Constants
        self.k_b = 1.380649e-23
        self.eta = 2.8e-4  # dynamic viscosity of solvent (water)
        self.R = 8.31446261815324

        # Diffusion coefficient parameters (originally nano-meters)
        self.d_pore = 15.0 / 1.0e9  # pore diameter
        self.r_xy = 0.34 / 1.0e9  # xylose molecular radius
        self.r_xo = 1.02 / 1.0e9  # xylooligomer molecular radius
        self.r_f = 0.28 / 1.0e9  # furfural molecular radius

        # Thermal conductivities
        self.k_s = 0.12  # solid
        self.k_l = 0.63  # liquid
        self.k_g = 0.0197  # gas
        self.h = 5000.0  # convective heat transfer coeff

        # Specific heat capacities
        self.C_s = 2000.0  # solid
        self.C_l = 4350.0  # liquid
        self.C_g = 1970.0  # gas

        # Mass densities
        self.rho_s = 1000.0  # solid
        self.rho_l = 1000.0  # liquid
        self.rho_g = 2.01928202521  # gas
        self.rho_x = 730.0  # xlyan
        self.rho_gl = 1500.0  # glucan
        self.rho_li = 1300.0  # lignin
        self.rho_os = 0.5 * (self.rho_gl + self.rho_li)  # other solids

        # Molecular weights (originally g/mol)
        self.M_x = 132.0 / 1.0e3  # xylan
        self.M_xy = 150.0 / 1.0e3  # xylose
        self.M_xo = 450.0 / 1.0e3  # xylooligomer
        self.M_f = 96.0 / 1.0e3  # furfural
        self.M_w = 18.01528 / 1.0e3  # water

        # Latent heat of steam
        self.L_cond = 2200.0e3

        # Condensation and evaporation rate (s^-1)
        self.k_bar = 10.0

    # ================================================================

    def generate_mesh(self, nn=15):
        self.length = 5e-3 / 2.0
        self.diam = 2.0 * self.length / 20.0

        self.mesh = dolfinx.mesh.create_interval(MPI.COMM_WORLD, nn, [0.0, self.length])
        self.ndim = self.mesh.topology.dim

    # ================================================================

    def build_functions(self, degree=2):
        d1p3 = element("P", self.mesh.basix_cell(), degree)
        self.S = dolfinx.fem.functionspace(self.mesh, d1p3)

        self.num_comp = 7

        d1p3_mix = mixed_element([d1p3 for k in range(self.num_comp)])
        self.Q = dolfinx.fem.functionspace(self.mesh, d1p3_mix)

        self.spaces = []
        self.maps = []
        for i in range(self.num_comp):
            space_i, map_i = self.Q.sub(i).collapse()
            self.spaces.append(space_i)
            self.maps.append(map_i)

        self.v = ufl.TestFunctions(self.Q)

        self.u = dolfinx.fem.Function(self.Q)
        self.u_n = dolfinx.fem.Function(self.Q)

        self.c_s_id = 0
        self.eps_l_id = 1
        self.fstar_x_id = 2
        self.cstar_xo_id = 3
        self.cstar_xy_id = 4
        self.cstar_f_id = 5
        self.T_id = 6

    # ================================================================

    def build_problem(self, ve_params):

        self._define_reaction_rates()

        self._set_initial_conditions(ve_params)
        self._set_boundary_conditions()

        self._build_general_forms()

        self._build_convection_terms()
        self._build_diffusion_terms()
        self._build_reaction_terms()
        self._build_forcing_terms()

        self._build_nonlinear_form()

    # ================================================================

    def _define_reaction_rates(self):
        # Reaction rates (originally in A:cm^3/mol/s and E_A:kJ/mol)
        self.k_xo = 1.0e9 * ufl.exp(-110.0e3 / (self.R * self.u[self.T_id]))
        self.k_x1 = 8.0e11 * ufl.exp(-130.0e3 / (self.R * self.u[self.T_id]))
        self.k_x2 = 2.5e8 * ufl.exp(-110.0e3 / (self.R * self.u[self.T_id]))
        self.k_f = 7.0e5 * ufl.exp(-98.0e3 / (self.R * self.u[self.T_id]))

    # ================================================================

    def _set_initial_conditions(self, ve_params):

        if ve_params is None:
            # self.c_sbulk = 0.14 * 1e3   # bulk steam concentration (from paper, now set by lookup table)
            # self.c_acid0 = 0.1 * 1e3  # initial acid concentration: Convert mol/self.length to mol/m^3
            # self.f_x0 = 0.26  # xylan mass fraction
            # self.eps_p0 = 0.8  # initial porosity
            # self.f_is0 = 0.44  # initial solid fraction
            # self.T_s = 423.15  # steam temperature
            # self.glucan_solid_fraction_0 = 0.4
            self.c_acid0 = 0.008624*1e2  # initial acid concentration: Convert mol/self.length to mol/m^3
            self.f_x0 = 0.227  # xylan mass fraction
            self.eps_p0 = 0.8  # initial porosity
            self.f_is0 = 0.428195  # initial solid fraction
            self.T_s = 190.0 + 273.15  # steam temperature
            self.glucan_solid_fraction_0 = 0.416
        else:
            self.c_acid0 = ve_params.pt_in["initial_acid_conc"]*1e6 # converting mol/mL into mol/m^3
            self.f_x0 = ve_params.feedstock["xylan_solid_fraction"]
            self.eps_p0 = ve_params.feedstock["initial_porosity"]
            self.f_is0 = ve_params.pt_in["initial_solid_fraction"]
            self.T_s = ve_params.pt_in["steam_temperature"]
            self.glucan_solid_fraction_0 = ve_params.feedstock["glucan_solid_fraction"]

        self.c_acid0 = dolfinx.fem.Constant(self.mesh, self.c_acid0)
        self.f_x0 = dolfinx.fem.Constant(self.mesh, self.f_x0)
        self.eps_p0 = dolfinx.fem.Constant(self.mesh, self.eps_p0)
        self.f_is0 = dolfinx.fem.Constant(self.mesh, self.f_is0)
        self.T_s = dolfinx.fem.Constant(self.mesh, self.T_s)
        self.glucan_solid_fraction_0 = dolfinx.fem.Constant(self.mesh, self.glucan_solid_fraction_0)

        # Obtain steam concentration from lookup table and add to dictionary
        pt_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
        steam_datafile = os.path.join(pt_path, "lookup_tables", "sat_steam_table.csv")
        steam_data = np.loadtxt(steam_datafile, delimiter=",", skiprows=1)

        # build interpolator interp_steam = interp.interp1d(temp_in_K, dens_in_kg/m3)
        interp_steam = interp1d(steam_data[:, 2], steam_data[:, 4])
        steam_density = interp_steam(self.T_s.value)

        # Convert c_sbulk to (mol/m^3) => density (kg/m^3) / molecular weight (kg/mol)
        self.c_sbulk = steam_density / self.M_w

        self.eps_lt = self.eps_p0  # initial liquid+gas (= 1-solid) volumetric fraction
        self.eps_l0 = 0.25  # initial liquid volumetric fraction
        self.T_0 = 300.0

        ic = dolfinx.fem.Function(self.S)

        # u[self.eps_l_id]
        ic.x.array[:] = self.eps_l0
        self.u.x.array[self.maps[self.eps_l_id]] = ic.x.array[:]
        self.u_n.x.array[self.maps[self.eps_l_id]] = ic.x.array[:]
        # assign(self.u.sub(self.eps_l_id), ic)
        # assign(self.u_n.sub(self.eps_l_id), ic)

        # u[self.fstar_x_id]
        self.fstar_x0 = self.f_x0.value * (1.0 - self.eps_p0.value)
        ic.x.array[:] = self.fstar_x0
        self.u.x.array[self.maps[self.fstar_x_id]] = ic.x.array[:]
        self.u_n.x.array[self.maps[self.fstar_x_id]] = ic.x.array[:]

        # u[self.T_id]
        ic.x.array[:] = self.T_0
        self.u.x.array[self.maps[self.T_id]] = ic.x.array[:]
        self.u_n.x.array[self.maps[self.T_id]] = ic.x.array[:]

    # ================================================================

    def _set_boundary_conditions(self):
        def xlo_wall(x):
            return np.isclose(x[0], 0.0)

        def xhi_wall(x):
            return np.isclose(x[0], self.length)

        xlo_facets = dolfinx.mesh.locate_entities_boundary(
            self.mesh, self.ndim - 1, xlo_wall
        )
        xhi_facets = dolfinx.mesh.locate_entities_boundary(
            self.mesh, self.ndim - 1, xhi_wall
        )

        xhi_c_s_dofs = dolfinx.fem.locate_dofs_topological(
            self.Q.sub(self.c_s_id), self.ndim - 1, xhi_facets
        )
        xhi_T_dofs = dolfinx.fem.locate_dofs_topological(
            self.Q.sub(self.T_id), self.ndim - 1, xhi_facets
        )

        self.bcs = []
        self.bcs.append(
            dolfinx.fem.dirichletbc(
                dolfinx.fem.Constant(self.mesh, self.c_sbulk),
                xhi_c_s_dofs,
                self.Q.sub(self.c_s_id),
            )
        )

        self.bcs.append(
            dolfinx.fem.dirichletbc(
                # dolfinx.fem.Constant(self.mesh, self.T_s),
                self.T_s,
                xhi_T_dofs,
                self.Q.sub(self.T_id),
            )
        )

    # ================================================================

    def _build_general_forms(self):

        aa = self.rho_x / self.rho_os
        # aa = 2.0 # TODO: This is set as 2.0 in the paper/OG code
        bb = self.f_x0 / (1.0 - self.f_x0) + aa
        cc = 1.0 - self.eps_p0

        ss1 = aa**2 * cc**2
        ss2 = 2.0 * (aa - 2.0) * bb * cc * self.u_n[self.fstar_x_id]
        ss3 = bb**2 * self.u_n[self.fstar_x_id] ** 2

        ss = ufl.sqrt(ss1 - ss2 + ss3)

        self.eps_p = -(aa * cc + bb * (self.u_n[self.fstar_x_id] - 2.0) + ss) / (
            2.0 * bb
        )
        self.eps_g = self.eps_p - self.u_n[self.eps_l_id]

        # --------------------------------

        self.c_acid = self.c_acid0 * self.eps_l0 / self.u_n[self.eps_l_id]

        # --------------------------------

        self.D_s = (
            self.eps_g
            * self.d_pore
            * ufl.sqrt((8.0 * self.R * self.u_n[self.T_id]) / (ufl.pi * self.M_w))
            / 3.0
        )

        self.D_xy_0 = (self.k_b * self.u_n[self.T_id]) / (
            6.0 * ufl.pi * self.eta * self.r_xy
        )
        self.D_xo_0 = (self.k_b * self.u_n[self.T_id]) / (
            6.0 * ufl.pi * self.eta * self.r_xo
        )
        self.D_f_0 = (self.k_b * self.u_n[self.T_id]) / (
            6.0 * ufl.pi * self.eta * self.r_f
        )

        self.D_xy = self.u_n[self.eps_l_id] * self.D_xy_0
        self.D_xo = self.u_n[self.eps_l_id] * self.D_xo_0
        self.D_f = self.u_n[self.eps_l_id] * self.D_f_0

        # --------------------------------

        self.rho_eff_c_eff = (
            (1.0 - self.eps_p) * self.rho_s * self.C_s
            + self.u_n[self.eps_l_id] * self.rho_l * self.C_l
            + self.eps_g * self.rho_g * self.C_g
        )
        self.k_eff = (
            (1.0 - self.eps_p) * self.k_s
            + self.u_n[self.eps_l_id] * self.k_l
            + self.eps_g * self.k_g
        )

        # self.heat_flux = (
        #     self.h
        #     / self.k_eff
        #     * (self.u_n[self.T_id] - self.T_s)
        #     * self.v[self.T_id]
        #     * self.ds(1)
        # )  # TODO: fix heat flux specification, currently not used in formulation

        # --------------------------------

        self.k_cond = dolfinx.fem.Function(self.S)
        self.k_cond.x.array[:] = self.k_bar

        self.k_evap = dolfinx.fem.Function(self.S)
        self.k_evap.x.array[:] = 0.0

    # ================================================================

    def _build_convection_terms(self):
        self.conv = [0 for k in range(self.num_comp)]

        self.conv[self.c_s_id] = 0
        self.conv[self.eps_l_id] = 0
        self.conv[self.fstar_x_id] = 0
        self.conv[self.cstar_xo_id] = (
            self.D_xo_0 / self.u_n[self.eps_l_id] * ufl.Dx(self.u_n[self.eps_l_id], 0)
        )
        self.conv[self.cstar_xy_id] = (
            self.D_xy_0 / self.u_n[self.eps_l_id] * ufl.Dx(self.u_n[self.eps_l_id], 0)
        )
        self.conv[self.cstar_f_id] = (
            self.D_f_0 / self.u_n[self.eps_l_id] * ufl.Dx(self.u_n[self.eps_l_id], 0)
        )
        # self.conv[self.cstar_xo_id] = self.D_xo_0 / self.u_n[self.eps_l_id] * ufl.grad(self.u_n[self.eps_l_id])
        # self.conv[self.cstar_xy_id] = self.D_xy_0 / self.u_n[self.eps_l_id] * ufl.grad(self.u_n[self.eps_l_id])
        # self.conv[self.cstar_f_id] = self.D_f_0 / self.u_n[self.eps_l_id] * ufl.grad(self.u_n[self.eps_l_id])
        self.conv[self.T_id] = 0

    # ================================================================

    def _build_diffusion_terms(self):
        self.diff = [0 for k in range(self.num_comp)]

        self.diff[self.c_s_id] = self.D_s
        self.diff[self.eps_l_id] = 0
        self.diff[self.fstar_x_id] = 0
        self.diff[self.cstar_xo_id] = self.D_xo_0
        self.diff[self.cstar_xy_id] = self.D_xy_0
        self.diff[self.cstar_f_id] = self.D_f_0
        self.diff[self.T_id] = self.k_eff / self.rho_eff_c_eff

    # ================================================================

    def _build_reaction_terms(self):
        self.reac = [0 for k in range(self.num_comp)]

        self.reac[self.c_s_id] = -self.k_cond * self.eps_g
        self.reac[self.eps_l_id] = (
            -self.M_w * self.k_cond * self.u_n[self.c_s_id] / self.rho_l - self.k_evap
        )
        self.reac[self.fstar_x_id] = -(self.k_xo + self.k_x1) * self.c_acid * self.eps_p
        self.reac[self.cstar_xo_id] = -self.k_x2 * self.c_acid
        self.reac[self.cstar_xy_id] = -self.k_f * self.c_acid
        self.reac[self.cstar_f_id] = 0
        self.reac[self.T_id] = -4.0 * self.h / (self.rho_eff_c_eff * self.diam)
        # self.reac[self.T_id] = -4.0 * self.h / self.diam

    # ================================================================

    def _build_forcing_terms(self):
        self.forc = [0 for k in range(self.num_comp)]

        self.forc[self.c_s_id] = (
            self.k_evap
            * (self.u_n[self.eps_l_id] - self.eps_lt)
            * self.rho_l
            / self.M_w
        )
        self.forc[self.eps_l_id] = (
            self.k_cond * self.u_n[self.c_s_id] * self.eps_p * self.M_w / self.rho_l
            + self.k_evap * self.eps_lt
        )
        self.forc[self.fstar_x_id] = 0
        self.forc[self.cstar_xo_id] = (
            self.k_xo
            * self.rho_s
            * self.c_acid
            * self.eps_p
            * self.u_n[self.fstar_x_id]
            / self.M_xo
        )
        self.forc[self.cstar_xy_id] = (
            self.k_x1
            * self.rho_s
            * self.c_acid
            * self.eps_p
            * self.u_n[self.fstar_x_id]
            / self.M_xy
            + self.k_x2
            * self.u_n[self.cstar_xo_id]
            * self.c_acid
            * self.M_xo
            / self.M_xy
        )
        self.forc[self.cstar_f_id] = (
            self.k_f * self.u_n[self.cstar_xy_id] * self.c_acid * self.M_xy / self.M_f
        )
        self.forc[self.T_id] = (
            self.L_cond
            / self.rho_eff_c_eff
            * (
                self.k_cond * self.u_n[self.c_s_id] * self.M_w * (self.eps_g)
                - self.k_evap * self.rho_l * (self.u_n[self.eps_l_id] - self.eps_lt)
            )
            + 4.0 * self.h / (self.rho_eff_c_eff * self.diam) * self.T_s
        )
        # self.forc[self.T_id] = (
        #     (self.L_cond / self.rho_eff_c_eff)
        #     * (
        #         self.k_cond * self.u_n[self.c_s_id] * self.M_w * (self.eps_g)
        #         - self.k_evap * self.rho_l * (self.u_n[self.eps_l_id] - self.eps_lt)
        #     )
        #     + 4.0 * self.h / self.diam * self.T_s
        # )

    # ================================================================

    def _build_nonlinear_form(self):
        self.F = 0

        # Use a placeholder value of dt to build the dolfinx.fem.form, it will
        # be updated with the correct value in solve(...)
        self.dt_c = dolfinx.fem.Constant(self.mesh, 1.0)

        for k in range(self.num_comp):

            self.F += (1.0 / self.dt_c) * (self.u[k] - self.u_n[k]) * self.v[k] * ufl.dx

            if self.conv[k] != 0:
                # BEST ONE (move the ufl.grad onto the test space and swap sign)
                self.F -= self.conv[k] * self.u[k] * ufl.Dx(self.v[k], 0) * ufl.dx

                # self.F += ufl.Dx(self.conv[k]*self.u[k], 0) * self.v[k] * ufl.dx
                # self.F -= ufl.inner(self.conv[k]*self.u[k], ufl.grad(self.v[k]))*ufl.dx
                # self.F += dot(self.conv[k], ufl.grad(self.u[k]))*self.v[k]*ufl.dx

            if self.diff[k] != 0:
                # BEST ONE (move the ufl.grad onto the test space and swap sign)
                self.F += (
                    self.diff[k] * ufl.Dx(self.u[k], 0) * ufl.Dx(self.v[k], 0) * ufl.dx
                )

                # self.F += self.diff[k] * ufl.inner(ufl.grad(self.u[k]), ufl.grad(self.v[k])) * ufl.dx
                # self.F += self.diff[k] * dot(ufl.grad(self.u[k]), ufl.grad(self.v[k])) * ufl.dx
                # self.F += ufl.inner(self.diff[k] * ufl.grad(self.u[k]), ufl.grad(self.v[k])) * ufl.dx

            if self.reac[k] != 0:
                self.F -= self.reac[k] * self.u[k] * self.v[k] * ufl.dx

            if self.forc[k] != 0:
                self.F -= self.forc[k] * self.v[k] * ufl.dx

        # TODO: fix heat flux specification, currently not used in formulation
        # self.F += self.heat_flux

    # ================================================================

    def solve(self, t_final=30.0, dt=1.0, save_every=2.0):

        # Time stepping parameters
        self.dt = dt
        self.dt_c.value = self.dt

        self.t_final = t_final
        t_steps = int(self.t_final / self.dt)
        save_every_n = int(save_every / self.dt)

        tt = 0
        self.write_solution(tt)

        du = ufl.TrialFunction(self.Q)
        J = ufl.derivative(self.F, self.u, du)

        nonlinear_problem = dolfinx.fem.petsc.NonlinearProblem(self.F, self.u, self.bcs)
        nonlinear_solver = dolfinx.nls.petsc.NewtonSolver(MPI.COMM_WORLD, nonlinear_problem)

        for k in range(t_steps):
            # Solve the nonlinear problem
            nonlinear_solver.solve(self.u)

            # Extract a numpy vector for temp and eps_l
            np_T = self.u.x.array[self.maps[self.T_id]]
            np_eps_l = self.u.x.array[self.maps[self.eps_l_id]]

            self._update_k_cond(np_T)
            self._update_k_evap(np_T, np_eps_l)

            # Update previous solution
            self.u_n.x.array[:] = self.u.x.array[:]

            tt = (k + 1) * self.dt

            if (k + 1) % save_every_n == 0 or (k+1) == t_steps:
                self.write_solution(tt)

    # ================================================================

    def _update_k_cond(self, temp):

        T_tol = 0.001
        del_T = T_tol * self.T_s.value
        out_vec = linstep(temp, self.T_s.value - del_T, self.T_s.value + del_T, self.k_bar, 0.0)

        self.k_cond.x.array[:] = out_vec

    # ================================================================

    def _update_k_evap(self, temp, liq):

        T_tol = 0.001
        del_T = T_tol * self.T_s.value
        out_vec = linstep(temp, self.T_s.value - del_T, self.T_s.value + del_T, 0.0, self.k_bar)

        eps_lt_tol = 0.1
        del_eps_lt = eps_lt_tol * self.eps_lt.value
        out_vec_2 = linstep(liq, self.eps_lt.value, self.eps_lt.value + del_eps_lt, 0.0, 1.0)
        out_vec *= out_vec_2
        # TODO: needed to soften the transition to nonzero k_evap on
        # only one side of eps_lt

        relaxation_method = 3

        if relaxation_method == 1:
            k_evap_scaling = 0.05
            k_evap_new = k_evap_scaling * out_vec

        elif relaxation_method == 2:
            k_evap_old = self.k_evap.vector()[:]

            relaxation_time = 10.0
            delta_k_evap = 1.0 / relaxation_time * (out_vec - k_evap_old)

            k_evap_new = k_evap_old + delta_k_evap

            # k_evap_new *= 0.4

        elif relaxation_method == 3:
            k_evap_new = out_vec

        self.k_evap.x.array[:] = k_evap_new

    # ================================================================

    def write_solution(self, tt):

        self._save_spatially_varying_quantities(tt)
        self._write_integrated_quantities(tt)

    def _save_spatially_varying_quantities(self, tt):
        output_array = []

        data_header = "X-position (cm),"
        xpts = self.S.tabulate_dof_coordinates()[:, 0]
        xpts = xpts * 100.0
        xpts_id = np.argsort(xpts)
        output_array.append(xpts[xpts_id])

        # Steam (M)
        data_header += "Steam (M),"
        c_s_vec = self.u.x.array[self.maps[self.c_s_id]]
        c_s_vec = c_s_vec / 1000.0
        output_array.append(c_s_vec[xpts_id])

        # Liquid Fraction
        data_header += "Liquid Fraction,"
        eps_l_vec = self.u.x.array[self.maps[self.eps_l_id]]
        output_array.append(eps_l_vec[xpts_id])

        # Temperature (K)
        data_header += "Temperature (K),"
        T_vec = self.u.x.array[self.maps[self.T_id]]
        output_array.append(T_vec[xpts_id])

        # Acid (M)
        data_header += "Acid (M),"
        acid_fn = dolfinx.fem.Function(self.S)
        acid_fn = project(self.c_acid, acid_fn)
        acid_vec = acid_fn.x.array[:] / 1000.0
        output_array.append(acid_vec[xpts_id])

        # Xylan Fraction
        data_header += "Xylan Fraction,"
        f_x_fn = dolfinx.fem.Function(self.S)
        f_x_fn = project(self.u[self.fstar_x_id] / (1.0 - self.eps_p), f_x_fn)
        f_x_vec = f_x_fn.x.array[:]
        output_array.append(f_x_vec[xpts_id])

        # Xylooligomers (M)
        data_header += "Xylooligomer (M),"
        c_xo_vec = self.u.x.array[self.maps[self.cstar_xo_id]]
        c_xo_vec = c_xo_vec / eps_l_vec / 1000.0
        output_array.append(c_xo_vec[xpts_id])

        # Xylose (M)
        data_header += "Xylose (M),"
        c_xy_vec = self.u.x.array[self.maps[self.cstar_xy_id]]
        c_xy_vec = c_xy_vec / eps_l_vec / 1000.0
        output_array.append(c_xy_vec[xpts_id])

        # Furfural (mM)
        data_header += "Furfural (mM)"
        c_f_vec = self.u.x.array[self.maps[self.cstar_f_id]]
        c_f_vec = c_f_vec / eps_l_vec
        output_array.append(c_f_vec[xpts_id])

        output_array = np.array(output_array).T
        output_array_filename = os.path.join(
            self.path_to_data_files, f"data_t{tt:05.0f}s.csv"
        )

        np.savetxt(
            output_array_filename,
            output_array,
            delimiter=",",
            header=data_header,
        )

    def _write_integrated_quantities(self, tt):
        eps_l_bar = dolfinx.fem.assemble_scalar(
            dolfinx.fem.form(self.u[self.eps_l_id] * ufl.dx)
        )
        eps_p_bar = dolfinx.fem.assemble_scalar(
            dolfinx.fem.form((1.0 - self.eps_p) * ufl.dx)
        )

        # This is equivalent to dolfinx.fem.assemble_scalar(...) since fstar_x0 is scalar
        fstar_x0_bar = self.fstar_x0 * self.length / eps_p_bar
        fstar_x_bar = (
            dolfinx.fem.assemble_scalar(
                dolfinx.fem.form(self.u[self.fstar_x_id] * ufl.dx)
            )
            / eps_p_bar
        )
        cstar_xo_bar = (
            dolfinx.fem.assemble_scalar(
                dolfinx.fem.form(self.u[self.cstar_xo_id] * ufl.dx)
            )
            / eps_l_bar
        )
        cstar_xy_bar = (
            dolfinx.fem.assemble_scalar(
                dolfinx.fem.form(self.u[self.cstar_xy_id] * ufl.dx)
            )
            / eps_l_bar
        )
        cstar_f_bar = (
            dolfinx.fem.assemble_scalar(
                dolfinx.fem.form(self.u[self.cstar_f_id] * ufl.dx)
            )
            / eps_l_bar
        )

        fis = self.rho_s * eps_p_bar / (self.rho_s * eps_p_bar + self.rho_l * eps_l_bar)

        converted_mass = fstar_x0_bar - fstar_x_bar
        conversion_percent = converted_mass / fstar_x0_bar


        # Get the total mass of the solid and liquid phase
        mass_of_solid_phase = self.rho_s * eps_p_bar
        mass_of_liquid_phase = self.rho_l * eps_l_bar

        # Get the mass of each of the soluble species (xylog, xylose, furfural)
        cstar_xo_mass = cstar_xo_bar * self.M_xo * eps_l_bar
        cstar_xy_mass = cstar_xy_bar * self.M_xy * eps_l_bar
        cstar_f_mass = cstar_f_bar * self.M_f * eps_l_bar

        # Get the initial and current mass of xylan
        fstar_x0_mass = fstar_x0_bar*mass_of_solid_phase
        fstar_x_mass = fstar_x_bar*mass_of_solid_phase

        # Sum the soluble species' mass, total_mass_of_products should be <= fstar_x0_mass
        total_mass_of_products = cstar_xo_mass + cstar_xy_mass + cstar_f_mass

        if self.verbose:
            print(f"Mass of xylan initial (g): {fstar_x0_mass}")
            print(f"Mass of xylan (g): {fstar_x_mass}")

            print(f"Unadjusted mass of xylog (g): {cstar_xo_mass}")
            print(f"Unadjusted mass of xylose (g): {cstar_xy_mass}")
            print(f"Unadjusted mass of furfural (g): {cstar_f_mass}")

            print(f"Unadjusted mass of products (g): {total_mass_of_products}")
            print(f"Unadjusted mass of products / initial xylan (g): {total_mass_of_products/fstar_x0_mass}")

        # Since mass conservation may not be exactly enforced, scale the soluble species masses
        # using the actual reacted mass, i.e., (fstar_x0_mass * conversion_percent)
        try:
            adjusted_cstar_xo_mass = cstar_xo_mass/total_mass_of_products * fstar_x0_mass * conversion_percent
            adjusted_cstar_xy_mass = cstar_xy_mass/total_mass_of_products * fstar_x0_mass * conversion_percent
            adjusted_cstar_f_mass = cstar_f_mass/total_mass_of_products * fstar_x0_mass * conversion_percent
        except:
            adjusted_cstar_xo_mass = np.nan
            adjusted_cstar_xy_mass = np.nan
            adjusted_cstar_f_mass = np.nan

        # Compute the new soluble species sum after the adjustment, now guaranteed <= fstar_x0_mass
        adjusted_total_mass_of_products = adjusted_cstar_xo_mass + adjusted_cstar_xy_mass + adjusted_cstar_f_mass

        # Compute the fractional conversion amounts on a molar basis
        frac_conv_xylog = (adjusted_cstar_xo_mass/self.M_xo) / (fstar_x0_mass/self.M_x)
        frac_conv_xylose = (adjusted_cstar_xy_mass/self.M_xy) / (fstar_x0_mass/self.M_x)
        frac_conv_furfural = (adjusted_cstar_f_mass/self.M_f) / (fstar_x0_mass/self.M_x)

        if self.verbose:
            print(f"Adjusted mass of xylog (g): {adjusted_cstar_xo_mass}")
            print(f"Adjusted mass of xylose (g): {adjusted_cstar_xy_mass}")
            print(f"Adjusted mass of furfural (g): {adjusted_cstar_f_mass}")

            print(f"Adjusted mass of products (g): {adjusted_total_mass_of_products}")
            print(f"Adjusted mass of products / initial xylan (g): {adjusted_total_mass_of_products/fstar_x0_mass}")

            print(f"Fractional conversion of xylog {frac_conv_xylog}")
            print(f"Fractional conversion of xylose {frac_conv_xylose}")
            print(f"Fractional conversion of furfural {frac_conv_furfural}")


        glucan_solid_fraction = self.glucan_solid_fraction_0.value / (
            1.0 - self.f_x0.value * conversion_percent
        )

        # Olga's version from intuition (very close to above!)
        # glucan_f = (1.0 - fstar_x_bar)*(self.glucan_solid_fraction_0)/(1.0 - self.f_x0)

        if not hasattr(self, "integrated_quantities"):
            self.integrated_quantities = {}

        self.integrated_quantities["fis_0"] = fis
        self.integrated_quantities["conv"] = conversion_percent
        self.integrated_quantities["X_X"] = fstar_x_bar
        self.integrated_quantities["X_G"] = glucan_solid_fraction
        self.integrated_quantities["rho_x"] = (
            cstar_xo_bar * self.M_xo + cstar_xy_bar * self.M_xy
        )
        self.integrated_quantities["rho_f"] = cstar_f_bar * self.M_f


        self.integrated_quantities["frac_conv_xylog"] = frac_conv_xylog
        self.integrated_quantities["frac_conv_xylose"] = frac_conv_xylose
        self.integrated_quantities["frac_conv_furfural"] = frac_conv_furfural


        if not hasattr(self, "integrated_quantities_array"):
            self.integrated_quantities_array = []

        data_header = [
            "Time (s)",
            "Xylan (g/g)",
            "Xylog (g/L)",
            "Xylose (g/L)",
            "Furfural (g/L)",
            "FIS (g/g)",
            "rho_x (g/g)",
            "X_G (g/g)",
            "conv",
            "frac_conv_xylog",
            "frac_conv_xylose",
            "frac_conv_furfural"
        ]

        self.integrated_quantities_array.append(
            [
                tt,
                fstar_x_bar,
                cstar_xo_bar * self.M_xo,
                cstar_xy_bar * self.M_xy,
                cstar_f_bar * self.M_f,
                fis,
                cstar_xo_bar * self.M_xo + cstar_xy_bar * self.M_xy,
                glucan_solid_fraction,
                conversion_percent,
                frac_conv_xylog,
                frac_conv_xylose,
                frac_conv_furfural
            ]
        )

        if self.verbose:
            for label, value in zip(data_header, self.integrated_quantities_array[-1]):
                print(f"{label}: {value}")
            print()

        output_array = np.array(self.integrated_quantities_array)
        output_array_filename = os.path.join(
            self.path_to_output, f"spatially_averaged_timeseries.csv"
        )

        np.savetxt(
            output_array_filename,
            output_array,
            delimiter=",",
            header=",".join(data_header),
        )

    # ================================================================
