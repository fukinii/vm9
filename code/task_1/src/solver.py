import numpy as np
from src.euler_model import Euler
from src.utils import convert_conserved_to_primitive, calc_speed_of_sound, calc_divided_difference
from abc import ABCMeta, abstractmethod


class Solver:
    def __init__(self, model: Euler, cfl=0.5):
        self.model = model
        self.cfl = cfl

    def build_boundary_conditions_1p(self, u, ulr):
        if self.model.boundary_condition == 'transparent':
            u[0] = u[1]
            u[-1] = u[-2]
            u[:, 0] = u[:, 1]
            u[:, -1] = u[:, -2]

            ulr[1, 1:-1, :, 0] = u[0]
            ulr[-2, 1:-1, :, 2] = u[-1]
            ulr[1:-1, 1, :, 3] = u[:, 0]
            ulr[1:-1, -2, :, 1] = u[:, -1]

            ulr[0, 1:-1, :, 2] = u[0]
            ulr[-1, 1:-1, :, 0] = u[-1]
            ulr[1:-1, 0, :, 1] = u[:, 0]
            ulr[1:-1, -1, :, 3] = u[:, -1]

        indent = 1
        return u, ulr, indent

    @staticmethod
    def make_reconstruction_1p(u):
        # a = 2
        b = u[:], u[:], u[:], u[:]
        # c = np.transpose(b)

        # d = 4
        return np.transpose(b)

    @staticmethod
    def calc_rhs(flux_left, flux_right, h):
        return -(flux_right - flux_left) / h

    @staticmethod
    def calc_fhs_2d(flux_up, flux_down, flux_left, flux_right, h):
        return -(flux_up + flux_right + flux_down + flux_left) / h

    def calc_delta_time(self, u, h):
        min_value = 1e9
        for i in range(len(u)):
            v_i = convert_conserved_to_primitive(u[i])
            velocity = v_i[1]
            down = calc_speed_of_sound(u[i]) + np.abs(velocity)
            min_value = min(min_value,
                            h / down)
        return self.cfl * min_value

    def solve_problem(self, position_nodes_number, t_final, xl=0., xr=1., yl=0., yr=1.):
        u = np.zeros((1, position_nodes_number, position_nodes_number, self.model.dim_sol))

        # flux = np.zeros((position_nodes_number + 1, self.model.dim_sol))
        flux = np.zeros((position_nodes_number, position_nodes_number, self.model.dim_sol, 4))
        ulr = np.zeros((position_nodes_number + 2, position_nodes_number + 2, self.model.dim_sol, 4))
        # ulr =
        rhs = np.zeros((position_nodes_number, position_nodes_number, self.model.dim_sol))

        h_x = (xr - xl) / position_nodes_number
        position_x = np.linspace(xl + h_x / 2, xr - h_x / 2, position_nodes_number)

        h_y = (yr - yl) / position_nodes_number
        position_y = np.linspace(yl + h_y / 2, yr - h_y / 2, position_nodes_number)

        u[0] = self.model.build_initial_conditions(position_x, position_y)
        time = 0
        t = 0
        while time < t_final:
            # delta_time = self.calc_delta_time(u[t][0], h_x)
            delta_time = h_x * 0.3
            u[t], ulr, indent = self.build_boundary_conditions(u[t], ulr)

            for cell_index_y in range(indent - 1, np.shape(position_y)[0]):
                for cell_index_x in range(indent - 1, np.shape(position_x)[0]):
                    ulr_cell_index_x = cell_index_x + 1
                    ulr_cell_index_y = cell_index_y + 1
                    ulr[ulr_cell_index_y][ulr_cell_index_x] = self.reconstruct(
                        u[t][cell_index_y][cell_index_x])

            for cell_index_y in range(indent - 1, np.shape(position_y)[0]):
                for cell_index_x in range(indent - 1, np.shape(position_x)[0]):
                    print(cell_index_y, cell_index_x)
                    # if cell_index_y == 0:
                    #     flux[cell_index_y, cell_index_y, :, :] = 0
                    # if cell_index_x != 0 and cell_index_y != 0 and cell_index_x != np.shape(position_y)[0] - 1 \
                    #         and cell_index_y != np.shape(position_x)[0] - 1:

                    ulr_cell_index_x = cell_index_x + 1
                    ulr_cell_index_y = cell_index_y + 1

                    if cell_index_x == 1 and cell_index_y == 1:
                        debug = 1

                    flux[cell_index_y, cell_index_x, :, 2] = self.model.flux_riemann_solver_y(
                        ulr[ulr_cell_index_y + 1, ulr_cell_index_x, :, 0],
                        ulr[ulr_cell_index_y, ulr_cell_index_x, :, 2],
                    )
                    assert not np.isnan(flux).any()

                    flux[cell_index_y, cell_index_x, :, 1] = self.model.flux_riemann_solver_x(
                        ulr[ulr_cell_index_y, ulr_cell_index_x + 1, :, 3],
                        ulr[ulr_cell_index_y, ulr_cell_index_x, :, 1],
                    )
                    assert not np.isnan(flux).any()

                    flux[cell_index_y, cell_index_x, :, 0] = self.model.flux_riemann_solver_y(
                        ulr[ulr_cell_index_y - 1, ulr_cell_index_x, :, 2],
                        ulr[ulr_cell_index_y, ulr_cell_index_x, :, 0],
                    )
                    assert not np.isnan(flux).any()

                    flux[cell_index_y, cell_index_x, :, 3] = self.model.flux_riemann_solver_x(
                        ulr[ulr_cell_index_y, ulr_cell_index_x - 1, :, 1],
                        ulr[ulr_cell_index_y, ulr_cell_index_x, :, 3],
                    )
                    assert not np.isnan(flux).any()
                # else:
                #     flux[cell_index_y, cell_index_x, :, :] = 0

                    rhs[cell_index_y, cell_index_x] = self.calc_fhs_2d(
                        flux[cell_index_y, cell_index_x, :, 0],
                        flux[cell_index_y, cell_index_x, :, 1],
                        flux[cell_index_y, cell_index_x, :, 2],
                        flux[cell_index_y, cell_index_x, :, 3],
                        h_x)
                a = 1
            u_next = u[t] + delta_time * rhs
            print(u_next.shape)
            u_next = np.reshape(u_next, (1, len(u_next), len(u_next[0]), self.model.dim_sol))
            print(u_next.shape)
            print(u.shape)
            u = np.concatenate((u, u_next), axis=0)
            t += 1
            time += delta_time
            print(time)
        return u, position_x, position_y, time

    @abstractmethod
    def reconstruct(self, u):
        pass

    @abstractmethod
    def build_boundary_conditions(self, u, ulr):
        pass
