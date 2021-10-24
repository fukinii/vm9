import numpy as np
from euler_model import Euler
from utils import convert_conserved_to_primitive, calc_speed_of_sound, calc_divided_difference
from abc import ABCMeta, abstractmethod


class Solver:
    def __init__(self, model: Euler, cfl=0.5):
        self.model = model
        self.cfl = cfl

    def build_boundary_conditions_1p(self, u, ulr):
        if self.model.boundary_condition == 'transparent':
            u[0] = u[1]
            u[-1] = u[-2]
            ulr[0, :, 0] = u[0]
            ulr[-1, :, 1] = u[-1]
        indent = 0
        return u, ulr, indent

    @staticmethod
    def make_reconstruction_1p(u):
        return u[:], u[:]

    @staticmethod
    def calc_rhs(flux_left, flux_right, h):
        return -(flux_right - flux_left) / h

    def calc_delta_time(self, u, h):
        min_value = 1e9
        for i in range(len(u)):
            v_i = convert_conserved_to_primitive(u[i])
            velocity = v_i[1]
            down = calc_speed_of_sound(u[i]) + np.abs(velocity)
            min_value = min(min_value,
                            h / down)
        return self.cfl * min_value

    def solve_problem(self, position_nodes_number, t_final, xl=0., xr=1.):
        u = np.zeros((1, position_nodes_number, 3))

        flux = np.zeros((position_nodes_number + 1, 3))
        ulr = np.zeros((position_nodes_number + 1, 3, 2))

        rhs = np.zeros((position_nodes_number, 3))
        h = (xr - xl) / position_nodes_number

        position = np.linspace(xl + h / 2, xr - h / 2, position_nodes_number)

        u[0] = self.model.build_initial_conditions(position)
        time = 0
        t = 0
        while time < t_final:

            u[t], ulr, indent = self.build_boundary_conditions(u[t], ulr)
            delta_time = self.calc_delta_time(u[t], h)

            for cell_index in range(indent - 1, np.shape(position)[0] - indent):
                ulr[cell_index, :, 1], ulr[cell_index + 1, :, 0] = self.reconstruct(u[t], position,
                                                                                    cell_index, order=3)

                flux[cell_index, :] = self.model.flux_riemann_solver(ulr[cell_index, :, 0],
                                                                     ulr[cell_index, :, 1])
                if cell_index > 1:
                    rhs[cell_index - 1, :] = self.calc_rhs(flux[cell_index - 1, :], flux[cell_index, :], h)
            u_next = u[t] + delta_time * rhs
            u_next = np.reshape(u_next, (1, len(u_next), 3))
            u = np.concatenate((u, u_next), axis=0)
            t += 1
            time += delta_time
            print(time)
        return u, position, time

    @abstractmethod
    def reconstruct(self, u, position, cell_index, order):
        pass

    @abstractmethod
    def build_boundary_conditions(self, u, ulr):
        pass
