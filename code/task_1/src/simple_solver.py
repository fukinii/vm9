import numpy as np
from euler_model import Euler
from utils import convert_conserved_to_primitive, calc_speed_of_sound, calc_divided_difference
from solver import Solver


class SimpleSolver(Solver):

    def reconstruct(self, u, position, cell_index, order):
        return self.make_reconstruction_1p(u[cell_index, :])

    def build_boundary_conditions(self, u, ulr):
        return self.build_boundary_conditions_1p(u, ulr)
