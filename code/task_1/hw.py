import numpy as np
from src.utils import convert_conserved_to_primitive
from src.euler_model import Euler
from src.simple_solver import SimpleSolver

import matplotlib.pyplot as plt

''' Configuration 11 '''
# w = [rho, u, v, p]
w_2 = np.array([0.5313, 0.8276, 0, 0.4])
w_1 = np.array([1, 0.1, 0, 1])
w_3 = np.array([0.8, 0.1, 0, 0.4])
w_4 = np.array([0.5313, 0.1, 0.7276, 0.4])

position_nodes_number = 50

x0 = 0.5
y0 = 0.5
t_final = 0.03

euler_model = Euler(initial_conditions_=[w_1, w_2, w_3, w_4, x0, y0], wave_speed_estimator_='pressure_based', dim=2)

simple_solver = SimpleSolver(model=euler_model, cfl=0.15)

u, space_x, space_y, t_final = simple_solver.solve_problem(position_nodes_number=position_nodes_number,
                                                           t_final=t_final)

time_nodes_number = np.shape(u)[0]

v = np.zeros((time_nodes_number, position_nodes_number, position_nodes_number, 4))
v_simple = np.zeros((time_nodes_number, position_nodes_number, position_nodes_number, 4))

for t in range(time_nodes_number):
    for y in range(position_nodes_number):
        for x in range(position_nodes_number):
            v[t, x, y, :] = convert_conserved_to_primitive(u[t, x, y, :])

fig = plt.figure(figsize=(14, 10))
fig.suptitle('Тест, ' + str(t_final), fontsize=16)

Z = v[-1][:, :, 0]
print(Z)
norm = plt.Normalize(Z.min(), Z.max())
# C = plt.cm.Blues_r(norm(Z)/2)

ax_1 = fig.add_subplot(111, projection='3d')
ax_1.plot_surface(space_x, space_y, Z)

plt.show()
# ax_1.scatter(space, v[-1][:, 0], s=5, c='blue', label='HLLC + ENO')
