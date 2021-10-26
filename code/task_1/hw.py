import numpy as np
from src.utils import convert_conserved_to_primitive
from src.euler_model import Euler
from src.simple_solver import SimpleSolver
from mpl_toolkits.mplot3d import Axes3D
import pylab
import pickle

import matplotlib.pyplot as plt

''' Configuration 11 '''
# w = [rho, u, v, p]
w_2 = np.array([0.5313, 0.8276, 0, 0.4])
w_1 = np.array([1, 0.1, 0, 1])
w_3 = np.array([0.8, 0.1, 0, 0.4])
w_4 = np.array([0.5313, 0.1, 0.7276, 0.4])

position_nodes_number = 100

x0 = 0.5
y0 = 0.5
t_final = 0.01

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

out = {'u': u, 'v': v, 'space_x': space_x, 'space_y': space_y, 't_final': t_final}

with open('v_50_50.pickle', 'wb') as f:
    pickle.dump(out, f)

# with open('v_50_50.pickle', 'rb') as f:
#      out = pickle.load(f)

u = out['u']
v = out['v']
space_x = out['space_x']
space_y = out['space_y']
t_final = out['t_final']

Z = v[-1][:, :, 0]
Z_1d = Z.ravel()
space_x_1d = np.tile(space_x, len(space_x))
space_y_1d = np.repeat(space_y, len(space_y))
a = 1

# print(space_x.shape)
# print(space_y.shape)
# print(Z.shape)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(space_x_1d, space_y_1d, Z_1d, marker='o')

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()
