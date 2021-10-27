import numpy as np
import pickle
import matplotlib.pyplot as plt

with open('v_100_100_15_changed.pickle', 'rb') as f:
     out = pickle.load(f)

u = out['u']
v = out['v']
space_x = out['space_x']
space_y = out['space_y']
t_final = out['t_final']


# v[-1][:, :, a]
# a = 0 - плотность
# a = 1 - скорость по x
# a = 2 - скорость по y
# a = 3 - давление

Z = v[-1][:, :, 0]
Z_1d = Z.ravel()
space_x_1d = np.tile(space_x, len(space_x))
space_y_1d = np.repeat(space_y, len(space_y))

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(space_x_1d, space_y_1d, Z_1d)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()
