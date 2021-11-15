import numpy as np
import pandas as pd
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

rho = v[-1][:, :, 0]
vel_x = v[-1][:, :, 1]
vel_y = v[-1][:, :, 2]
p = v[-1][:, :, 3]

out = np.zeros((np.shape(rho)[0], np.shape(rho)[1], 6))

print(space_x)
print(space_y)

print(np.shape(space_y))
print(np.shape(rho))
print(np.shape(out))

# with open("f.csv", "w") as f:
#     np.savetxt(f, delimiter=',', header="rho, vel_x, vel_y, p")

# marker = True
#
# for i in range(len(out)):
#     for j in range(len(out[0])):
#         out[i, j, 0] = space_x[j]
#         out[i, j, 1] = space_y[i]
#         out[i, j, 2] = rho[i, j]
#         out[i, j, 3] = vel_x[i, j]
#         out[i, j, 4] = vel_y[i, j]
#         out[i, j, 5] = p[i, j]
#
#         res = np.zeros((6, 1), dtype=np.double)
#         res[0][0] = space_x[j]
#         res[1][0] = space_y[i]
#         res[2][0] = rho[i, j]
#         res[3][0] = vel_x[i, j]
#         res[4][0] = vel_y[i, j]
#         res[5][0] = p[i, j]
#
#         res = res.transpose()
#
#         if marker:
#             with open("f.csv", "w") as f:
#                 np.savetxt(f, res, delimiter=',', header="x,y,rho,vel_x,vel_y,p")
#         else:
#             with open("f.csv", "a+") as f:
#                 np.savetxt(f, res, delimiter=',')
#
#         marker = False

    # marker = False
    # out[i, 4] = space_x[i]
    # out[i, 5] = space_y

# a = bd
# with open("f.csv", "w") as f:
#     np.savetxt(f, out, delimiter=',', header="rho, vel_x, vel_y, p")
