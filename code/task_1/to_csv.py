import numpy as np
import pickle

with open('v_100_100_05_2d_01cfl.pickle', 'rb') as f:
    out = pickle.load(f)

u = out['u']
v = out['v']
space_x = out['space_x']
space_y = out['space_y']
t_final = out['t_final']

rho = v[-1][:, :, 0]
vel_x = v[-1][:, :, 1]
vel_y = v[-1][:, :, 2]
p = v[-1][:, :, 3]

out = np.zeros((np.shape(rho)[0], np.shape(rho)[1], 6))


# with open("v_100_100_15_1d_x_01cfl.csv", "w") as f:
#     np.savetxt(f, delimiter=',', header="rho, vel_x, vel_y, p")

marker = True

for i in range(len(out)):
    for j in range(len(out[0])):
        out[i, j, 0] = space_x[j]
        out[i, j, 1] = space_y[i]
        out[i, j, 2] = rho[i, j]
        out[i, j, 3] = vel_x[i, j]
        out[i, j, 4] = vel_y[i, j]
        out[i, j, 5] = p[i, j]

        res = np.zeros((6, 1), dtype=np.double)
        res[0][0] = space_x[j]
        res[1][0] = space_y[i]
        res[2][0] = rho[i, j]
        res[3][0] = vel_x[i, j]
        res[4][0] = vel_y[i, j]
        res[5][0] = p[i, j]

        res = res.transpose()

        if marker:
            with open("v_100_100_05_2d_01cfl.csv", "w") as f:
                np.savetxt(f, res, delimiter=',', header="x,y,rho,vel_x,vel_y,p")
        else:
            with open("v_100_100_05_2d_01cfl.csv", "a+") as f:
                np.savetxt(f, res, delimiter=',')

        marker = False



# print(space_x)
# print(space_y)
#
# space_x_1d = np.tile(space_x, len(space_x))
# space_y_1d = np.repeat(space_y, len(space_y))
#
# print(space_x_1d)
# print(space_y_1d)
#
# Z = v[-1][:, :, 3]
# for j in range(space_y):
#     for i in range(space_x):
#         out = np.array([space_x[]])
# np.savetxt("data3.csv", Z,
#            delimiter=",")
# with open('v_100_100_15_1d_x_01cfl.csv', 'wb') as f:
#     pickle.dump(out, f)
