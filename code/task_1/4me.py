# import matplotlib.pyplot as plt
# import numpy as np
#
# from src.utils import convert_primitive_to_conserved, convert_conserved_to_primitive
#
# w = np.array([0.5313, 0.8276, 0, 0.4])
# u = convert_primitive_to_conserved(w)
# w_new = convert_conserved_to_primitive(u=u)
#
# print(w)
# print(u)
# print(w_new)
#
# w = np.array([1, 0.1, 0, 1])
# u = convert_primitive_to_conserved(w)
# w_new = convert_conserved_to_primitive(u=u)
#
# print(w)
# print(u)
# print(w_new)
#
# w = np.array([0.8, 0.1, 0, 0.4])
# u = convert_primitive_to_conserved(w)
# w_new = convert_conserved_to_primitive(u=u)
#
# print(w)
# print(u)
# print(w_new)
#
# w = np.array([0.5313, 0.1, 0.7276, 0.4])
# u = convert_primitive_to_conserved(w)
# w_new = convert_conserved_to_primitive(u=u)
#
# print(w)
# print(u)
# print(w_new)

from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

def gen(n):
    phi = 0
    while phi < 2*np.pi:
        yield np.array([np.cos(phi), np.sin(phi), phi])
        phi += 2*np.pi/n

def update(num, data, line):
    line.set_data(data[:2, :num])
    line.set_3d_properties(data[2, :num])

N = 100
data = np.array(list(gen(N))).T
line, = ax.plot(data[0, 0:1], data[1, 0:1], data[2, 0:1])

print(data, line)

# Setting the axes properties
ax.set_xlim3d([-1.0, 1.0])
ax.set_xlabel('X')

ax.set_ylim3d([-1.0, 1.0])
ax.set_ylabel('Y')

ax.set_zlim3d([0.0, 10.0])
ax.set_zlabel('Z')

ani = animation.FuncAnimation(fig, update, N, fargs=(data, line), interval=10000/N, blit=False)
#ani.save('matplot003.gif', writer='imagemagick')
plt.show()