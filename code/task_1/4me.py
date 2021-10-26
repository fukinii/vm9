import matplotlib.pyplot as plt
import numpy as np

from src.utils import convert_primitive_to_conserved, convert_conserved_to_primitive

w = np.array([0.5313, 0.8276, 0, 0.4])
u = convert_primitive_to_conserved(w)
w_new = convert_conserved_to_primitive(u=u)

print(w)
print(u)
print(w_new)

w = np.array([1, 0.1, 0, 1])
u = convert_primitive_to_conserved(w)
w_new = convert_conserved_to_primitive(u=u)

print(w)
print(u)
print(w_new)

w = np.array([0.8, 0.1, 0, 0.4])
u = convert_primitive_to_conserved(w)
w_new = convert_conserved_to_primitive(u=u)

print(w)
print(u)
print(w_new)

w = np.array([0.5313, 0.1, 0.7276, 0.4])
u = convert_primitive_to_conserved(w)
w_new = convert_conserved_to_primitive(u=u)

print(w)
print(u)
print(w_new)
