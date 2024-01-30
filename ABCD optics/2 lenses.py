# Imports
from math import *
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


sns.set_style("dark")

# Define variables
length = 0  # Set starting length to 0
spacing = 0.001
x_up = 10
x_low = 0
# y_up =10
# y_low =0
# length scales

nm = 10 ** -9
um = 10 ** -6
mm = 10 ** -3
km = 10 ** 3
Mm = 10 ** 6
Gm = 10 ** 9

'''1/q = -j*lambda/(pi*w^2) + 1/R'''

# Lambda = WL
WL = 679 * nm  # Wavelength

w = 700  * um  # Width of the beam at the source

R = 100000  # Radius of curvature at the source

# complex beam parameter
q_inv = ((-1j * WL / (pi * w ** 2)) + 1 / R)
q = q_inv ** (-1)
q_matrix = np.array([[q], [1]])
print("q = " + str(q))


# print(q_matrix)

# Objects


def none():
    """No object / Path"""
    ident = np.array([[1, 0], [0, 1]])
    return ident


def space(d):
    """d is the distance"""
    air = np.array([[1, d], [0, 1]])
    global length
    length = length + d
    return air


def lens(f):
    """f is the focal length of the lens"""
    lens_path = np.array([[1, 0], [-1 / f, 1]])
    return lens_path


def thick_lens(n_1, n_2, r_1, r_2, t):
    """n1 = refractive index outside the lens.
n2 = refractive index of the lens itself (inside the lens).
r1 = Radius of curvature of First surface.
r2 = Radius of curvature of Second surface.
t = center thickness of lens."""
    mat_1 = np.array([[1, 0], [-((n_2 - n_1) / (r_2 * n_1)), (n_2 / n_1)]])
    mat_2 = np.array([[1, t], [0, 1]])
    mat_3 = np.array([[1, 0], [-((n_2 - n_1) / (r_1 * n_2)), (n_1 / n_2)]])
    mat_4 = np.matmul(mat_1, mat_2)
    mat_5 = np.matmul(mat_4, mat_3)
    return mat_5


def curved_mirror(m):
    """n1 = initial refractive index
n2 = final refractive index."""

    r_e = m  # r_e = R cos(x) and/or R/cos(x)
    mir = np.array([[1, 0], [-2 / r_e, 1]])
    return mir


def refraction_flat(n_1, n_2):
    """n1 = initial refractive index
n2 = final refractive index."""

    mat_ref = np.array([[1, 0], [0, n_1 / n_2]])

    return mat_ref


def single_prism(k, d, n):
    """k=(cos psi /cos phi )
    is the beam expansion factor, where phi
   is the angle of incidence,
    psi  is the angle of refraction,
     d = prism path length,
     n = refractive index of the prism material.
      This matrix applies for orthogonal beam exit."""
    mat_prism = np.array([[k, d / (n * k)], [0, 1 / k]])
    return mat_prism


# Components
start = 100 * mm
middle = 100 * mm
end = 270 * mm
focal_1 = 1000 * mm
focal_2 = 300 * mm
first = space(start)
second = lens(focal_1)
third = space(middle)
fourth = lens(focal_2)
fifth = space(end)

# Combinations 1st on the right, last on the left

step_1 = np.matmul(second, first)
step_2 = np.matmul(third, step_1)
step_3 = np.matmul(fourth, step_2)
step_4 = np.matmul(fifth, step_3)

# More steps can be added

"""q' = (Aq + B)/(Cq + D)"""

q_prime_matrix = np.matmul(step_4, q_matrix)
normal_q_prime_matrix = q_prime_matrix / q_prime_matrix[1]

q_prime = normal_q_prime_matrix[0]
R_prime = (np.real(1 / q_prime)) ** -1
w_prime = (-np.imag(1 / q_prime) * pi / WL) ** (-1 / 2)

print('q\' = ' + str(q_prime))
print('R\' = ' + str(R_prime))
print('w\' = ' + str(w_prime))

# Intensity
# need a value for r
r = 1  # Should be calculated for you what is this??
P = 1  # Power of the laser

I_prime = (P / (np.absolute(q_prime) ** 2)) * exp(-(2 * r ** 2) / (w_prime ** 2))

print('The intensity at ' + str(r) + ' m is proportional to: ' + str(I_prime))
# These are the problem!! I need to use a different array!!!


x_array_before_lens_1 = np.arange(0, start, spacing)

x_array_after_lens_1 = np.arange(0, middle, spacing)

x_array_after_lens_2 = np.arange(0, end, spacing)

plot_array = []
plot_R_array = []

for x in x_array_before_lens_1:
    first = space(x)
    q_prime_matrix = np.matmul(first, q_matrix)
    normal_q_prime_matrix = q_prime_matrix / q_prime_matrix[1]

    q_prime = normal_q_prime_matrix[0]
    inv_q_prime = q_prime ** (-1)
    w_prime = (-1 * np.imag(inv_q_prime) * np.pi / WL) ** (-1 / 2)
    R_prime = (np.real(inv_q_prime)) ** (-1)

    plot_array.append(float(w_prime))
    plot_R_array.append(float(R_prime))

first = space(start)
plot_array_after_1 = []
plot_R_array_after_1 = []
plot_array_after_2 = []
plot_R_array_after_2 = []

for x in x_array_after_lens_1:
    third = space(x)
    step_2 = np.matmul(third, step_1)
    q_prime_matrix = np.matmul(step_2, q_matrix)
    normal_q_prime_matrix = q_prime_matrix / q_prime_matrix[1]  #
    q_prime = normal_q_prime_matrix[0]
    inv_q_prime = q_prime ** (-1)
    w_prime = (-1 * np.imag(inv_q_prime) * np.pi / WL) ** (-1 / 2)
    R_prime = (np.real(inv_q_prime)) ** (-1)

    plot_array_after_1.append(float(w_prime))
    plot_R_array_after_1.append(float(R_prime))

# reset
first = space(start)
second = lens(focal_1)
third = space(middle)
fourth = lens(focal_2)
step_1 = np.matmul(second, first)
step_2 = np.matmul(third, step_1)
step_3 = np.matmul(fourth, step_2)

for x in x_array_after_lens_2:
    fifth = space(x)
    step_4 = np.matmul(fifth, step_3)
    q_prime_matrix = np.matmul(step_4, q_matrix)
    normal_q_prime_matrix = q_prime_matrix / q_prime_matrix[1]

    q_prime = normal_q_prime_matrix[0]
    inv_q_prime = q_prime ** (-1)
    w_prime = (-1 * np.imag(inv_q_prime) * np.pi / WL) ** (-1 / 2)
    R_prime = (np.real(inv_q_prime)) ** (-1)

    plot_array_after_2.append(float(w_prime))
    plot_R_array_after_2.append(float(R_prime))

total = start + middle + end
x_array = np.arange(0, total, spacing)
x_array_after_lens_1 = np.arange(0, middle, spacing) + start
x_array_after_lens_2 = np.arange(0, end, spacing) + start + middle
negative_plot_1 = []
negative_plot_after_1 = []
negative_plot_after_2 = []
for y in plot_array:
    negative_plot_1.append(-1 * float(y))
for y in plot_array_after_1:
    negative_plot_after_1.append(-1 * float(y))

# After lens


for y in plot_array_after_2:
    negative_plot_after_2.append(-1 * float(y))

plt.plot(x_array_before_lens_1, plot_array, color="blue")
plt.plot(x_array_after_lens_1, plot_array_after_1, color="red")
plt.plot(x_array_before_lens_1, negative_plot_1, color="blue")
plt.plot(x_array_after_lens_1, negative_plot_after_1, color="red")

plt.plot(x_array_after_lens_2, negative_plot_after_2, color="green")
plt.plot(x_array_after_lens_2, plot_array_after_2, color="green")

plt.axvline(start, color="black")
plt.axvline(middle + start, color="black")
# plt.Circle((0,0),5)

# constrain axis
# plt.xlim(x_low, x_up)
# plt.ylim(y_low, y_up)


plt.show()
