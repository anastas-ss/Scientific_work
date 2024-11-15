import math

import numpy as np
from sympy.solvers import solve
import matplotlib.pyplot as plt
from sympy import Symbol
from mpmath import *

# Задаем параметры
Lz = 1
h = 10**(-3)
eps1, eps2 = 10**(-9), 10**(-12)
eps = 10**(-9)
i = complex(0, 1)
A1_next_n = A1_next = 1
A3_next_n = A3_next = 0
y = 1
N = Lz*10**3 - 1
delta_k = 1
z_all = np.arange(0, 5+h, h)
A_all = np.exp(-i * z_all)
z_values, A1_values, A3_values = [], [], []
I1, I3 = [], []
p3_array, cn_array = [], []
I30_array = []
b = delta_k/y


for z in z_all:
    while True:
        flag = 0
        A1_n_n = ((-1) * i * y * h * ((
                (A1_next_n + A1_next)/2) * ((abs(A1_next_n))**2 + abs(A1_next)**2)/2 + 2 * (
                (A1_next_n + A1_next)/2) * ((abs(A3_next_n))**2 + (abs(A3_next))**2)/2 + (
                (A3_next_n + A3_next)/2) * ((A1_next_n.conjugate())**2 + A1_next_n.conjugate() * A1_next.conjugate() + A1_next.conjugate()**2)/3) + A1_next
                  + i*h*(A1_next_n + A1_next)*delta_k/4)
        A3_n_n = ((-1) * i * y * h *((
                6*(A3_next_n + A3_next)/2) * (abs(A1_next_n) ** 2 + abs(A1_next) ** 2) / 2 + 3 * (
                (A3_next_n + A3_next) / 2) * ((abs(A3_next_n) ** 2 + (abs(A3_next)) ** 2) / 2) + (
                A1_next_n**3 + A1_next**3) / 2) + A3_next
                  + i*h*(A3_next_n + A3_next)*delta_k/4)

        if (abs(A1_n_n - A1_next_n) < eps1 * abs(A1_next_n) + eps2 and
                abs(A3_n_n - A3_next_n) < eps1 * abs(A3_next_n) + eps2):
            z_values.append(z)
            A1_values.append(A1_n_n)
            A3_values.append(A3_n_n)

            # проверка 1ого инварианта
            I1.append(abs(A1_n_n) ** 2 + abs(A3_n_n) ** 2)
            # проверка 3ого инварианта
            I3.append(3 / 2 * y * (
                    abs(A1_n_n) ** 4 + 4 * abs(A1_n_n) ** 2 * abs(A3_n_n) ** 2 + abs(A3_n_n) ** 4) + y * (
                              A3_n_n * (A1_n_n.conjugate()) ** 3 + (A3_n_n.conjugate()) * A1_n_n ** 3) - (
                                  3 * abs(A1_n_n) ** 2 + abs(A3_n_n) ** 2) * delta_k / 2)

            break

        A1_next_n = A1_n_n
        A3_next_n = A3_n_n

    # Шаг
    A1_next = A1_n_n
    A3_next = A3_n_n

I30 = 2 * I3[0] - 3 * y * I1[0] ** 2 + 3 * delta_k * I1[0]
print("z = ", z_all[0])
print("I1 = ", I1[0])
print("I3 = ", I3[0])

# решение многочлена
a = Symbol('a')
p = solve(-13 * a ** 4 + a ** 3 * (30 + 6 * b) - a ** 2 * (21 + 6 * b + b ** 2 - I30 * 3 / y) + a * (
    4 + I30 * b / y + I30 * 3 / y) - I30 ** 2 / 4 * y ** 2, a)
for x in p:
    print(x)


#print("ppppppppppppppppppppp ", type(p[0]))

print("I30 = ", I30)
# перевод корней в привычные типы - float и complex
p_Re = []
for rt in p:
    if (abs(rt - abs(rt)) < eps) or (abs(rt + abs(rt)) < eps):
        p_Re.append(float(rt))
x1 = p_Re[0]
x2 = p_Re[1]

x3 = float(p[2] + p[3]) / 2
x4 = float((p[3] - p[2]) / (2j))
print("x1 = ", x1)
print("x2 = ", x2)
print("x3 = ", x3)
print("x4 = ", x4)
c = math.sqrt((x2 - x3) ** 2 + x4 ** 2)
d = math.sqrt((x1 - x3) ** 2 + x4 ** 2)
#m = (c*d+(x3-x1)*(x2-x3)-x4**2)/(2*c*d)
m = math.sqrt(((d + c) ** 2 - (x2 ** 2 + x1 ** 2)) / c * d) * 1 / 2


for i in range(0, len(z_values)):
    cn = ellipfun('cn', -math.sqrt(13 * c * d) * y * z_values[i], m**2)
    cn_array.append(cn)
    p3 = (c * x2 + d * x1 + cn * (c * x2 + d * x1)) / (cn * (c - d) + c + d)
    #p3 = (c * x2 + d * x1 + cn * (-c * x2 + d * x1)) / (-cn * (c - d) + c + d)
    #print("p3 = ", p3)
    p3_array.append(p3)



plt.figure(figsize=(8, 6))
plt.plot(z_values, abs(np.array(A1_values)**2), color='r', linestyle='-')
plt.plot(z_values, abs(np.array(A3_values)**2), color='b', linestyle='-', linewidth= 3)
#plt.plot(z_all, abs(np.array(A_all)), color='m', linestyle='-')
plt.plot(z_values, abs(np.array(p3_array)), color='g', linestyle='-')
plt.xlabel('X')
plt.ylabel('Y')
#plt.title('График третьего инварианта')
plt.title('График решений: A1^2 - red, A3^2 - blue, p3 - green')
plt.grid(True)
plt.show()