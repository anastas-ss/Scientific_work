import math

from sympy.solvers import solve
from mpmath import *
from sympy import Symbol

# Задаем параметры
Lz = 1
h = 10**(-3)
i = complex(0, 1)
A1_next_n = A1_next = 1
A3_next_n = A3_next = 0
z = 0.0
y = 1
delta_k = 1
b = delta_k/y
eps = 10**(-9)

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
print("z = ", z)
print("A1_n_n = ", A1_n_n)
print("A3_n_n = ", A3_n_n)
#нахождение 1ого инварианта
I1 = abs(A1_n_n)**2 + abs(A3_n_n)**2
#нахождение 3ого инварианта
I3 = 3/2*y*(
        abs(A1_n_n)**4 + 4*abs(A1_n_n)**2*abs(A3_n_n)**2 + abs(A3_n_n)**4) + y*(
        A3_n_n*(A1_n_n.conjugate())**3 + (A3_n_n.conjugate())*A1_n_n**3) - (3*abs(A1_n_n)**2 + abs(A3_n_n)**2)*delta_k/2
print("I1 = ", I1)
print("I3 = ", I3)
I30 = 2*I3 - 3*y*I1**2 + 3*delta_k*I1
print("I30 = ", I30)

#решение многочлена
a = Symbol('a')
p = solve(-13*a**4 + a**3*(30 + 6*b) - a**2*(21 + 6*b + b**2 - I30*3/y) + a*(4 + I30*b/y + I30*3/y) - I30**2/4*y**2, a)

#перевод корней в привычные типы - float и complex
p_Re = []
for rt in p:
    if (abs(rt-abs(rt)) < eps) or (abs(rt+abs(rt)) < eps):
        p_Re.append(float(rt))
x1 = p_Re[0]
x2 = p_Re[1]
x3 = float(p[2]+p[3])/2
x4 = float((p[3]-p[2])/(2j))

c = math.sqrt((x2 - x3)**2 + x4**2)
d = math.sqrt((x1 - x3)**2 + x4**2)
m = math.sqrt(((d + c)**2 - (x2**2 + x1**2))/c*d)*1/2
cn = ellipfun('cn', -math.sqrt(13*c*d)*y*z, m)
p3 = (c*x2 + d*x1 + cn*(-c*x2 + d*x1))/(-cn*(c - d) + c + d)

print("Корни:")
print("x1 = ", x1)
print("x2 = ", x2)
print("x3 = ", x3)
print("x4 = ", x4)

print("c = ", c)
print("d = ", d)
print("m = ", m)
print("cn = ", cn)
print("p3 = ", p3)







