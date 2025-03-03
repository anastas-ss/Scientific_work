import sympy as sp
import numpy as np

b, s, p3 = sp.symbols('b s p3')
b_value = 1

polynom = -13*p3**4 + (30 + 6*b)*p3**3 - (21 + 6*b + b**2 + 3*s)*p3**2 + (4 + b*s + 3*s)*p3 - s**2/4
sturm_sequence = sp.polys.polytools.sturm(polynom, p3)

print("Последний многочлен системы Штурма:")
print(sturm_sequence[-1])

for i, poly in enumerate(sturm_sequence):
    print(f"Polynom {i}: {poly}")

print("\n")
# print(sturm_sequence[-1].factor()) #представить в виде произведения скобок
# print(type(sturm_sequence[-1].factor()))

#print("s=0.5, b=1: ", polynom.subs(s,  -0.431271619546922).subs(b, 1))
#p = sp.poly(polynom.subs(s,  -0.431271619546922).subs(b, 1), p3).nroots()

intervals = [(-10, -5), (-5, 0), (0, 5), (5, 10)]


for interval in intervals:
    lower, upper = interval
    sign_count_lower = sum(sp.sign(poly.subs(p3, lower)) for poly in sturm_sequence)
    sign_count_upper = sum(sp.sign(poly.subs(p3, upper)) for poly in sturm_sequence)

    num_roots = sign_count_lower - sign_count_upper
    print(f"Количество корней на промежутке {interval}: {num_roots}")

