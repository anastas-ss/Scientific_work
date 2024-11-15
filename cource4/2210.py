# нашли производную третьего инварианта (измененного):
# рассматриваем два случая: cos = 1 и cos = -1
#  приравниваем к 0 и находим корни

import sympy as sp
import numpy as np

from sympy import symbols

a3, b = symbols('a3 b')

equation1 = 12*a3 - 24*a3*a3*a3 + 4*1*(3*(1-a3**2)**(1/2)*a3*a3 + (1-a3*a3)**(3/2)) + 4*b*a3
equation2 = 12*a3 - 24*a3*a3*a3 + 4*(-1)*(3*(1-a3**2)**(1/2)*a3*a3 + (1-a3*a3)**(3/2)) + 4*b*a3

rs1, rs2 = [], []

b_values = np.arange(-1, 1.1, 0.1)

#print("cos = 1:")
for b_val in b_values:
    solution1 = sp.solve(equation1.subs(b, b_val), a3)
    real_solutions1 = [sol.evalf() for sol in solution1 if sol.is_real]
    print((b_val, real_solutions1[0]))

#print("cos = -1:")
for b_val in b_values:
    solution2 = sp.solve(equation2.subs(b, b_val), a3)
    real_solutions2 = [sol.evalf() for sol in solution2 if sol.is_real]
    print((b_val, real_solutions2[0]))

