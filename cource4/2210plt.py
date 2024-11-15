# нашли производную третьего инварианта (измененного):
# рассматриваем два случая: cos = 1 и cos = -1
#  приравниваем к 0 и находим корни
# строим график полученных решений

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# Определяем переменные
a3, b = sp.symbols('a3 b')

# Определяем уравнения
equation1 = 12*a3 - 24*a3**3 + 4*(3*(1-a3**2)**(1/2)*a3**2 + (1-a3**2)**(3/2)) + 4*b*a3
equation2 = 12*a3 - 24*a3**3 + 4*(-1)*(3*(1-a3**2)**(1/2)*a3**2 + (1-a3**2)**(3/2)) + 4*b*a3

# Списки для хранения значений
b_values = np.arange(-1, 1.1, 0.1)
a3_values1 = []
a3_values2 = []

# Решение для cos = 1
for b_val in b_values:
    solution1 = sp.solve(equation1.subs(b, b_val), a3)
    real_solutions1 = [sol.evalf() for sol in solution1 if sol.is_real]
    if real_solutions1:  # Проверяем, есть ли решения
        a3_values1.append(real_solutions1[0])
    else:
        a3_values1.append(None)  # Если нет решения, добавляем None

# Решение для cos = -1
for b_val in b_values:
    solution2 = sp.solve(equation2.subs(b, b_val), a3)
    real_solutions2 = [sol.evalf() for sol in solution2 if sol.is_real]
    if real_solutions2:  # Проверяем, есть ли решения
        a3_values2.append(real_solutions2[0])
    else:
        a3_values2.append(None)  # Если нет решения, добавляем None

# Построение графика
plt.figure(figsize=(10, 6))
plt.plot(b_values, a3_values1, label='cos = 1', marker='o')
plt.plot(b_values, a3_values2, label='cos = -1', marker='x')
plt.title('Зависимость b от a3')
plt.xlabel('b')
plt.ylabel('a3')
plt.axhline(0, color='black', linewidth=0.5, linestyle='--')
plt.axvline(0, color='black', linewidth=0.5, linestyle='--')
plt.grid()
plt.legend()
plt.show()
