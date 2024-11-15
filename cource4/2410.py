# рассматриваем случай f' = 0 при cos = 1, где f это I3 modified
# находим значения a3 при различных значениях b [-1, 1], оставляем из них только действительные значения
# добавляем к полученным значения в точках 0 и 1, ищем min и max

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# Определяем символы
a3, b = sp.symbols('a3 b')

# Определяем уравнение
equation1 = 12*a3 - 24*a3**3 + 4*(3*(1-a3**2)**(1/2)*a3**2 + (1-a3**2)**(3/2)) + 4*b*a3

# Список для хранения результатов
results = []
b_values = np.arange(-1, 1.1, 0.1)

# Решаем уравнение для различных значений b
for b_val in b_values:
    solution1 = sp.solve(equation1.subs(b, b_val), a3)
    real_solutions1 = [sol.evalf() for sol in solution1 if sol.is_real]
    
    if real_solutions1:
        a3_val = real_solutions1[0]  # Берем первое реальное решение
        
        # Подставляем в функцию f
        f = 6*a3_val**2 - 6*a3_val**4 + 4*a3_val*(1-a3_val**2)**(3/2) + 2*b_val*a3_val**2
        results.append((b_val, f))

# Значения функции при a3 = 0 и a3 = 1
f_a3_0 = [6*0**2 - 6*0**4 + 4*0*(1-0**2)**(3/2) + 2*b_val*0**2 for b_val in b_values]
f_a3_1 = [6*1**2 - 6*1**4 + 4*1*(1-1**2)**(3/2) + 2*b_val*1**2 for b_val in b_values]

# Объединяем результаты
all_b_values = np.concatenate((b_values, b_values, b_values))
all_f_values = np.concatenate((np.array([f for _, f in results]), f_a3_0, f_a3_1))

# Находим максимум и минимум
max_f = max(all_f_values)
min_f = min(all_f_values)

# Строим график
plt.figure(figsize=(10, 6))
plt.plot(b_values, [f for _, f in results], label='f(a3)', marker='o')
plt.plot(b_values, f_a3_0, label='f(a3=0)', linestyle='--')
plt.plot(b_values, f_a3_1, label='f(a3=1)', linestyle='--')

# Отмечаем максимум и минимум
plt.scatter([b_val for b_val, f in results if f == max_f], [max_f], color='red', zorder=5, label='Max')
plt.scatter([b_val for b_val, f in results if f == min_f], [min_f], color='blue', zorder=5, label='Min')

# Настройки графика
plt.title('Зависимость f от b')
plt.xlabel('b')
plt.ylabel('f')
plt.axhline(0, color='black',linewidth=0.5, ls='--')
plt.axvline(0, color='black',linewidth=0.5, ls='--')
plt.legend()
plt.grid()
plt.show()

print(f"Максимум: {max_f}, Минимум: {min_f}")
