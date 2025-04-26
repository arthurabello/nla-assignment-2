import numpy as np
import matplotlib.pyplot as plt

x = np.array([1, 2, 3])
y = np.array([1, 2, 2])

m, b = np.polyfit(x, y, 1)

x_line = np.linspace(x.min() - 0.5, x.max() + 0.5, 100)
y_line = m * x_line + b

plt.figure(figsize=(6, 4))
plt.scatter(x, y, label='Data')
plt.plot(x_line, y_line, label=f'Best line')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Least Squares Idea')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
