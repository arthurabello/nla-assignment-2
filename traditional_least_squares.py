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


x = np.array([1, 2, 3])
y = np.array([1, 2, 2])

m, b = np.polyfit(x, y, 1)

y_pred = m * x + b
residuals = y - y_pred

x_line = np.linspace(x.min() - 0.5, x.max() + 0.5, 100)
y_line = m * x_line + b

plt.figure(figsize=(6, 4))
plt.scatter(x, y, label='Data')
plt.plot(x_line, y_line, label='Best line')

for xi, yi, ypi, err in zip(x, y, y_pred, residuals):
    plt.plot([xi, xi], [ypi, yi], color='red')
    mid_y = (yi + ypi) / 2
    plt.text(xi + 0.05, mid_y, f'error_{xi}', va='center')

plt.xlabel('x')
plt.ylabel('y')
plt.title('The errors to be minimized')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

