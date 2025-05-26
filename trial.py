# plot sin(kx) for k = 1, 2, 3, 4, 5
import numpy as np
import matplotlib.pyplot as plt
# using a for loop and not def function

for k in range(1, 3):
    x = np.linspace(-10, 10, 100)
    y = np.sin(k * x)
    plt.plot(x, y, label=f'k={k}')
    plt.title('Sin(kx) for k = 1, 2, 3, 4, 5')
    plt.xlabel('x')
    plt.ylabel('sin(kx)')
    plt.legend()
    plt.grid()
    plt.savefig(f'sin_kx_plot_{k}.png')
    plt.show()
