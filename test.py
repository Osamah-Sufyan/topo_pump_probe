import numpy as np 

import matplotlib.pyplot as plt

x = np.linspace(0, 10, 1000)
y = np.sin(x)

subset_start = 0  # Start of the subset
subset_end = 2*np.pi    # End of the subset

# Select the subset of x and y
x_subset = x[(x >= subset_start) & (x <= subset_end)]
y_subset = y[(x >= subset_start) & (x <= subset_end)]

# Perform numerical integration over the subset
integral = np.trapz(y_subset, x_subset)
print(f"Integral over the subset [{subset_start}, {subset_end}]: {integral}")

# Plot the sine wave and highlight the subset
plt.plot(x, y, label="Sine Wave")
plt.fill_between(x_subset, y_subset, alpha=0.3, label="Integration Subset")
plt.title("Sine Wave with Integration Subset")
plt.legend()
plt.show()