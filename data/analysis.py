import pandas as pd
import matplotlib.pyplot as plt

# Import the data
file_path = 'pump_probe_cd_t2.csv'
data = pd.read_csv(file_path)
data = data.sort_values(by='t2')


# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(data['t2'], data['int1'], label='HO 1', marker='o')
plt.plot(data['t2'], data['int5'], label='HO 5', marker='o')
plt.plot(data['t2'], data['int7'], label='HO 7', marker='o')

# Add labels, title, and legend
plt.xlabel(r'$t_2$', fontsize=14)
plt.ylabel('Dichroism', fontsize=14)
plt.title('Intensity vs t2', fontsize=16)
plt.legend(fontsize=12)

plt.grid(True)

# Show the plot
plt.savefig('../plots/pump_probe_cd_t2.png', dpi=300)
plt.show()