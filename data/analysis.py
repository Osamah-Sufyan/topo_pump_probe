import pandas as pd
import matplotlib.pyplot as plt

# Import the data
file_1 = 'current_A0_left_t2_0.csv'
file_2 = 'current_A0_left_t2_neg015.csv'
file_3 = 'current_A0_left_t2_01.csv'


data1 = pd.read_csv(file_1)
data1 = data1.sort_values(by='A0')
data2 = pd.read_csv(file_2)
data2 = data2.sort_values(by='A0')
data3 = pd.read_csv(file_3)
data3 = data3.sort_values(by='A0')

# delete duplicate A0 values
# now plot the data vs A0 in  a single plot
plt.figure(figsize=(10, 6))
plt.plot(data1['A0'], data1['current_y'], label='t2=0', color='blue')
plt.plot(data2['A0'], data2['current_y'], label='t2=0.1', color='orange')
plt.plot(data3['A0'], data3['current_y'], label='t2=-0.1', color='green')
plt.xlabel('A0 (a.u.)')
plt.ylabel('edge current (a.u.)')
plt.title('Current vs A0 for different t2 values (Left Polarization)')
plt.legend()
plt.grid()

# Show the plot
plt.savefig('../plots/current_A0_left.png', dpi=300)
plt.show()