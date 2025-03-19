import numpy as np
import matplotlib.pyplot as plt

# File path to your data file
file_path = 'nea1x1.txt'
file_path2 = 'nea2x2.txt'
file_path3 = 'nea4x4.txt'
# file_path3 = 'nea1x2.txt'

# Read the data file using numpy (skipping the first two header lines)
data = np.genfromtxt(file_path, skip_header=2)
data2 = np.genfromtxt(file_path2, skip_header=2)
data3 = np.genfromtxt(file_path3, skip_header=2)

# Extract columns for time and relative power
time = data[:, 0]              # First column: time
rel_power = data[:, 1]         # Second column: relative power
rel_power2 = data2[:, 1]       # Second column: relative power
rel_power3 = data3[:, 1]       # Second column: relative power

# Second set of data for comparison
time2 = [0.5, 1, 1.5, 2, 2.5]
rel_power_ref = [0.9997826086956519, 0.2541784302653869, 0.2012718802936193, 
0.19472896668548834, 0.18545454545454554]

# Generate the full time list from 0 to 60.0 in increments of 0.25
result = []  # List for new time points
new_rel_power = []  # List for new relative power values
rel_power_index = 0  # Index to track the corresponding rel_power2 values

for t in np.arange(0, 2.5, 0.5):
    if t in time2:
        result.append(t)
        new_rel_power.append(rel_power_ref[rel_power_index])
        rel_power_index += 1
    else:
        result.append(None)
        new_rel_power.append(None)

# Plot for Relative Power
plt.figure(figsize=(12, 7))

# Plot Time vs. rel_power1
plt.plot(time[:len(rel_power)], rel_power, color='blue', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='OpenNode (1x1x1 mesh size)')
plt.plot(time[:len(rel_power2)], rel_power2, color='red', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='OpenNode (2x2x1 mesh size)')
plt.plot(time[:len(rel_power3)], rel_power3, color='green', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='OpenNode (4x4x1 mesh size)')

# Plot Time vs. rel_power2 (second dataset)
# plt.plot(result[:len(new_rel_power)], new_rel_power, color='cyan', marker='*', linestyle='None', linewidth=1.0, markersize=8, label='SPANDEX')

# Enhancing the plot
plt.title('Time vs. Max. Fuel Centerline Temperature Comparison', fontsize=18, fontweight='bold', color='navy')
plt.xlabel('Time (s)', fontsize=14, fontweight='bold')
plt.ylabel('Max. Fuel Centerline Temperature', fontsize=14, fontweight='bold')

# Adding grid, legend, and fine-tuning axis
plt.grid(True, linestyle='--', linewidth=0.6, alpha=0.7)
plt.xticks(np.arange(0, 5.5, 0.5), fontsize=12)  # Time axis ticks every 5 seconds
plt.yticks(np.arange(1640, 1710, 5), fontsize=12)  # Power axis ticks every 0.2
plt.legend(fontsize=12, loc='upper right')

# Set x-axis limits to start from 0 without extra space
plt.xlim(0, 5)

# Saving the enhanced plot
plt.savefig('enhanced_relative_power_plot.png', dpi=300)
plt.show()
