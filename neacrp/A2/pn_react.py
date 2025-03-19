import numpy as np
import matplotlib.pyplot as plt

# File path to your data file
file_path = 'react4x4.txt'
file_path2 = 'react2x2.txt'
file_path3 = 'react1x1.txt'
#file_path4 = 'nea1x1x2.txt'
#file_path5 = 'nea1x1x4.txt'
#file_path3 = 'nea1x2.txt'

# Read the data file using numpy (skipping the first two header lines)
data = np.genfromtxt(file_path, skip_header=2)
data2 = np.genfromtxt(file_path2, skip_header=2)
data3 = np.genfromtxt(file_path3, skip_header=2)
#data4 = np.genfromtxt(file_path4, skip_header=2)
#data5 = np.genfromtxt(file_path5, skip_header=2)

# Extract columns for time and relative power
time = data[:, 0]              # First column: time
rel_power = data[:, 1]         # Second column: relative power
rel_power2 = data2[:, 1]         # Second column: relative power
rel_power3 = data3[:, 1]         # Second column: relative power
#rel_power4 = data4[:, 1]         # Second column: relative power
#rel_power5 = data5[:, 1]         # Second column: relative power

# Second set of data for comparison
time2 = [0.5, 1, 1.5, 2, 2.5]
#rel_power_ref = [0.9997826086956519, 0.2541784302653869, 0.2012718802936193, 
#0.19472896668548834,0.19420783645655876]

time_ref = [
    0.008771929824561403,
    0.014619883040935672,
    0.02631578947368421,
    0.049707602339181284,
    0.07894736842105263,
    0.195906432748538,
    0.347953216374269,
    0.5087719298245614,
    0.6286549707602339,
    0.8654970760233918,
    1.0584795321637426,
    1.4269005847953216,
    1.722222222222222,
    1.95906432748538
]
rel_power_ref = [
    0.014869888475836434,
    0.030483271375464686,
    0.044423791821561344,
    0.06282527881040893,
    0.0684014869888476,
    0.06301115241635688,
    0.05520446096654276,
    0.047955390334572495,
    0.04330855018587361,
    0.03624535315985131,
    0.031784386617100376,
    0.02546468401486989,
    0.021189591078066918,
    0.01858736059479554
]

time_ref1 = [
    0.008771929824561403,
    0.023391812865497075,
    0.03508771929824561,
    0.043859649122807015,
    0.13157894736842105,
    0.29239766081871343,
    0.347953216374269,
    0.5146198830409356,
    0.7368421052631579,
    0.956140350877193,
    1.1023391812865497,
    1.3216374269005846,
    1.4678362573099415,
    1.6228070175438596,
    1.7690058479532162,
    1.915204678362573,
    1.9824561403508771
]
rel_power_ref1 = [
    0.018401486988847585,
    0.03940520446096655,
    0.05223048327137547,
    0.060594795539033464,
    0.06635687732342008,
    0.05817843866171005,
    0.05520446096654276,
    0.04776951672862454,
    0.039591078066914503,
    0.03382899628252788,
    0.031040892193308554,
    0.027323420074349444,
    0.024907063197026024,
    0.022862453531598517,
    0.02100371747211896,
    0.01951672862453532,
    0.018773234200743498
]

# Generate the full time list from 0 to 60.0 in increments of 0.25
result = []  # List for new time points
new_rel_power = []  # List for new relative power values
rel_power_index = 0  # Index to track the corresponding rel_power2 values

#for t in np.arange(0, 3.0, 0.5):
#    if t in time2:
#        result.append(t)
#        new_rel_power.append(rel_power_ref[rel_power_index])
#        rel_power_index += 1
#    else:
#        result.append(None)
#        new_rel_power.append(None)

# Plot for Relative Power
plt.figure(figsize=(12, 7))

# Plot Time vs. rel_power1
plt.plot(time[:len(rel_power3)], rel_power3, color='green', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Radial Node 1x1')
plt.plot(time[:len(rel_power2)], rel_power2, color='red', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Radial Node 2x2')
plt.plot(time[:len(rel_power)], rel_power, color='blue', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Radial Node 4x4')
#plt.plot(time[:len(rel_power4)], rel_power4, color='yellow', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Axial Node 1x1x2')
#plt.plot(time[:len(rel_power5)], rel_power5, color='magenta', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Axial Node 1x1x4')
plt.xlim(0, max(time))  # Start x-axis at 0

# Plot Time vs. rel_power2 (second dataset)
#plt.plot(result[:len(new_rel_power)], new_rel_power, color='cyan', marker='*', linestyle='None', linewidth=1.0, markersize=8, label='PANTHER')

# Plot reference results
plt.plot(time_ref, rel_power_ref, color='cyan', marker='*', linestyle='None', linewidth=1.0, markersize=8, label='PARCS')
plt.plot(time_ref1, rel_power_ref1, color='crimson', marker='*', linestyle='None', linewidth=1.0, markersize=8, label='CORCA-K')

# Enhancing the plot
plt.title('', fontsize=18, fontweight='bold', color='navy')
plt.xlabel('Time (s)', fontsize=14, fontweight='bold')
plt.ylabel('Reactivity ($)', fontsize=14, fontweight='bold')

# Adding grid, legend, and fine-tuning axis
plt.grid(True, linestyle='--', linewidth=0.6, alpha=0.7)
plt.xticks(np.arange(0, 5.5, 0.5), fontsize=12)  # Time axis ticks every 5 seconds
plt.yticks(np.arange(0, 0.08, 0.01), fontsize=12)  # Power axis ticks every 0.2
plt.ylim(0, 0.07)  # Limit y-axis range
plt.legend(fontsize=12, loc='upper right')
plt.tight_layout()

# Saving the enhanced plot
plt.savefig('enhanced_relative_power_plot.png', dpi=500)
plt.show()


