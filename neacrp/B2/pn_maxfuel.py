import numpy as np
import matplotlib.pyplot as plt

# File path to your data file
file_path = 'maxf4x4.txt'
file_path2 = 'maxf2x2.txt'
file_path3 = 'maxf1x1.txt'
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

time_ref3 = [
    0.22693997071742314,0.3879941434846267,
    0.4758418740849195,
    0.5051244509516838,
    0.5710102489019033,
    0.746705710102489,
    1.0322108345534406,
    1.4494875549048316,
    2.4158125915080526,
    2.957540263543192,
    3.49194729136164,
    4.062957540263543,
    4.773060029282577,4.992679355783309
]
rel_power_ref3 = [
    284.47081784386614,287.2589219330855,
    303.9875464684015,
    330.939219330855,
    381.12509293680296,
    446.18085501858735,
    498.2254646840148,
    549.3407063197026,
    628.3369888475836,
    657.1473977695167,
    678.5228624535316,
    698.0395910780669,
    717.5563197026022,723.1325278810409
]


time_ref1 = [
    0.1390922401171303,
    0.3733528550512445,
    0.5124450951683748,
    0.5710102489019033,
    0.7320644216691069,
    1.2884333821376281,
    1.896046852122987,
    2.8184480234260616,
    3.7042459736456808,
    4.136163982430454,
    4.882869692532943,4.985358711566618
]
rel_power_ref1 = [
    284.47081784386614,
    285.40018587360595,
    323.5042750929368,
    367.1845724907063,
    429.45223048327136,
    513.0953531598512,
    569.7868029739777,
    622.760780669145,
    653.4299256505576,
    666.4410780669145,
    683.1697026022305,685.0284386617101
]

time_ref2 = [0.022779043280182234,0.38724373576309795,
    0.5694760820045558,
    0.6719817767653758,
    0.7858769931662871,
    1.0820045558086562,
    1.5148063781321186,
    2.5170842824601367,
    3.029612756264237,
    4.009111617312073,
    4.4874715261959,
    4.9886104783599095
]
rel_power_ref2 = [285.8569051580699,287.52079866888516,
    300.83194675540767,
    350.7487520798669,
    391.5141430948419,
    453.0782029950083,
    503.82695507487523,
    580.3660565723794,
    606.1564059900167,
    640.2662229617305,
    653.5773710482529,
    664.3926788685524
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
#plt.plot(time_ref, rel_power_ref, color='cyan', marker='*', linestyle='None', linewidth=1.0, markersize=8, label='PARCS')
#plt.plot(time_ref1, rel_power_ref1, color='crimson', marker='*', linestyle='None', linewidth=1.0, markersize=8, label='CORCA-K')

plt.plot(time_ref1, rel_power_ref1, color='cyan', marker='*', linestyle='--', linewidth=1.0, markersize=8, label='PARCS')
plt.plot(time_ref3, rel_power_ref3, color='crimson', marker='*', linestyle='--', linewidth=1.0, markersize=8, label='CORCA-K')
plt.plot(time_ref2, rel_power_ref2, color='teal', marker='*', linestyle='--', linewidth=1.0, markersize=8, label='PANTHER')

# Enhancing the plot
plt.title('', fontsize=18, fontweight='bold', color='navy')
plt.xlabel('Time (s)', fontsize=14, fontweight='bold')
plt.ylabel('Maximum Fuel Temperature (°C)', fontsize=14, fontweight='bold')

# Adding grid, legend, and fine-tuning axis
plt.grid(True, linestyle='--', linewidth=0.6, alpha=0.7)
#plt.xticks(np.arange(0, 5.5, 0.5), fontsize=12)  # Time axis ticks every 5 seconds
plt.yticks(np.arange(300, 755, 50), fontsize=12)  # Power axis ticks every 0.2
plt.ylim(280, 750)  # Limit y-axis range
plt.legend(fontsize=12, loc='lower right')
plt.tight_layout()

# Saving the enhanced plot
plt.savefig('enhanced_relative_power_plot.png', dpi=500)
plt.show()


