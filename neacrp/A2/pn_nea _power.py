import numpy as np
import matplotlib.pyplot as plt

# File path to your data file
file_path = 'pow1x1x1.txt'
file_path2 ='pow1x1x2.txt'
file_path3 ='pow2x2.txt'
file_path4 ='pow4x4.txt'

# Read the data file using numpy (skipping the first two header lines)
data = np.genfromtxt(file_path, skip_header=2)
data2 = np.genfromtxt(file_path2, skip_header=2)
data3 = np.genfromtxt(file_path3, skip_header=2)
data4 = np.genfromtxt(file_path4, skip_header=2)

# Extract columns for time and relative power
time = data[:, 0]              # First column: time
rel_power = data[:, 1]         # Second column: relative power
rel_power2 = data2[:, 1]       # Second column: relative power
rel_power3 = data3[:, 1]       # Second column: relative power
rel_power4 = data4[:, 1]       # Second column: relative power

# Second set of data for comparison
time2 = [0.5, 1, 1.5, 2, 2.5]
rel_power_ref = [0.9997826086956519, 0.2541784302653869, 0.2012718802936193, 
0.19472896668548834, 0.18545454545454554]

# Generate the full time list from 0 to 60.0 in increments of 0.25
result = []  # List for new time points
new_rel_power = []  # List for new relative power values
rel_power_index = 0  # Index to track the corresponding rel_power2 values

#for t in np.arange(0, 2.5, 0.5):
#    if t in time2:
#        result.append(t)
#        new_rel_power.append(rel_power_ref[rel_power_index])
#        rel_power_index += 1
#    else:
#        result.append(None)
#        new_rel_power.append(None)


time1 = [
    0.005847953216374269,
    0.017543859649122806,
    0.023391812865497075,
    0.04093567251461988,
    0.05555555555555555,
    0.09649122807017543,
    0.30409356725146197,
    0.47953216374269003,
    0.7280701754385964,
    0.9064327485380117,
    1.1374269005847952,
    1.3567251461988303,
    1.657894736842105,
    1.8859649122807016,
    1.9678362573099415
]
rel_power_ref1 = [
    1.0143122676579925,
    1.0282527881040893,
    1.0407063197026023,
    1.0650557620817844,
    1.0771375464684017,
    1.0860594795539034,
    1.079925650557621,
    1.0734200743494424,
    1.0665427509293681,
    1.0622676579925652,
    1.0581784386617101,
    1.054646840148699,
    1.050743494423792,
    1.0485130111524164,
    1.0475836431226766
]

time3 = [
    0.005847953216374269,
    0.011695906432748537,
    0.02046783625730994,
    0.038011695906432746,
    0.05555555555555555,
    0.07017543859649122,
    0.2807017543859649,
    0.391812865497076,
    0.5847953216374269,
    1.0350877192982455,
    1.2456140350877192,
    1.3976608187134503,
    1.6286549707602338,
    1.7719298245614035,
    1.9269005847953216,
    1.9883040935672514
]
rel_power_ref3 = [
    1.0120817843866172,
    1.024907063197026,
    1.0401486988847584,
    1.0605947955390336,
    1.0771375464684017,
    1.0842007434944239,
    1.0808550185873607,
    1.0763940520446098,
    1.0700743494423792,
    1.0598513011152417,
    1.0566914498141264,
    1.054460966542751,
    1.0514869888475837,
    1.05,
    1.0485130111524164,
    1.0477695167286245
]

time2 = [
    0.029082774049217,
    0.03467561521252796,
    0.039149888143176735,
    0.04586129753914989,
    0.05704697986577181,
    0.07046979865771812,
    0.07829977628635347,
    0.10626398210290827,
    0.15436241610738255,
    0.23937360178970918,
    0.3232662192393736,
    0.38814317673378074,
    0.42953020134228187,
    0.4686800894854586,
    0.4977628635346756
]
rel_power_ref2 = [
    1.04,
    1.0458064516129033,
    1.0506912442396314,
    1.056405529953917,
    1.0658986175115208,
    1.0719815668202766,
    1.076405529953917,
    1.0774193548387097,
    1.0767741935483872,
    1.074746543778802,
    1.0726267281105992,
    1.0708755760368665,
    1.0696774193548388,
    1.068479262672811,
    1.067741935483871
]

time4 = [
    0.599078341013825,
1.0023041474654377,
1.5207373271889402	,2.2350230414746544	,3.0069124423963136,3.5944700460829493,4.274193548387097,4.688940092165899	,4.965437788018433]
rel_power_ref4 = [
    1.0650292397660819,
1.0559844054580898,1.0499668615984407,1.043261208576998,1.0385185185185186,1.0369590643274855,1.035711500974659,1.034775828460039,1.034775828460039]


# Plot for Relative Power
plt.figure(figsize=(12, 7))

# Plot Time vs. rel_power1
plt.plot(time[:len(rel_power3)], rel_power3, color='red', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Radial Node 2x2')
plt.plot(time[:len(rel_power4)], rel_power4, color='green', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Radial Node 4x4')
plt.plot(time[:len(rel_power)], rel_power, color='blue', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Axial Node 1x1')
plt.plot(time[:len(rel_power2)], rel_power2, color='magenta', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Axial Node 1x2')

# Plot Time vs. rel_power2 (second dataset)
# plt.plot(result[:len(new_rel_power)], new_rel_power, color='cyan', marker='*', linestyle='None', linewidth=1.0, markersize=8, label='SPANDEX')
plt.plot(time1, rel_power_ref1, color='cyan', marker='*', linestyle='--', linewidth=0.6, markersize=8, label='PARCS')
plt.plot(time3, rel_power_ref3, color='crimson', marker='*', linestyle='--', linewidth=0.6, markersize=8, label='CORCA-K')
plt.plot(time2, rel_power_ref2, color='teal', marker='*', linestyle='--', linewidth=0.6, markersize=8, label='PANTHER')
plt.plot(time4, rel_power_ref4, color='teal', marker='*', linestyle='--', linewidth=0.6, markersize=8)

# Enhancing the plot
#plt.title('', fontsize=18, fontweight='bold', color='navy')
plt.xlabel('Time (s)', fontsize=14, fontweight='bold')
plt.ylabel('Relative Power', fontsize=14, fontweight='bold')

# Adding grid, legend, and fine-tuning axis
plt.grid(True, linestyle='--', linewidth=0.6, alpha=0.7)
plt.xticks(np.arange(0, 5.5, 0.5), fontsize=12)  # Time axis ticks every 5 seconds
plt.yticks(np.arange(0, 1.10, 0.01), fontsize=12)  # Power axis ticks every 0.2
plt.ylim(1, 1.09)  # Limit y-axis range
plt.legend(fontsize=12, loc='upper right')

# Set x-axis limits to start from 0 without extra space
plt.xlim(0, 5)

# Saving the enhanced plot
plt.savefig('enhanced_relative_power_plot.png', dpi=500)
plt.show()
