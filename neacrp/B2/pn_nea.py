import numpy as np
import matplotlib.pyplot as plt

# Relative power file
#file_path = 'nea1x1.txt'
file_path2 = 'nea2x2.txt'
file_path3 = 'nea4x4.txt'
#file_path4 = 'nea1x1x2.txt'
#file_path5 = 'nea1x1x4.txt'
#file_path3 = 'nea1x2.txt'

# Read the data file using numpy (skipping the first two header lines)
#data = np.genfromtxt(file_path, skip_header=2)
data2 = np.genfromtxt(file_path2, skip_header=2)
data3 = np.genfromtxt(file_path3, skip_header=2)
#data4 = np.genfromtxt(file_path4, skip_header=2)
#data5 = np.genfromtxt(file_path5, skip_header=2)

# Extract columns for time and relative power
time = data2[:, 0]              # First column: time
#rel_power = data[:, 1]         # Second column: relative power
rel_power2 = data2[:, 1]         # Second column: relative power
rel_power3 = data3[:, 1]         # Second column: relative power
#rel_power4 = data4[:, 1]         # Second column: relative power
#rel_power5 = data5[:, 1]         # Second column: relative power

# Second set of data for comparison
#time2 = [0.5, 1, 1.5, 2, 2.5]
#rel_power_ref = [0.9997826086956519, 0.2541784302653869, 0.2012718802936193, 
#0.19472896668548834,0.19420783645655876]


time2 = [
    0.2713178294573643,
    0.35193236714975845,
    0.3896135265700483,
    0.41086956521739126,
    0.42826086956521736,
    0.44565217391304346,
    0.46207729468599035,
    0.4717391304347826,
    0.4823671497584541,
    0.4900966183574879,
    0.5026570048309178,
    0.5171497584541063,
    0.53743961352657,
    0.5442028985507246,
    0.5548309178743961,
    0.5664251207729468,
    0.5818840579710145,
    0.6021739130434782,
    0.6166666666666667,
    0.6495169082125605,
    0.6833333333333333
]
rel_power_ref2 = [
    0.015789473684210527,
    0.06942148760330578,
    0.19834710743801653,
    0.40661157024793393,
    0.6545454545454545,
    1.006611570247934,
    1.4578512396694217,
    1.775206611570248,
    2.0677685950413225,
    2.266115702479339,
    2.5338842975206615,
    2.63801652892562,
    2.424793388429752,
    2.300826446280992,
    2.0429752066115703,
    1.7652892561983473,
    1.4528925619834712,
    1.1553719008264463,
    1.0165289256198349,
    0.8280991735537191,
    0.7289256198347108
]
time3 = [
    0.7235142118863048,
    0.8010335917312661,
    0.934031007751938,
    1.124031007751938,
    1.6666666666666665,
    2.235142118863049,
    2.9586563307493536,
    3.4237726098191215,
    4.0310077519379846,
    4.521963824289405,
    4.961240310077519
]
rel_power_ref3 = [
    0.6578947368421052,
    0.5368421052631579,
    0.43005263157894735,
    0.38505263157894735,
    0.3600578947368421,
    0.3368421052631579,
    0.3368421052631579,
    0.3368421052631579,
    0.33157894736842103,
    0.33157894736842103,
    0.33157894736842103
]


time1 = [
    0.34883720930232553,
    0.38759689922480617,
    0.40051679586563305,
    0.42051679586563305,
    0.43051679586563305,
    0.4434366925064599,
    0.450051679586563305,
    0.4563565891472868,
    0.4763565891472868,
    0.4892764857881137,
    0.5180361757105943,
    0.5380361757105943,
    0.5467958656330749,
    0.5813953488372092,
    0.5813953488372092,
    0.6201550387596899,
    0.7235142118863048,
    0.9173126614987079,
    1.214470284237726,
    1.8863049095607234,
    2.622739018087855,
    3.229974160206718,
    3.7209302325581395,
    4.315245478036175,
    4.8966408268733845
]
rel_power_ref1 = [
    0.09473684210526316,
    0.3631578947368421,
    0.6210526315789473,
    0.8789473684210526,
    1.1052631578947367,
    1.3473684210526315,
    1.694736842105263,
    1.9263157894736842,
    2.263157894736842,
    2.484210526315789,
    2.3105263157894735,
    2.057894736842105,
    1.7578947368421052,
    1.431578947368421,
    1.2263157894736842,
    0.9736842105263158,
    0.6842105263157895,
    0.4423157894736842,
    0.3905263157894737,
    0.35789473684210527,
    0.33157894736842103,
    0.33157894736842103,
    0.33157894736842103,
    0.33157894736842103,
    0.33157894736842103
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
#plt.plot(time[:len(rel_power)], rel_power, color='blue', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Radial Node 1x1')
plt.plot(time[:len(rel_power2)], rel_power2, color='red', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Radial Node 2x2')
plt.plot(time[:len(rel_power3)], rel_power3, color='green', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Radial Node 4x4')
#plt.plot(time[:len(rel_power4)], rel_power4, color='darkgoldenrod', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Axial Node 1x2')
#plt.plot(time[:len(rel_power5)], rel_power5, color='magenta', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Axial Node 1x4')
plt.xlim(0, max(time))  # Start x-axis at 0

# Plot reference results
#plt.plot(time3, rel_power_ref3, color='crimson', marker='*', linestyle='--', linewidth=0.6, markersize=8, label='CORCA-K')
plt.plot(time2, rel_power_ref2, color='teal', marker='*', linestyle='', linewidth=0.6, markersize=8, label='PANTHER')
plt.plot(time3, rel_power_ref3, color='teal', marker='*', linestyle='', linewidth=0.6, markersize=8)
plt.plot(time1, rel_power_ref1, color='cyan', marker='*', linestyle='', linewidth=0.6, markersize=8, label='QUANDRY')
# Plot Time vs. rel_power2 (second dataset)
#plt.plot(result[:len(new_rel_power)], new_rel_power, color='teal', marker='*', linestyle='--', linewidth=1.0, markersize=8, label='PANTHER')


# Enhancing the plot
#plt.title('', fontsize=18, fontweight='bold', color='navy')
plt.xlabel('Time (s)', fontsize=14, fontweight='bold')
plt.ylabel('Relative Power', fontsize=14, fontweight='bold')

# Adding grid, legend, and fine-tuning axis
plt.grid(True, linestyle='--', linewidth=0.6, alpha=0.7)
plt.xticks(np.arange(0, 5.5, 0.5), fontsize=12)  # Time axis ticks every 5 seconds
plt.yticks(np.arange(0, 2.9, 0.4), fontsize=12)  # Power axis ticks every 0.2
plt.ylim(-0.005, 2.8)  # Limit y-axis range
plt.legend(fontsize=12, loc='upper right')

# Saving the enhanced plot
plt.savefig('enhanced_relative_power_plot.png', dpi=10)
plt.show()


