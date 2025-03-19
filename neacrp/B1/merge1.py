import numpy as np
import matplotlib.pyplot as plt

# Define file paths
#power_files1 = ['nea1x1.txt']
#power_files_1 = ['nea1x1x2.txt', 'nea1x1x4.txt']
power_files2 = ['nea2x2.txt', 'nea4x4.txt']
#reactivity_files1 = ['react1x4.txt']
reactivity_files2 = ['react2x2.txt', 'react4x4.txt']

#power_names1 = ['Radial Node 1x1']
#power_names_1 = ['Axial Node 1x2', 'Axial Node 1x4']
power_names2 = ['Radial Node 2x2', 'Radial Node 4x4']
#react_names1 = ['Axial Node 1x4']
react_names2 = ['Radial Node 2x2', 'Radial Node 4x4']

# Initialize dictionaries to store data
reactivity_data1 = {}
reactivity_data2 = {}
power_data1 = {}
power_data_1 = {}
power_data2 = {}


colors1 = ['black']
colors_1 = ['dimgray','black']
colors__1 = ['dimgray','black']
colors2 = ['blue']
colors_2 = ['blue','dodgerblue']



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
    0.52043961352657,
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
    2.384793388429752,
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




# Read power files
for idx, file_path in enumerate(power_files2):
    try:
        name = power_names2[idx]
        data = np.genfromtxt(file_path, skip_header=2)
        time = data[:, 0]  # Time column
        relative_power2 = data[:, 1]  # Relative power column
        power_data2[name] = (time, relative_power2)
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")

# Read reactivity files
for idx_, file_path_ in enumerate(reactivity_files2):
    try:
        name = react_names2[idx_]
        data = np.genfromtxt(file_path_, skip_header=2)
        time = data[:, 0]  # Time column
        reactivity2 = data[:, 1]  # Reactivity column
        reactivity_data2[name] = (time, reactivity2)
    except Exception as e:
        print(f"Error reading file {file_path_}: {e}")

# Plotting with dual y-axes
plt.figure(figsize=(12, 7))


# Plot relative power on the left y-axis
ax1 = plt.gca()  # Get the current axes
color_index_ = 0  # Initialize color index
for key, (time, rel_power1) in power_data1.items():
    ax1.plot(time, rel_power1, label=key, linestyle='--', marker='o',linewidth=1.2, markersize=1,
    color=colors1[color_index_ % len(colors1)])
    color_index_ += 1  # Update color index
color_index_3 = 0  # Initialize color index
for key, (time, rel_power_1) in power_data_1.items():
    ax1.plot(time, rel_power_1, label=key, linestyle='-.', marker='o',linewidth=1.2, markersize=1,
    color=colors__1[color_index_3 % len(colors__1)])
    color_index_3 += 1  # Update color index
color_index_2 = 0  # Initialize color index
for key, (time, rel_power2) in power_data2.items():
    ax1.plot(time, rel_power2, label=key, linestyle='-', marker='o',linewidth=1.2, markersize=1,
    color=colors_1[color_index_2 % len(colors_1)])
    color_index_2 += 1  # Update color index

ax1.set_xlabel('Time (s)', fontsize=14, fontweight='bold')
ax1.set_ylabel('Relative Power', fontsize=14, fontweight='bold', color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax1.grid(True, linestyle='--', alpha=0.7)
# Set the limits of the x and y axes to start from 0
ax1.set_xlim(left=0)  # x-axis starts from 0
ax1.set_ylim(bottom=0)  # Left y-axis starts from 0
ax1.set_ylim(top=3)  # Left y-axis starts from 0
ax1.set_xlim(right=5)  # x-axis starts from 0

# Create a second y-axis for reactivity
ax2 = ax1.twinx()
color_index = 0  # Initialize color index
for key, (time, reactivity1) in reactivity_data1.items():
    ax2.plot(time, reactivity1, label=key, linestyle='--', marker='x',linewidth=1.2, markersize=1, 
    color=colors2[color_index % len(colors2)])
    color_index += 1  # Update color index

color_index2 = 0  # Initialize color index
for key, (time, reactivity2) in reactivity_data2.items():
    ax2.plot(time, reactivity2, label=key, linestyle='-', marker='x',linewidth=1.2, markersize=1, 
    color=colors_2[color_index2 % len(colors_2)])
    color_index2 += 1  # Update color index


ax2.set_ylabel('Reactivity ($)', fontsize=14, fontweight='bold', color='blue')
ax2.tick_params(axis='y', labelcolor='blue')
# Set the right y-axis to start from 0
ax2.set_ylim(bottom=0)
ax2.set_ylim(top=1.2)  # Left y-axis starts from 0

# Set fontweight='bold' for both x and y ticks on ax1
for label in ax1.get_xticklabels() + ax1.get_yticklabels():
    label.set_fontweight('bold')

# Set fontweight='bold' for both x and y ticks on ax2
for label in ax2.get_yticklabels():
    label.set_fontweight('bold')


# Plot reference results of relative power
ax1.plot(time2, rel_power_ref2, color='black', marker='+', linestyle=':', linewidth=0.6, label='PANTHER', markersize=8)
ax1.plot(time3, rel_power_ref3, color='black', marker='+', linestyle=':', linewidth=0.6, markersize=8)
ax1.plot(time1, rel_power_ref1, color='black', marker='o', linestyle=':', linewidth=0.6, label='QUANDRY', markersize=6)
# Plot reference results of reactivity

# Add a combined legend
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc='upper right', fontsize=9)


# Title and save the plot
#plt.title('Relative Power and Reactivity Over Time', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig('dual_y_axes_plot.png', dpi=500)
plt.show()
