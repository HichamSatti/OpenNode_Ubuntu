import numpy as np
import matplotlib.pyplot as plt

# Define file paths
power_files1 = ['pow1x1x1.txt']
power_files_1 = ['pow1x1x2.txt']
power_files2 = ['pow2x2.txt', 'pow4x4.txt']
reactivity_files1 = ['react1x1.txt']
reactivity_files2 = ['react2x2.txt', 'react4x4.txt']

power_names1 = ['Radial Node 1x1']
power_names_1 = ['Axial Node 1x2']
power_names2 = ['Radial Node 2x2', 'Radial Node 4x4']
react_names1 = ['Radial Node 1x1']
react_names2 = ['Radial Node 2x2', 'Radial Node 4x4']

# Initialize dictionaries to store data
reactivity_data1 = {}
reactivity_data2 = {}
power_data1 = {}
power_data_1 = {}
power_data2 = {}


colors1 = ['black']
colors__1 = ['dimgray']
colors_1 = ['dimgray','black']
colors2 = ['blue']
colors_2 = ['blue','dodgerblue']


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




time_ref_ = [
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
rel_power_ref_ = [
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
time_ref_1 = [
    0.008771929824561403,
    0.023391812865497075,
    0.03508771929824561,
    0.043859649122807015,
    0.083859649122807015,
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
rel_power_ref_1 = [
    0.018401486988847585,
    0.03940520446096655,
    0.05223048327137547,
    0.060594795539033464,
    0.067594795539033464,
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

# Read power files
for idx, file_path in enumerate(power_files1):
    try:
        name = power_names1[idx]
        data = np.genfromtxt(file_path, skip_header=2)
        time = data[:, 0]  # Time column
        relative_power1 = data[:, 1]  # Relative power column
        power_data1[name] = (time, relative_power1)
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
for idx, file_path in enumerate(power_files1):
    try:
        name = power_names_1[idx]
        data = np.genfromtxt(file_path, skip_header=2)
        time = data[:, 0]  # Time column
        relative_power_1 = data[:, 1]  # Relative power column
        power_data_1[name] = (time, relative_power_1)
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
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
for idx_, file_path_ in enumerate(reactivity_files1):
    try:
        name = react_names1[idx_]
        data = np.genfromtxt(file_path_, skip_header=2)
        time = data[:, 0]  # Time column
        reactivity1 = data[:, 1]  # Reactivity column
        reactivity_data1[name] = (time, reactivity1)
    except Exception as e:
        print(f"Error reading file {file_path_}: {e}")
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
#ax1.set_ylim(bottom=0)  # Left y-axis starts from 0
#ax1.ax1.set_ylim(top=1.4)  # Left y-axis starts from 0
ax1.set_xlim(right=5)  # x-axis starts from 0
ax1.set_ylim(1.0, 1.10)
# Set y-axis ticks (200 to 500 by step of 50)
yticks = np.arange(1, 1.11, 0.02)  # 501 is used to include 500
ax1.set_yticks(yticks)


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
#ax2.set_ylim(top=1.2)  # Left y-axis starts from 0

ax2.set_ylim(0, 0.10)
# Set y-axis ticks (200 to 500 by step of 50)
yticks = np.arange(0, 0.11, 0.02)  # 501 is used to include 500
ax2.set_yticks(yticks)

# Set fontweight='bold' for both x and y ticks on ax1
for label in ax1.get_xticklabels() + ax1.get_yticklabels():
    label.set_fontweight('bold')

# Set fontweight='bold' for both x and y ticks on ax2
for label in ax2.get_yticklabels():
    label.set_fontweight('bold')


# Plot reference results of relative power
ax1.plot(time1, rel_power_ref1, color='black', marker='+', linestyle=':', linewidth=0.6, label='PARCS', markersize=8)
ax1.plot(time3, rel_power_ref3, color='black', marker='*', linestyle=':', linewidth=0.6, label='CORCA-K', markersize=8)
ax1.plot(time2, rel_power_ref2, color='black', marker='o', linestyle=':', linewidth=0.6, label='PANTHER', markersize=6)
ax1.plot(time4, rel_power_ref4, color='black', marker='o', linestyle=':', linewidth=0.6, markersize=6)
# Plot reference results of reactivity
ax2.plot(time_ref_, rel_power_ref_, color='blue', marker='+', linestyle=':', linewidth=0.6, label='PARCS', markersize=8)
ax2.plot(time_ref_1, rel_power_ref_1, color='blue', marker='*', linestyle=':', linewidth=0.6, label='CORCA-K', markersize=8)

# Add a combined legend
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc='upper right', fontsize=9)


# Title and save the plot
#plt.title('Relative Power and Reactivity Over Time', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig('dual_y_axes_plot.png', dpi=500)
plt.show()
