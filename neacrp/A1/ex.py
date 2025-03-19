import matplotlib.pyplot as plt
import numpy as np

# Original data files
file_path = 'nea4x4.txt'
file_path_ = 'react4x4.txt'

# Read the data files
data = np.genfromtxt(file_path, skip_header=2)
data_ = np.genfromtxt(file_path_, skip_header=2)

# Extract columns for time and values
time1 = data[:, 0]
rel_power_ref1 = data[:, 1]

time = data_[:, 0]
react = data_[:, 1]

# Provided reference data
time3 = [
    0.0182648401826484, 0.0365296803652968, 0.0578386605783866, 
    0.076103500761035, 0.1035007610350076, 0.2343987823439878, 
    0.3835616438356164, 0.4687975646879756, 0.5814307458143074, 
    0.6605783866057838, 0.7640791476407914, 0.8706240487062404, 
    1.0076103500761033, 1.1933028919330289, 1.4033485540334856, 
    1.5981735159817352, 1.7990867579908674, 1.9969558599695585
]
rel_power_ref3 = [
    0.008035714285714285,
    0.17142857142857143,
    0.36428571428571427,
    0.5544642857142857,
    0.7714285714285715,
    0.9857142857142858,
    1.1491071428571429,
    1.0232142857142856,
    0.8946428571428572,
    0.7178571428571429,
    0.58125,
    0.43125,
    0.3080357142857143,
    0.24107142857142858,
    0.2142857142857143,
    0.20089285714285715,
    0.1982142857142857,
    0.1982142857142857
]
time_ref_1 = [
    0.0547945205479452, 0.1948249619482496, 0.3013698630136986, 0.3744292237442922,
    0.4170471841704718, 0.4535768645357686, 0.471841704718417, 0.487062404870624,
    0.502283105022831, 0.5144596651445966, 0.5296803652968036, 0.5449010654490106,
    0.5692541856925418, 0.5905631659056316, 0.60882800608828, 0.6301369863013698,
    0.669710806697108, 0.7275494672754946, 0.8310502283105022, 1.0624048706240488,
    1.4977168949771689, 1.7716894977168949, 1.9634703196347032
]
rel_power_ref_1 = [
    0, 0, 0.005405405405405405, 0.04594594594594594, 0.14864864864864863,
    0.3648648648648648, 0.5702702702702702, 0.7513513513513512, 0.9324324324324323,
    1.1189189189189188, 1.2297297297297296, 1.2594594594594593, 1.1729729729729728,
    0.9594594594594593, 0.8243243243243242, 0.6756756756756757, 0.5135135135135135,
    0.40810810810810805, 0.3216216216216216, 0.24864864864864863, 0.2135135135135135,
    0.2081081081081081, 0.2054054054054054
]

# Create the figure and the first axis
fig, ax1 = plt.subplots(figsize=(10, 6))

# Plot relative power and its reference
ax1.plot(time1, rel_power_ref1, color='cyan', marker='*', linestyle='--',
         linewidth=0.6, markersize=8, label='PARCS - Relative Power')
ax1.plot(time3, rel_power_ref3, color='blue', marker='x', linestyle=':',
         linewidth=0.6, markersize=6, label='Reference 3 - Relative Power')

ax1.set_xlabel('Time (s)', fontsize=14)
ax1.set_ylabel('Relative Power', fontsize=14, color='cyan')
ax1.tick_params(axis='y', labelcolor='cyan')
ax1.legend(loc='upper left', fontsize=12)

# Create the second y-axis for reactivity
ax2 = ax1.twinx()
ax2.plot(time, react, color='magenta', marker='o', linestyle='-',
         linewidth=0.6, markersize=6, label='Reactivity')
ax2.plot(time_ref_1, rel_power_ref_1, color='red', marker='s', linestyle=':',
         linewidth=0.6, markersize=6, label='Reference 1 - Reactivity')

ax2.set_ylabel('Reactivity ($)', fontsize=14, color='magenta')
ax2.tick_params(axis='y', labelcolor='magenta')
ax2.legend(loc='upper right', fontsize=12)

# Add grid and title
ax1.grid(True, linestyle='--', linewidth=0.7)
plt.title('Dual Y-Axis Plot with Reference Data: Relative Power and Reactivity',
          fontsize=16, fontweight='bold')

# Save or show the plot
plt.savefig('dual_y_axis_plot_with_provided_references.png', dpi=300)  # Save the plot
plt.show()  # Display the plot
