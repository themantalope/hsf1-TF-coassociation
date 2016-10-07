import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt

data_dir = os.path.join("..", "data_files")
oscillations_file = os.path.join(data_dir, "hsf1_reg_oscillatory_genes.xlsx")

oscillations = pd.read_excel(open(oscillations_file))

osc_norm = oscillations[[c for c in oscillations.columns.tolist() if "Normalized" in c or "Symbol" in c]]
osc_norm.set_index("Symbol", inplace=True)

# print osc_norm
# print osc_norm.index
osc_norm = osc_norm.T
osc_norm["time"] = osc_norm.index
gammas = sns.load_dataset("gammas")

# print gammas


# print osc_norm

# sns.tsplot(data=osc_norm.values, time=osc_norm.index, condition=osc_norm.columns)

for i in range(0, len(osc_norm.columns), 5):
    osc_norm.iloc[:, i:i+5].plot()
    plt.show()

# osc_norm.iloc[:, 0:10].plot()
# plt.show()