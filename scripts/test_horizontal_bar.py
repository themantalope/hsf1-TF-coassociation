import mattplots
import os
import pandas as pd
import numpy as np
import math

data_dir = os.path.join("..", "data_files")
data_file = os.path.join(data_dir, "hsf1_features_cctf.xlsx")

data = pd.read_excel(open(data_file))

plotrows = data.loc[["Fisher Exact Odds Ratio", "Fisher Exact P-value"], :]
# print plotrows

logv = 10000

ytick_labels = plotrows.columns.tolist()
xlabel_left = plotrows.index[0]
xlabel_right = "-log{n}".format(n=logv) + " " + plotrows.index[1]
leftdata = plotrows.iloc[0, :].values
rightdata = [math.log(x,logv) for x in plotrows.iloc[1, :].values]


fig, plt = mattplots.back_to_back_bar(leftdata[::-1],
                                      rightdata[::-1],
                                      xaxis_left_label=xlabel_left,
                                      xaxis_right_label=xlabel_right,
                                      yaxis_tick_labels=ytick_labels[::-1],
                                      n_ticks_left=2)

plt.show()