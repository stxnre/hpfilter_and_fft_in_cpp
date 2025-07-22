import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()


data = pd.read_csv("fft_hpfilter/season_trend_decomp.csv", header=0)

data['Quarter'] = data['Quarter'].str.replace(r'(\d+) (Q\d)', r'\1-\2', regex=True)
data['Quarter'] = pd.PeriodIndex(data['Quarter'], freq='Q').to_timestamp()

datamelt = data.melt('Quarter', var_name='cols', value_name='vals')

fig,axs = plt.subplots(nrows=3)
fig.suptitle("Season Trend Decomposition of Quarterly Australian Electricity Production")

sns.lineplot(data, ax=axs[0], x="Quarter", y = "Original")
sns.lineplot(data, ax=axs[1], x="Quarter", y = "Trend")
sns.lineplot(data, ax=axs[2], x="Quarter", y = "Seasonal")
plt.show()

plt.savefig("Season-Trend Decomposition",format="svg")