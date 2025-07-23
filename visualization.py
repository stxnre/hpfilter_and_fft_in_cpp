import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()


data = pd.read_csv("fft_hpfilter/season_trend_decomp.csv", header=0)

data['Quarter'] = data['Quarter'].str.replace(r'(\d+) (Q\d)', r'\1-\2', regex=True)
data['Quarter'] = pd.PeriodIndex(data['Quarter'], freq='Q').to_timestamp()


# Seasonal-Trend Decomposition
fig,axs = plt.subplots(nrows=3, figsize=(12,8))
fig.suptitle("Quarterly Beer Production in Australia")

sns.lineplot(data, ax=axs[0], x="Quarter", y = "Original")
axs[0].set_title("Original Series")
sns.lineplot(data, ax=axs[1], x="Quarter", y = "Trend",color="red")
axs[1].set_title("Trend")
sns.lineplot(data, ax=axs[2], x="Quarter", y = "Seasonal", color="green")
axs[2].set_title("Seasonal")
plt.tight_layout()
plt.savefig("australian_beer_decomp.svg",format="svg")
plt.close()


# Periodogram
period = pd.read_csv("fft_hpfilter/periodogram.csv",names=["PSD"],header=None)
period['Freq'] = np.arange(0,period.shape[0],1) / period.shape[0]

sns.lineplot(period, x = 'Freq', y = "PSD")
plt.tight_layout()
plt.savefig("aus_prod_periodogram.svg",format="svg")
plt.close()