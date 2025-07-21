import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


data = pd.read_csv("fft_hpfilter/season_trend_decomp.csv", header=0)

data['Quarter'] = data['Quarter'].str.replace(r'(\d+) (Q\d)', r'\1-\2', regex=True)
data['Quarter'] = pd.PeriodIndex(data['Quarter'], freq='Q').to_timestamp()

datamelt = data.melt('Quarter', var_name='cols', value_name='vals')

sns.catplot(x="Quarter", y="vals", hue='cols', data=datamelt, kind='line')
plt.show()