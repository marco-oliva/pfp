#!/usr/bin/env python

# Imports
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob, os, itertools

# Set
base_folder = 'experiments_2_22_2021'
lengths = ['64', '128', '256', '512', '1024']


# Read csv files
csv_files = glob.glob('{}/*.csv'.format(base_folder))
df = pd.concat([pd.read_csv(fp, index_col=False).assign(length=os.path.basename(fp).split('_')[4].split('.')[0][1:]) for fp in csv_files])

w_p = ['{},{}'.format(row['w'], row['p']) for index,row in df.iterrows()]
df['w,p'] = w_p

measures = ['parse_lenght', 'dict_phrases', 'dict_tot_length']
measures_names = {'parse_lenght' : 'Parse Length', 'dict_phrases':'Dict Phrases', 'dict_tot_length':'Dict Tot Length'}

fig, axs = plt.subplots(ncols=len(lengths), nrows=len(measures), figsize=(40, 24))
fig.subplots_adjust(top=0.8)
for r, measure in enumerate(measures):
    for c, l in enumerate(lengths):
        to_plot = df.loc[df['length'] == l].pivot(index='w,p', columns='f', values=measure)
        sns.heatmap(to_plot, ax=axs[r][c], linewidths=.5, cmap='YlGnBu')
        axs[r][c].set_title('Measure: {}, DS size: {}'.format(measures_names[measure], l))

fig.savefig('chromosome_19_results.png')




