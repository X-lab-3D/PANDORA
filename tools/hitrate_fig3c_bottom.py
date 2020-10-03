import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)

axis_font = {'fontname':'Arial', 'size':'14'}

#######################################################
# Plot the hitrate of fig 3c bottom
#######################################################

def hitrate(score, target):
    '''Compute hitrate.'''
    idx = np.argsort(score)
    hit = np.cumsum(target[idx])
    if hit.max()!= 0:
        return hit/hit.max()
    else:
        return None

# read the data
fname = '../../output/merged/all/all.rawdata.tsv'
df = pd.read_csv(fname, sep='\t')
df_test = df[df['label']=='Test']

# case names
cases = np.unique(df_test['caseID'].to_numpy())

#read all the data
hitrate_dr = []
hitrate_hs = []
for c in cases:

    print(c)
    df_case = df_test[df_test['caseID']==c]

    target = df_case['target'].to_numpy()
    dr = df_case['DR'].to_numpy()
    hs = df_case['HS'].to_numpy()

    dr = hitrate(dr,target)
    hs = hitrate(hs,target)
    if dr is not None:
        hitrate_dr.append(dr)
        hitrate_hs.append(hs)

# cut min length
len_max = np.array([len(x) for x in hitrate_dr]).min()
data_hs = np.array([x[:len_max] for x in hitrate_hs])
data_dr = np.array([x[:len_max] for x in hitrate_dr])

# hs data
hs_mean = np.quantile(data_hs, 0.5, axis=0)
hs_a = np.quantile(data_hs,0.25, axis=0)
hs_b = np.quantile(data_hs,0.75, axis=0)

# dr data
dr_mean = np.quantile(data_dr, 0.5, axis=0)
dr_a = np.quantile(data_dr,0.25, axis=0)
dr_b = np.quantile(data_dr,0.75, axis=0)

# xaxis
x = np.arange(len_max)

# figure
fig = plt.figure(num=None, figsize=(5, 3), dpi=80, facecolor='w', edgecolor='k')

N = 1000

# main plot
ax1  = fig.add_subplot()

# hs data
ax1.fill_between(x[:N],hs_a[:N], hs_b[:N], color='#3279a8', alpha=0.5)
ax1.plot(x[:N], hs_mean[:N], color='#3279a8', label='Haddock')

#dr data
ax1.fill_between(x[:N],dr_a[:N], dr_b[:N], color='#e88700', alpha=0.5)
ax1.plot(x[:N], dr_mean[:N], color='#e88700', label='DeepRank')

# legend
ax1.set_xlabel('Top N', **axis_font)
ax1.set_ylabel('Hitrate (All Models)', **axis_font)


# text
# ax1.text(650,0.86,'DeepRank',color = '#e88700', fontsize=14)
# ax1.text(1500,0.75,'Haddock',color = '#3279a8', fontsize=14)
ax1.text(700,0.84,'DeepRank',color = '#e88700', fontsize=14)
ax1.text(750,0.70,'Haddock',color = '#3279a8', fontsize=14)

# inset
ax2 = plt.axes([0,0,1,1])
ip = InsetPosition(ax1, [0.5,0.1,0.4,0.4])
ax2.set_axes_locator(ip)
mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec='0.5')

n = 100

# hs data
ax2.plot(x[:n], hs_mean[:n], color='#3279a8')
ax2.fill_between(x[:n],hs_a[:n], hs_b[:n], color='#3279a8', alpha=0.25)

# dr data
ax2.plot(x[:n], dr_mean[:n], color='#e88700')
ax2.fill_between(x[:n],dr_a[:n], dr_b[:n], color='#e88700', alpha=0.25)
ax2.set_yticks([0,0.1,0.2])
ax2.set_xticks([0,50,100])

# plot
plt.subplots_adjust(bottom=0.15)
plt.show()




# plt.plot(hitrate_dr, label='DeepRank' )
# plt.plot(hitrate_hs, label='Haddock' )
# plt.grid()
# # plt.xlim((0,1000))
# plt.xlabel('Top N', **axis_font)
# plt.ylabel('Hit Rate', **axis_font)
# plt.legend()
# plt.show()




