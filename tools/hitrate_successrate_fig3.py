import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


axis_font = {'fontname':'Arial', 'size':'14'}

#######################################################
# Plot the hitrate of fig 3c top
#######################################################

def hitrate(score, target):
    '''Compute the hitrate.'''
    idx = np.argsort(score) # https://github.com/numpy/numpy/issues/8757
    hit = np.cumsum(target[idx])
    if hit.max()!= 0:
        return hit/hit.max()
    else:
        return None
    
def successrate(score, target): # high <= 1;  1 < mid <= 1.5 , 1.5 <  acc <= 2, 2 < near-acc <= 2.5
    '''Compute the hitrate.'''
    idx = np.argsort(score)
    hit = np.cumsum(target[idx])
    success = []
    for x in hit:
        if x == 0:
            success.append(0)
        else:
            success.append(1)
    success = np.array(success)
    return success #gotta sum all of them position by position in "vertical"

# read the data
#fname = './it0.rawdata.tsv'
fname = './hitrate_rawdata.tsv'
df = pd.read_csv(fname, sep='\t')
#df_test = df[df['label']=='Test']
df_test = df

# case names
cases = np.unique(df_test['caseID'].to_numpy())

# read all the data
hitrate_dr = []
hitrate_hs = []
successrate_dr = []
successrate_hs = []
for c in cases:

    df_case = df_test[df_test['caseID']==c]

    target = df_case['target'].to_numpy()
    dr = df_case['DR'].to_numpy() #vector with all dr scores for given target
    hs = df_case['HS'].to_numpy()

    dr_hr = hitrate(dr,target)
    hs_hr = hitrate(hs,target)
    
    dr_sr = successrate(dr, target)
    hs_sr = successrate(hs, target)
    
    if dr_hr is not None:
        hitrate_dr.append(dr_hr)
        hitrate_hs.append(hs_hr)
    if dr_sr is not None:
        successrate_dr.append(dr_sr)
        successrate_hs.append(hs_sr)

successrate_hs = sum(successrate_hs)/len(successrate_hs)
successrate_dr = sum(successrate_dr)/len(successrate_dr)

# cut data to min length
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

# x asis
x = np.arange(1, len_max +1)

# fig
fig = plt.figure(num=None, figsize=(5, 3), dpi=80, facecolor='w', edgecolor='k')
N = 20

# plot hs
plt.fill_between(x[:N],hs_a[:N], hs_b[:N], color='#3279a8', alpha=0.5)
plt.plot(x[:N], hs_mean[:N], color='#3279a8', label='Haddock')

# plot dr
plt.fill_between(x[:N],dr_a[:N], dr_b[:N], color='#e88700', alpha=0.5)
plt.plot(x[:N], dr_mean[:N], color='#e88700', label='DeepRank')

# label and stuff
plt.xlabel('Top N', **axis_font)
plt.ylabel('Hitrate', **axis_font)
plt.title('Hitrate 2A')

# text
# plt.text(400,0.7,'DeepRank',color = '#e88700', fontsize=14)
# plt.text(2100,0.05,'Haddock',color = '#3279a8', fontsize=14)

#plt.text(650,0.80,'DeepRank',color = '#e88700', fontsize=14)
#plt.text(720,0.13,'Haddock',color = '#3279a8', fontsize=14)

# plot
plt.subplots_adjust(bottom=0.15)
plt.savefig('./Hitrate.png')
plt.show()
plt.clf()

#%%

# cut data to min length
#len_max = np.array([len(x) for x in successrate_dr]).min()
#data_hs = np.array([x[:len_max] for x in successrate_hs])
#data_dr = np.array([x[:len_max] for x in successrate_dr])

# hs data
#hs_mean = np.quantile(data_hs, 0.5, axis=0)
#hs_a = np.quantile(data_hs,0.25, axis=0)
#hs_b = np.quantile(data_hs,0.75, axis=0)

# dr data
#dr_mean = np.quantile(data_dr, 0.5, axis=0)
#dr_a = np.quantile(data_dr,0.25, axis=0)
#dr_b = np.quantile(data_dr,0.75, axis=0)

# x asis
x = np.arange(1, len(successrate_dr) + 1)

# fig
fig = plt.figure(num=None, figsize=(5, 3), dpi=80, facecolor='w', edgecolor='k')
N = 10

# plot hs
#plt.fill_between(x[:N],hs_a[:N], hs_b[:N], color='#3279a8', alpha=0.5)
plt.plot(x[:N], successrate_hs[:N], color='#3279a8', label='Haddock')

# plot dr
#plt.fill_between(x[:N],dr_a[:N], dr_b[:N], color='#e88700', alpha=0.5)
plt.plot(x[:N], successrate_dr[:N], color='#e88700', label='DeepRank')

# label and stuff
plt.xlabel('Top N', **axis_font)
plt.ylabel('Successrate', **axis_font)
plt.title('Successrate 2A')

# text
# plt.text(400,0.7,'DeepRank',color = '#e88700', fontsize=14)
# plt.text(2100,0.05,'Haddock',color = '#3279a8', fontsize=14)

#plt.text(650,0.80,'DeepRank',color = '#e88700', fontsize=14)
#plt.text(720,0.13,'Haddock',color = '#3279a8', fontsize=14)

# plot
plt.subplots_adjust(bottom=0.15)
plt.savefig('./Successrate.png')
plt.show()
plt.clf()



