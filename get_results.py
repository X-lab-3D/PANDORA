#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 19:04:24 2020

@author: rafaella
"""
###usage : run from pMHC_Modelling main folder 
# python get_results.py 
## the folders are hardcoded to make it easier to run
#outputs Z value and p value for each comparison 

#%%
import os
import sys
import pickle

#%%
### PANDORA
def pandora_rmsd (directory):
    os.chdir(directory)
    with open('all_best_RMSDs.pkl', 'rb') as infile:
        best_rmsds = pickle.load(infile)
    
    with open('../../benchmark/PANDORA_benchmark_dataset.tsv', 'r') as cvfile:
        for i, line in enumerate(cvfile):
            if i != 0:
                row = line.split('\t')
                try:
                    best_rmsds[row[0]].append(len(row[1]))
                except:
                    pass
    return best_rmsds


#%%
### APE-Gen
def apgen_rmsd(best_rmsds, directory):
    os.chdir(directory)
    bests = {}
    for fol in os.listdir('./'):
        if os.path.isdir('./' + fol):
            try:
                ID = fol[-4:]
               # length = best_rmsds[ID][-1]
                rmsds = []
                with open('./' + fol + '/rmsds_and_final_scores.tsv', 'r') as infile:
                    for i, line in enumerate(infile):
                        if i != 0:
                            rmsds.append(float(line.split('\t')[3]))
                bests[ID] = min(rmsds)
            except:
                pass

#%%
    with open('../../benchmark/' +'APE-Gen_dataset.csv', 'r') as agfile:
        ag = {}
        for i, line in enumerate(agfile):
            if i != 0:
                row = line.split(',')
                ag[row[0].upper()] = float(row[2])
    
    ag_comparison = {}
   # print('Missing cases with APE-Gen:')
    for key in ag:
        try:
            ag_comparison[key] = ((bests[key], ag[key]))
        except:
           # print(key)
           continue
            
    return ag_comparison


#%%
#AG_comparison = [x for x in bests if x[0] < 12]

#with open('../../benchmark/' + outdir + 'APE-Gen_comparison.txt', 'w') as outfile:
#    print(outfile)
#    for i in [9, 10, 11]:
#        l = [x[1] for x in bests if x[0] == i]
#        outfile.write('%s-mers AVG best CA RMSD: ' %str(i) + str(sum(l)/len(l)) + '\n')
#    outfile.write('Overall AVG best CA RMSD: ' + str(sum([x[1] for x in AG_comparison])/len([x[1] for x in AG_comparison])))
#    outfile.close()

#plt.plot([ag_comparison[x][0] for x in ag_comparison], [ag_comparison[y][1] for y in ag_comparison], 'o')
#if label:
 #   for key in ag_comparison:
 #       plt.annotate(key, (ag_comparison[key][0], ag_comparison[key][1]))
        #ag
#plt.plot([0, 2, 4], [0, 2, 4], 'g-')
#plt.xlabel('PANDORA')
#plt.ylabel(' APE-Gen ' )
#plt.savefig('../../benchmark/' + outdir + 'CA_PANDORA_vs_APE-Gen.pdf')
#plt.show()
#plt.clf()

#plt.hist([ag_comparison[x][0] - ag_comparison[x][1] for x in ag_comparison], bins=50)
#plt.xlabel('PANDORA - APE-Gen')
#plt.ylabel(' # of cases ' )
#plt.axvline(0, color= 'r')
#plt.savefig('../../benchmark/' + outdir + 'CA_PANDORA_vs_APE-Gen_hist.pdf')
#plt.show()
#plt.clf()


#%%
### DockTope
def docktope_rmsd(best_rmsds):
    os.chdir('/home/rafaellab/pMHC_Modelling')
                  
    with open('./benchmark/' +'DockTope_benchmark_dataset_calc.tsv', 'r') as dtfile:
        dt = {}
        for line in dtfile:
            row = line.split(' ')
            dt[row[1]] = row[4]
    
    dt_comparison = {}
    #print('Missing cases with DockTope:')
    for key in dt:
        try:
            dt_comparison[key] = ((best_rmsds[key][0], float(dt[key])), (best_rmsds[key][3], float(dt[key])))
        except:
            #print(key)
            continue
    return dt_comparison

#%%
#plt.plot([dt_comparison[x][0][0] for x in dt_comparison], [dt_comparison[y][0][1] for y in dt_comparison], 'o')
#if label:
#    for key in dt_comparison:
#        plt.annotate(key, (dt_comparison[key][0][0], dt_comparison[key][0][1]))
#plt.plot([0, 2, 4], [0, 2, 4], 'g-')
#plt.xlabel('PANDORA molpdf')
#plt.ylabel(' DockTope ' )
#plt.savefig('../../benchmark/' + outdir + 'CA_PANDORA_molpdf_vs_DockTope.pdf')
#plt.show()
#plt.clf()

#plt.hist([dt_comparison[x][0][0] - dt_comparison[x][0][1] for x in dt_comparison], bins=50)
#plt.xlabel('PANDORA molpdf - DockTope')
#plt.ylabel(' # of cases ' )
#plt.axvline(0, color= 'r')
#plt.savefig('../../benchmark/' + outdir + 'CA_PANDORA_molpdf_vs_DockTope_hist.pdf')
#plt.show()
#plt.clf()

#plt.plot([dt_comparison[x][1][0] for x in dt_comparison], [dt_comparison[y][1][1] for y in dt_comparison], 'o')
#if label:
#    for key in dt_comparison:
#        plt.annotate(key, (dt_comparison[key][1][0], dt_comparison[key][1][1]))
#plt.plot([0, 2, 4], [0, 2, 4], 'g-')
#plt.xlabel('PANDORA DOPE')
#plt.ylabel(' DockTope ' )
#plt.savefig('../../benchmark/' + outdir + 'CA_PANDORA_DOPE_vs_DockTope.pdf')
#plt.show()
#plt.clf()

#plt.hist([dt_comparison[x][1][0] - dt_comparison[x][1][1] for x in dt_comparison], bins=50)
#plt.xlabel('PANDORA molpdf - DockTope')
#plt.ylabel(' # of cases ' )
#plt.axvline(0, color= 'r')
#plt.savefig('../../benchmark/' + outdir + 'CA_PANDORA_DOPE_vs_DockTope_hist.pdf')
#plt.show()
#plt.clf()

#%%
### GradDock
def graddock_rmsd(best_rmsds):
    os.chdir('/home/rafaellab/pMHC_Modelling')
    
    GD = []
    with open('./benchmark/' + 'GradDock_results_table.txt', 'r') as gdinfile:
        for i, line in enumerate(gdinfile):
            if i > 107:
                GD.append(float(line.split('\t')[0]))
    
    
    '''
    GD_names = {}
    with open('../../benchmark/' + 'GradDock_benchmark_dataset.tsv', 'r') as gdfile:
        for i, line in enumerate(gdfile):
            row = line.split(' ')
            GD_names[(row[1].split('-')[0]).upper()] = GD[i]
    '''
    
    
    surrogate = [
    ('4f7t_D','3wlb_A'), ('2mha_C','1vad_A'), ('3rwh_A','3rwd_A'), ('3bwa_A','2nw3_A'),
    ('2fo4_A','1leg_A'), ('3p4n_D','3pab_A'), ('3fon_C','3fom_A'), ('4eup_A','1i1y_A'),
    ('1a9b_A','3mv7_A'), ('1ydp_A','3kyo_A'), ('2c7u_D','4jfp_A'), ('3bp7_A','1of2_A'),
    ('1ld9_A','1ldp_H'), ('4nqx_K','1w72_D'), ('1p7q_A','3d25_A'), ('3rwg_A','3rwf_A'),
    ('3tbw_G','1jpf_A'), ('1xr9_A','3c9n_A'), ('1q94_D','1qvo_A'), ('2vlr_F','1hhk_A'),
    ('2vab_A','3roo_A'), ('2cik_A','3lkr_A'), ('4l29_O','1qew_A'), ('3bxn_A','3bvn_A'),
    ('2bck_A','3i6l_D'), ('1fzo_A','1rjz_A'), ('1eez_D','1i1f_A'), ('3qeq_A','4l29_O'),
    ('1w72_D','4nqx_K'), ('1sys_A','4jqx_A'), ('3upr_C','3vrj_A'), ('2v2x_D','1hhi_A'),
    ('3lkr_A','2cik_A'), ('2xpg_A','3rl2_A'), ('3rwf_A','3rwg_A'), ('3fqt_A','2pye_A'),
    ('1jpg_A','1s7x_A'), ('1wbz_C','3p4m_A'), ('2x4s_A','1im3_A'), ('3ffc_F','1m05_A'),
    ('4l8d_C','3tbw_A'), ('1ldp_H','1ld9_A'), ('1ogt_A','4g8g_A'), ('3mrc_A','1qr1_A'),
    ('3czf_A','3bp7_A'), ('1i1y_A','4eup_A'), ('3gsw_A','3giv_A'), ('3to2_A','4gkn_A'),
    ('1xr8_A','3c9n_A'), ('1ce6_A','3quk_A'), ('3hg1_A','1eey_A'), ('2yf5_A','3bew_A'),
    ('2axf_A','3vft_A'), ('4g9f_A','1jge_A'), ('3giv_D','2gj6_A'), ('3vxr_A','4f7m_D'),
    ('1qvo_A','1q94_D'), ('4f7m_D','3vxr_A'), ('1i4f_A','1jht_A'), ('3vrj_A','3upr_C'),
    ('1n3n_A','4huv_A'), ('3utt_F','3mrp_A'), ('3mgt_A','2git_A'), ('1jgd_A','3bp4_A'),
    ('3wlb_A','4f7t_D'), ('2axg_A','3mv7_A'), ('2hn7_A','1q94_D'), ('3bew_D','2yf5_A'),
    ('2clr_D','3mgt_A'),
    ]
    
    GD_names = {}
    for i, case in enumerate(surrogate):
        GD_names[case[0][:4].upper()] = (GD[i], case)
    
    print()
    #print('Missing cases with GradDock:')
    gd_comparison = {}
    for key in GD_names:
        try:
            gd_comparison[key] = ((best_rmsds[key][2], float(GD_names[key][0])), (best_rmsds[key][4], float(GD_names[key][0])))
        except:
            #print(key)
            continue
            
            
    return gd_comparison


#%%
#plt.plot([gd_comparison[x][0][0] for x in gd_comparison], [gd_comparison[y][0][1] for y in gd_comparison], 'o')
#if label:
#    for key in gd_comparison:
#        plt.annotate(key, (gd_comparison[key][0][0], gd_comparison[key][0][1]))
#plt.plot([0, 2, 4], [0, 2, 4], 'g-')
#plt.xlabel('PANDORA molpdf')
#plt.ylabel(' GradDock ' )
#plt.savefig('../../benchmark/' + outdir + 'CA-BB-CB_PANDORA_molpdf_vs_GradDock.pdf')
#plt.show()
#plt.clf()

#plt.hist([gd_comparison[x][0][0] - gd_comparison[x][0][1] for x in gd_comparison], bins=50)
#plt.xlabel('PANDORA molpdf - GradDock')
#plt.ylabel(' # of cases ' )
#plt.axvline(0, color= 'r')
#plt.savefig('../../benchmark/' + outdir + 'CA-BB-CB_PANDORA_molpdf_vs_GradDock_hist.pdf')
#plt.show()
#plt.clf()

#plt.plot([gd_comparison[x][1][0] for x in gd_comparison], [gd_comparison[y][1][1] for y in gd_comparison], 'o')
#if label:
#    for key in gd_comparison:
#        plt.annotate(key, (gd_comparison[key][1][0], gd_comparison[key][1][1]))
#plt.plot([0, 2, 4], [0, 2, 4], 'g-')
#plt.xlabel('PANDORA DOPE')
#plt.ylabel(' GradDock ' )
#plt.savefig('../../benchmark/' + outdir + 'CA-BB-CB_PANDORA_DOPE_vs_GradDock.pdf')
#plt.show()
#plt.clf()

#plt.hist([gd_comparison[x][1][0] - gd_comparison[x][1][1] for x in gd_comparison], bins=50)
#plt.xlabel('PANDORA DOPE - GradDock')
#plt.ylabel(' # of cases ' )
#plt.axvline(0, color= 'r')
#plt.savefig('../../benchmark/' + outdir + 'CA-BB-CB_PANDORA_DOPE_vs_GradDock_hist.pdf')
#plt.show()
#plt.clf()


#%%

#plt.figure(num=None, figsize=(20, 20), dpi=80, facecolor='w', edgecolor='k')
#plt.plot([x[1] for x in best_rmsds.values()], [x[4] for x in best_rmsds.values()], 'o')
#if label:
#    for key in best_rmsds:
#        plt.annotate(key, (best_rmsds[key][1], best_rmsds[key][4]), size = 8)
#plt.xticks([0.25*x for x in range(20)])
#plt.yticks([0.25*x for x in range(20)])
#plt.axvline(2, color='r')
#plt.axhline(2, color='r')
#plt.xlabel('molpdf top 5')
#plt.ylabel('DOPE top 5')
#plt.plot([0, 2, 4, 5], [0, 2, 4, 5], 'g-')
#plt.savefig('../../benchmark/' + outdir + 'PANDORA_CA_molpdf_vs_DOPE.pdf')
#plt.show()
#plt.clf()

#plt.figure(num=None, figsize=(20, 20), dpi=80, facecolor='w', edgecolor='k')
#plt.plot([x for x in bests.values()], [x for x in bests.values()], 'o')
#if label:
#    for key in bests:
#        plt.annotate(key, (bests[key], bests[key]), size = 8)
#plt.xticks([0.25*x for x in range(20)])
#plt.yticks([0.25*x for x in range(20)])
#plt.axvline(2, color='r')
#plt.axhline(2, color='r')
#plt.xlabel('best RMSD model')
#plt.ylabel('best RMSD model')
#plt.plot([0, 2, 4, 5], [0, 2, 4, 5], 'g-')
#plt.savefig('../../benchmark/' + outdir + 'PANDORA_CA_best_rmsd.pdf')
#plt.show()
#plt.clf()

#plt.hist([x for x in bests.values()], bins = 50)
#plt.axvline(2, color='r')
#plt.xlabel('Best model i-RMSD')
#plt.ylabel('# of cases')
#plt.savefig('../../benchmark/' + outdir + 'PANDORA_CA_best_rmsd_hist.pdf')
#plt.show()
#plt.clf()
#%%
def Z_test (before_13, after_13):
    import statsmodels.stats.api as sms
    
    cm = sms.CompareMeans(sms.DescrStatsW(before_13),sms.DescrStatsW( after_13))
    z,pval = cm.ztest_ind( alternative = 'two-sided', usevar = 'unequal')
    
    return z, pval
    #H0: means of distributions are not statistically different
    #at 0.05: Reject H0 if Z < -1.960 or if Z > 1.960.
    #p>0.05: reject H0
   # print('z: {} , pval: {}'.format(z, pval))

#%%calling the function, calculating best rmsds

#directory1 = sys.argv[1]
#directory2 = sys.argv[2]

directory2 = '/home/dariom/PANDORA_master/outputs/bm_issue_21_20201229' 
directory1 = '/home/dariom/PANDORA_issue_21_before_issue_13/outputs/bm_issue_13_20201223'

best_rmsds_1= pandora_rmsd(directory1)
best_rmsds_2 = pandora_rmsd(directory2)


before_13= []
for i in best_rmsds_1:
    before_13.append (min(best_rmsds_1[i]))  

after_13 = []
for i in best_rmsds_2:
    after_13.append (min(best_rmsds_2[i]))

z,pval = Z_test (before_13, after_13)   
print('Best RMSDS: z:', z, 'p:', pval)


#%% ape-gen 
ag_comparison_1 = apgen_rmsd(best_rmsds_1, directory1)
before_13_ag = []
for x in ag_comparison_1:
   before_13_ag.append(ag_comparison_1[x][0] - ag_comparison_1[x][1])
   
ag_comparison_2 = apgen_rmsd(best_rmsds_2, directory2)
after_13_ag = []
for x in ag_comparison_2:
   after_13_ag.append(ag_comparison_2[x][0] - ag_comparison_2[x][1])
   
z,pval = Z_test (before_13_ag, after_13_ag)   
print('APEGEN comparison: z:', z, 'p:', pval)

#%%dock tope

dt_comparison_1 = docktope_rmsd(best_rmsds_1)
dt_comparison_2 = docktope_rmsd(best_rmsds_2)

before_13_dt = []
for x in dt_comparison_1:
    before_13_dt.append (dt_comparison_1[x][0][0] - dt_comparison_1[x][0][1])
after_13_dt = []
for x in dt_comparison_2:
    after_13_dt.append (dt_comparison_2[x][0][0] - dt_comparison_2[x][0][1] )

z,pval = Z_test (before_13_dt, after_13_dt)   
print('Pandora Molpdf - DockTope comparison: z:', z, 'p:', pval)

before_13_dt = []
for x in dt_comparison_1:
    before_13_dt.append (dt_comparison_1[x][1][0] - dt_comparison_1[x][1][1])
after_13_dt = []
for x in dt_comparison_2:
    after_13_dt.append (dt_comparison_2[x][1][0] - dt_comparison_2[x][1][1] )

z,pval = Z_test (before_13_dt, after_13_dt)   
print('Pandora Dope- DockTope comparison: z:', z, 'p:', pval)

#%%
gd_comparison_1 = graddock_rmsd(best_rmsds_1)
gd_comparison_2 = graddock_rmsd(best_rmsds_2)

before_13_gd = []
for x in gd_comparison_1:
    before_13_gd.append (gd_comparison_1[x][0][0] - gd_comparison_1[x][0][1])
after_13_gd = []
for x in gd_comparison_2:
    after_13_gd.append (gd_comparison_2[x][0][0] - gd_comparison_2[x][0][1] )

z,pval = Z_test (before_13_gd, after_13_gd)   
print('Pandora Molpdf - GradDock comparison: z:', z, 'p:', pval)

before_13_gd = []
for x in gd_comparison_1:
    before_13_gd.append (gd_comparison_1[x][1][0] - gd_comparison_1[x][1][1])
after_13_gd = []
for x in gd_comparison_2:
    after_13_gd.append (gd_comparison_2[x][1][0] - gd_comparison_2[x][1][1] )

z,pval = Z_test (before_13_gd, after_13_gd)   
print('Pandora Dope - GradDock comparison: z:', z, 'p:', pval)
