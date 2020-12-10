# -*- coding: utf-8 -*-
import urllib.request
import urllib.parse
import os
from copy import deepcopy
from pyparsing import nestedExpr
import time

print_outfiles = True

#def download_ids_alleles_imgt(out_tsv = 'auto_generated_IDs_alleles_from_IMGT.tsv', out_pkl = 'IDs_and_alleles_identity_percs_from_imgt.pkl', print_outfiles = False):
cwd = os.getcwd()

'''
params = { 'ReceptorType' : 'peptide/MH1',
        'type-entry': 'PDB'}
'''
params = { 'ReceptorType' : 'MH2',
    'type-entry': 'PDB'}

url = "http://www.imgt.org/3Dstructure-DB/cgi/3Dquery.cgi"

data = urllib.parse.urlencode(params)
data = data.encode('ascii') # data should be bytes
req = urllib.request.Request(url, data)

IDs_list = []


with urllib.request.urlopen(req) as response:
    text = response.read().decode('utf-8')
    text = text.splitlines()
    if print_outfiles: 
        temp_outfile = open('../outputs/test/test_download/test.html', 'w')
        temp_outfile.write(('\n').join(text))
        temp_outfile.close()

IDs_list = [x for x in text if 'href' in x and 'pdbcode' in x]
IDs_list = [x.split('"') for x in IDs_list]
IDs_list = [x[3][-4:] for x in IDs_list]

time.sleep(10)

ID_allele = {}
for ID in IDs_list[:10]:
    
    time.sleep(10)
    str_url = 'http://www.imgt.org/3Dstructure-DB/cgi/details.cgi?pdbcode=%s' %ID
    #str_url = "http://www.imgt.org/3Dstructure-DB/cgi/details.cgi?pdbcode=6ujo"
    
    req = urllib.request.Request(str_url)
    
    with urllib.request.urlopen(req) as response: #TODO: Handle exception "HTTPError: Internal Server Error"
        text = response.read().decode('utf-8')
        text = text.splitlines()
        if print_outfiles: 
            temp_outfile = open('../outputs/test/test_download/%s.html' %ID, 'w')
            temp_outfile.write(('\n').join(text))
            temp_outfile.close()
    
    domain_flag = False
    domains = []
    domain = []
    for line in text:
        if 'collier_perles' in line:
            continue
        if line == '<td class="titre_h" align="center" rowspan="6">G-DOMAIN</td>':
            domain_flag = True
            #print(line)
        elif line =='\n' or line =='':
            domain_flag = False
            if domain != []:
                domains.append(domain)
                domain = []
        if domain_flag:
            domain.append(line)
    
    for i, dom in enumerate(domains):
        flag = False
        descr = 0
        for j, data in enumerate(dom):
            if flag == True and j == (descr + 2):
                #print('Is this what are you looking for? :', data)
                flag = False
                domains[i] = data.split('title')[0]
                continue
            if 'IMGT gene and allele name' in data:
                descr = j
                flag = True
    
    percs = [{}, {}]
    clean_domains = deepcopy(domains)
    for i, domain in enumerate(domains):
        clean_domains[i] = domains[i].replace(',','').replace(';','').split('&nbsp')
    for i, domain in enumerate(clean_domains[:2]):
        for j, d in enumerate(clean_domains[i]):
            if '%' in d:
                percs[i][clean_domains[i][j-1]] = float((nestedExpr('(',')').parseString(d).asList())[0][0].replace('%',''))
        clean_domains[i] = [x for x in clean_domains[i] if '%' not in x]
    
    
    # TODO: in future: in case of ambiguity, assign allele to the templates according to patient genotype?
    
    if percs != [{}, {}]:
        set_percs = set(percs[0]).intersection(*percs)
        
            
    ID_allele[ID] = [set_percs, percs] #TODO: check if allele is in both domains

'''    
if out_pkl:
    with open('%s/data/csv_pkl_files/%s' %(cwd, out_pkl), 'wb') as outpkl:
        pickle.dump(ID_allele, outpkl)

if out_tsv:
    with open('%s/data/csv_pkl_files/%s' %(cwd, out_tsv), 'w') as outtsv:
        outtsv.write('PDB ID' + '\t' + 'Equal identity alleles' + '\n')
        for ID in ID_allele:
            if len(ID_allele[ID][0]) != 0:
                outtsv.write(ID + '\t' + (';').join([x for x in ID_allele[ID][0]]) +'\n')
            elif len(ID_allele[ID][0]) == 0:
                to_write = []
                for x in ID_allele[ID][1]:
                    to_write += list(x.keys())
                to_write = list(set(to_write))
                outtsv.write(ID + '\t' + (';').join(to_write) +'\n')
'''