# -*- coding: utf-8 -*-

import csv

infile = open('../outputs/benchmark/best_RMSDs.tsv', 'r')
r = csv.reader(infile, delimiter='\t')
qfile = open('../benchmark/PANDORA_benchmark_dataset.tsv', 'r')
q = csv.reader(qfile, delimiter='\t')
outfile = open('../outputs/benchmark/id_best_RMSDs.tsv', 'wt')
tw = csv.writer(outfile, delimiter='\t')

d = {}
for k, line in enumerate(q):
    if k == 0:
        pass
    else:
        query = 'query_' + str(k)
        d[query]=line[0]

for i, row in enumerate(r):
    if i == 0:
        tw.writerow(row)
    else:
        tw.writerow([d[row[0]], row[1], row[2]])
qfile.close()
outfile.close()
infile.close()