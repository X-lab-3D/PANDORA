import csv

flag = False
with open('GradDock_benchmark_dataset.tsv', 'wt') as out:
    tsv_writer = csv.writer(out, delimiter='\t')
    tsv_writer.writerow(['Organism', 'Allele', 'Structure ID'])
    with open('GradDock_benchmarck_dataset.txt', 'r') as f:
        for line in f:
            l = line.replace(',','')
            k = l.replace('\n', '')
            j = k.split(' ')
            if j[0] in ['Human', 'Mouse', 'Rat', 'Cow', 'Pig', 'Chicken', 'Monkey']:
                prev_line = j
                flag = True
            else:
                j = prev_line + j
            if flag == True:
                for ID in j[2:]:
                    tsv_writer.writerow([j[0], j[1], ID])
                    flag = False