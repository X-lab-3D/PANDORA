import csv

flag = False
with open('GradDock_complete_dataset.tsv', 'wt') as out:
    tsv_writer = csv.writer(out, delimiter='\t')
    tsv_writer.writerow(['Organism', 'Allele', 'Structure ID'])
    with open('GradDock_complete_dataset.txt', 'r') as f:
        prev_line = []
        lines = []
        for i, line in enumerate(f):
            l = line.replace(',','')
            k = l.replace('\n', '')
            j = k.split(' ')
            if j[0] in ['Human', 'Mouse', 'Rat', 'Cow', 'Pig', 'Chicken', 'Monkey']:
                '''
                if flag == True:
                    raise Exception('ok')
                    for ID in prev_line[2:]:
                        tsv_writer.writerow([prev_line[0], prev_line[1], ID])
                    flag = False
                '''
                if i != 0:
                    lines.append(prev_line)
                prev_line = j
                #flag = True
            else:
                prev_line = prev_line + j
        for line in lines:
            for ID in line[2:]:
                tsv_writer.writerow([line[0], line[1], ID])
                #flag = False