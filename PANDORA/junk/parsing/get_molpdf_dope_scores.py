import sys
import csv


modeller_output_file = sys.argv[1]


with open(modeller_output_file) as result_file:
    out_file=result_file.read()
    if '>> Summary of successfully produced loop models:' in out_file:
        parts= out_file.split('>> Summary of successfully produced loop models:')

        if '>> Summary of failed loop models:' in parts[-1]:
            new_parts= parts[-1].split('>> Summary of failed loop models:')
            models=new_parts[0].splitlines()[3:-2]
        else:
            models=parts[-1].splitlines()[3:-1]

        with open('molpdf_DOPE.tsv',"wt") as f:
            tsv_writer = csv.writer(f, delimiter='\t')
            tsv_writer.writerow(['MODEL', 'molpdf', 'DOPE'])
            to_be_written = []
            for line in models:
                row = [x for x in line.split(' ') if x != '']
                to_be_written.append((float(row[1]),row))
            to_be_written = sorted(to_be_written)
            for row in to_be_written:
                tsv_writer.writerow(row[1])
