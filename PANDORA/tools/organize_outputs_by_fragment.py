import os


def organize_one_protein_output(indir):
    cwd = os.getcwd()
    os.chdir(cwd + '/' + indir)
    frags = []
    frag_to_names = {}
    for name in os.listdir('./'):
        if name.startswith('mut_') or name.startswith('sel_'):
            for f in os.listdir('./' + name):
                if f.endswith('.ali'):
                    with open('./' + name + '/' + f) as infile:
                        L = ''
                        for line in infile:
                            L += line
                    frag = L.split('/')[-1].replace('*', '')
                    frags.append(frag)
                    try:
                        frag_to_names[frag].append(name)
                    except:
                        frag_to_names[frag] = [name]
    frags = list(set(frags))

    for frag in frags:
        try:
            os.mkdir(frag)
        except FileExistsError:
            pass
        for name in frag_to_names[frag]:
            os.system('mv %s ./%s'%(name, frag))
    os.chdir(cwd)


for fol in os.listdir('./'):
    if len(fol) == 6 and '.' not in fol:
        organize_one_protein_output(fol)
