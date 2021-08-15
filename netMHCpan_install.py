
"""
The user needs to manually download the netMHCpan or netMHCIIpan software, since it requires agreement to an
academic license agreement

1. Go to: https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1 and https://services.healthtech.dtu.dk/service.php?NetMHCIIpan-4.0
2. press the download button
3. Download the most recent version appropriate for your operating system
4. Untar the file and put it in the root directory of your PANDORA install
5. Follow the readme or simply run this file to configure netMHCpan/netMHCIIpan
"""


import os
wd = os.path.dirname(os.path.abspath(__file__))

# netMHC Pan install

try:
# find version
    netMHCpan = [i for i in os.listdir(wd) if i.startswith("netMHCpan") and os.path.isdir(i)][0]
    # Changing working directory
    os.chdir(wd + '/' + netMHCpan)
    # Downloading data folder for netMHCIIpan
    os.system('wget https://services.healthtech.dtu.dk/services/' + netMHCpan[0].upper() + netMHCpan[1:] + '/data.tar.gz ')
    # Uncompressing
    os.system('gunzip data.tar.gz')
    os.system('tar -xvf data.tar')
    # remove tar file
    os.system('rm data.tar')
    # Make sure the tmp directory is created
    if not os.path.exists('tmp'):
        os.mkdir('tmp')

    # Edit paths in the netMHCIIpan script
    netMHCpan_script = []
    with open('netMHCpan') as f:
        for line in f:
            netMHCpan_script.append(line)
    with open('netMHCpan', 'w') as f:
        for line in netMHCpan_script:
            if 'setenv\tNMHOME' in line:
                f.write('setenv\tNMHOME\t' + wd + '/' + netMHCpan + '\n')
                f.write('setenv\tTMPDIR\t' + wd + '/' + netMHCpan + '/tmp\n')
            else:
                f.write(line)

    # Changing working directory back
    os.chdir(wd)
except:
    print('Something went wrong with configuring netMHCpan')
    os.chdir(wd)

# netMHCII Pan install

try:
    # find version
    netMHCIIpan = [i for i in os.listdir(wd) if i.startswith("netMHCIIpan") and os.path.isdir(i)][0]
    # Changing working directory
    os.chdir(wd + '/' + netMHCIIpan)
    # Downloading data folder for netMHCIIpan
    os.system('wget https://services.healthtech.dtu.dk/services/' + netMHCIIpan[0].upper() + netMHCIIpan[1:] + '/data.tar.gz')
    # Uncompressing
    os.system('gunzip data.tar.gz')
    os.system('tar -xvf data.tar')
    # remove tar file
    os.system('rm data.tar')
    # Make sure the tmp directory is created
    if not os.path.exists('tmp'):
        os.mkdir('tmp')

    # Edit paths in the netMHCIIpan script
    netMHCIIpan_script = []
    with open('netMHCIIpan') as f:
        for line in f:
            netMHCIIpan_script.append(line)
    with open('netMHCIIpan', 'w') as f:
        for line in netMHCIIpan_script:
            if 'setenv\tNMHOME' in line:
                f.write('setenv\tNMHOME\t' + wd + '/' + netMHCIIpan + '\n')
                f.write('setenv\tTMPDIR\t' + wd + '/' + netMHCIIpan + '/tmp\n')
            else:
                f.write(line)

    # Changing working directory back
    os.chdir(wd)
except:
    print('Something went wrong with configuring netMHCIIpan')
    os.chdir(wd)

# os.popen("alias KEY_MODELLER='XXXX'").read()
# os.popen("conda install -y -c salilab modeller").read()
# os.popen("conda install -y -c conda-forge biopython").read()
# os.popen("conda install -y dill").read()
# os.popen("conda install -y -c bioconda muscle").read()
# os.popen("pip install pdb-tools").read()
# os.popen("pip install pdb2sql").read()
#
# os.popen("pip install -e ./").read()