import sys

modeller_config_file = '/usr/share/miniconda/envs/test/lib/modeller-10.2/modlib/modeller/config.py'
modeller_license_key = sys.argv[1]

modeller_config = []
with open(modeller_config_file, 'r') as f:
    for line in f:
        modeller_config.append(line)
with open(modeller_config_file, 'w') as f:
    print("Writing MODELLER config file")
    for line in f:
        f.write(line.replace("XXXX", modeller_license_key))
        print(line.replace("XXXX", modeller_license_key))
