import argparse

arg_parser = argparse.ArgumentParser(description="Sets modeller license key")
arg_parser.add_argument(
    "--license","-l",
    help="MODELLER license key code",
    required=True
)
a = arg_parser.parse_args()

modeller_config_file = '/usr/share/miniconda/envs/test/lib/modeller-10.2/modlib/modeller/config.py'

modeller_config = []
with open(modeller_config_file, 'r') as f:
    for line in f:
        modeller_config.append(line)
with open(modeller_config_file, 'w') as f:
    print("Writing MODELLER config file")
    for line in f:
        f.write(line.replace("XXXX", a.license))
        print(line.replace("XXXX", a.license))