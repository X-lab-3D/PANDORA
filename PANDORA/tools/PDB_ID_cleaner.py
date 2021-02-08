import pickle
import PANDORA

IDd = open(PANDORA.PANDORA_data+"/csv_pkl_files/IDs_ChainsCounts_dict.pkl", "rb")
idd = pickle.load(IDd)
IDd.close()

match = []
for key in idd:
    chaind = idd[key]
    value = chaind[min(chaind, key=chaind.get)]
    if value > 20:
        print('NoPeptide:', key, idd[key])
        match.append(key)
print(len(match))
print(match)