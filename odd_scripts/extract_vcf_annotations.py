import vcf
import sys
import pandas as pd

vcf_files = sys.argv[1:]

records = {}

columns = ["Position", "Gene", "AA Change", "Is synonymous?"]

for file in vcf_files:
	name = file.split("/")[-1].replace(".annotated.vcf", "")
	reader = vcf.Reader(open(file, "r"))
	records[name] = {}
	for record in reader:
		if record.INFO["IsSynonymous"] != 9:
			records[name][str(record.INFO["CodonPosition"])] = [str(record.INFO["CodonPosition"]), str(record.INFO["Product"][0]).replace("[space]", " "), record.INFO["AminoAcidChange"][0], str(record.INFO["IsSynonymous"]).replace("1", "Y").replace("0", "N")]

dataframes = {}
names = []
df_list = []
			
for file in vcf_files:
	name = file.split("/")[-1].replace(".annotated.vcf", "")
	names.append(name)
	# dataframes[name] = pd.DataFrame(records[name], columns=columns)
	df_list.append(pd.DataFrame.from_dict(records[name], columns=columns, orient="index"))

df = df_list[0] #<- list of dataframes which correspond to a list of names (sample IDs)

for index, name in enumerate(names[1:]): #Join them together iteratively so I can use each name as a suffix
	if index == 0:
		df = df.join(df_list[index + 1], lsuffix=" " + names[0], rsuffix=" " + name, how="outer", on="Position") #For the first two add a suffix to both
	else:
		join = df_list[index + 1].add_suffix(" " + name)
		df = df.join(join, how="outer", on="Position")#For the 3rd onwards only add a suffix to the new one

for name in names:
	df.drop(labels="Position " + name, inplace=True, axis=1)

df["Position"] = pd.to_numeric(df["Position"], downcast="integer")
df.sort_values(inplace=True, by="Position")
df.to_csv("output.csv", index=False)