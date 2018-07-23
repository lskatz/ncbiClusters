#!/usr/local/bin/python
#|__Khushbu Patel
#|__Load module Python/3.4
# Usage: ./linelist.py


from Bio import Entrez
import glob
import pandas as pd



temp = []
line2 = []
biosample_acc = []
wgs_id = []
bio_SAMN = {}
bio1 = []
bio2 = []
gdist = []		# stores all values in compatible_distance column 

# ------------------------------------------------------------------------------------------------------------------------- #
# Reading the line list file

data = pd.read_excel('1804MLJMP-1_Montevideo_analreq180709.xlsx')
temp = data['WGS_id'].values.tolist() # NOTE: make sure WGS Id column header in all line lists docs is named the same way

for i in temp:
	i = str(i)
	if i.startswith('PNU'): 
		wgs_id.append(i)
	else:
		continue

for id in wgs_id:
	Entrez.email = 'oix2@cdc.gov'
	handle = Entrez.efetch('BioSample', id = id, retmode='text')
	line1 = handle.readline()
	line2 = (handle.readline().split())	# contains the biosample ID
	line3 = handle.readline()

	id = line2[2] # third element is the biosample accession
	biosample_acc.append(id[:-1]) # to slice off last character from "SAMN00039977;"


print("Printing WGS_id with their corresponding SAMN accession - \n")
for id,samn in zip(wgs_id,biosample_acc):
	print(id,samn)
	

# Reading the metadata file
data = pd.read_csv('PDG000000002.1169.reference_target.SNP_distances.tsv',sep = '\t',error_bad_lines = False)

data['biosample_acc_1'] = data['biosample_acc_1'].astype('|S')
data['biosample_acc_2'] = data['biosample_acc_2'].astype('|S')

bio1 = data['biosample_acc_1'].values.tolist()
bio2 = data['biosample_acc_2'].values.tolist()
gdist = data['compatible_distance'].values.tolist()	# genomic distances added to the list


print("Printing matches between linelist and metadata file - \n")
for a1,a2 in zip(bio1,bio2):
	for x in biosample_acc:
		if(x == a1):
			print(x,a2)		# Linking WGS SAMN's with bio acc 2 in tsv
		elif( x == a2):
			print(x,a1)		# Linking WGS SAMN's with bio acc 1 in tsv

