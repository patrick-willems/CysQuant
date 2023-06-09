#!/usr/bin/env python

"""Script to calculate the cysteine oxidation degree from the MaxQuant peptides.txt output.
This script is only sensible if light and heavy iodoacetamide labels were used within samples. 
Please see xxx for more information."""

__author__      = "Patrick Willems"
__email__       = "willems533@gmail.com"
__copyright__   = "Copyright 2023, Ghent University"

import re, sys, os, argparse
parser = argparse.ArgumentParser()
parser.add_argument("-in", "--input_file", required=True, help="Location of output 'peptides.txt' of labeled MaxQuant run (Cys IAM labels have to be set!).")
parser.add_argument("-fasta", "--fasta_file", required=True, help="UniProtKB FASTA file used in the search.")
parser.add_argument("-ox", "--label_ox", required=True, choices=['L', 'H'], help="Specify the label (either 'L' or 'H' for light/heavy) that was used to label oxidized cysteines (i.e. after DTT/TCEP reduction) .")
args = parser.parse_args()

#Store UniProt protein sequences
print('#Reading FASTA')
prot = {}
try: fasta = open(args.fasta_file, 'r')
except OSError: sys.exit("Could not open/read fasta file provided to script (-fasta): "+args.fasta_file)
for line in fasta:
    line = line.rstrip()
    if '>' in line:
        protein = re.match(r'>\w+\|(\S+)\|\S+_ARATH',line).group(1)
        prot.update({protein: ''})
    else: prot[protein] += line
print('  Stored '+str(len(prot))+' proteins from '+args.fasta_file+'.')

print('#Reading MaxQuant peptides.txt output')
try: f = open(args.input_file, 'r')
except OSError: sys.exit("Could not open peptides.txt file location provided to script (-in): "+args.input_file)
#Get i of required columns from the header line - quite long but should be more robust
ind = {}       
header = next(f)
headers = header.split('\t')
samples = {}    #Keep track of all samples in the report
summaries = {}  #Oxidation statistics per sample
for i in range(len(headers)):
    if 'C Count' in headers[i]: ind.update({'C': i})
    elif 'Sequence' in headers[i]: ind.update({'seq': i})
    elif 'Gene names' in headers[i]: ind.update({'gene': i})
    elif 'Protein names' in headers[i]: ind.update({'name': i})
    elif 'Proteins' == headers[i]: ind.update({'prot': i})
    elif 'Ratio H/L' in headers[i] and len(headers[i].split()) == 3 and headers[i] not in ['Ratio H/L normalized','Ratio H/L count','Ratio H/L iso-count','Ratio H/L type']:
        sample = (headers[i].split())[2]
        ind.update({sample: i})
        samples.update({sample: True})
        summaries.update({sample: {'nr': 0, 'ox': 0}})

data = {}       #Here all actual cysteine data will be stored
for line in f:
    line = line.rstrip()
    
    #All required data is stored to variables based on the header indices parsed above
    col = line.split('\t')
    seq,nr_cys,prot,gene,name = col[ind['seq']],col[ind['C']],col[ind['prot']],col[ind['gene']],col[ind['name']]
    if nr_cys == '1' and '__REV' not in line and '__CON' not in line and prot != '':    #No contaminant,reverse hits
        descr = gene+'; '+name
        for pattern in ['_ARATHDbj','_ARATHIsamemberofthePF','_ARATHGb']:  #Sometimes incorrect protein name encountered
            prot = prot.replace(pattern,'')
        data.update({seq: {'prot': prot, 'descr': descr, 'samples':{}}})
        for sample in samples:
            if args.label_ox == 'H': ox = 1 - (1 / ( 1 + float(col[ind[sample]])))
            else: ox = (1 / ( 1 + float(col[ind[sample]])))
            if ox == ox:  #NA test
                data[seq]['samples'].update({sample: ox})
                summaries[sample]['ox'] += ox
                summaries[sample]['nr'] += 1
    else: continue

#Print out summaries
print('#Writing out cysteine oxidation quantification summaries per sample:')
for run in summaries:
    average = 100*round((summaries[run]['ox']/summaries[run]['nr']),5)
    print('  '+run+': '+str(summaries[run]['nr'])+' ratios, '+str(average)[0:6]+'% Cys OX')

#Print out reports
outfile = 'ms1quant_' + args.input_file
print('#Writing out MS1 cysteine oxidation report to '+outfile)
out = open(outfile,'w')
samples_s = {k: samples[k] for k in sorted(samples)}    #Sort the sample alphabetically
out.write("Peptide\tProteins\tDescription_first\tProteins_position\t" + '\t'.join(samples_s) + '\n')

#Print quantifications
for seq in data:
    proteins = data[seq]['prot'].split(';')

    #Determine cysteine position(s) in protein
    position = ''
    for protein in proteins:
        if protein not in prot: continue
        posCys = str(seq.find('C') + prot[protein].find(seq) + 1)
        position += posCys + ';'

    #First general peptide precursor fields
    to_write = seq+'\t'+data[seq]['prot']+'\t'+data[seq]['descr']+'\t'+position[:-1]
    
    #Per sample MS1 quantification
    valid_values = 0
    for sample in samples_s: 
        if sample in data[seq]['samples']:
            to_write += '\t' + str(data[seq]['samples'][sample])
            valid_values += 1
        else: to_write += '\tNA'
    
    if valid_values > 0: out.write(to_write+'\n')