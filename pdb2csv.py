""" 
################################################################################
Created on Fri Nov 22 16:28:47 2013

@author: Spyros Charonis

This program takes as input a PDB file and extracts data such as K:R ratio,
D:E ratio, electrostatics-based metrics, sequence composition properties, etc.

It also extracts the sequence of the protein by converting 3-letter
AA Codes into their corresponding 1-letter code

If requested by the user, the script will write out the extracted data to 
a csv file

This script is invoked as: python pdb2csv.py *.pdb

################################################################################
""" 

#!/usr/bin/python

# Imports the pdb file as a command line argument called pdb
from sys import argv
script, pdbfi = argv

with open(pdbfi, 'rb') as fi:
    splitfi=fi.read()
    splitfi=splitfi.split("\n")[0:]

# Use TITLE field to extract protein source
title=[]
for i in range(len(splitfi)):
    if splitfi[i][0:5]=="TITLE":
        title.append(splitfi[i])
title = ''.join(title).strip('TITLE').strip()

## Use EXPDTA field to report structure determination method
expdta=[]
for i in range(len(splitfi)):
    if splitfi[i][0:6]=="EXPDTA":
        expdta.append(splitfi[i])
expdta = ''.join(expdta).strip('EXPDTA').strip()

## For NMR structures, check if multiple models (e.g. 2KH2)
## Use ATOM field to find alpha-carbons and extract all amino acid residues
aa_seq=[]
annot='nmr'.upper()

def extract_single_model():
    """ For NMR structures, this function extracts ATOM fields from
        only one MODEL; use ENDMDL field as a delimiter to terminate
        sequence list at the end of the first MODEL """
    for i in range(len(splitfi)):
        if splitfi[i][0:4]=="ATOM" and " CA " in splitfi[i]:
            aa_seq.append(splitfi[i])
        elif splitfi[i][0:6]=="ENDMDL":
            break
    return aa_seq

# if NMR structure, call function extract_single_model()
if annot in expdta:
    extract_single_model()
# if X-ray crystal structure, go straight for ATOM fields
else:
    for i in range(len(splitfi)):
        if splitfi[i][0:4]=="ATOM" and " CA " in splitfi[i]:
            aa_seq.append(splitfi[i])

## Check protein sequence for presence of unknown amino acids
def screen_unknown_res(seq):
    """ Screen PDB file for UNK amino acids and report how many if
        any are found """
    unk = []
    for i in range(len(aa_seq)):
        if 'UNK' in aa_seq[i]:
            unk.append(aa_seq[i])
    if len(unk):
        print '\n', 'WARNING: Protein sequence contains', len(unk), \
        'unkwnown amino acids', '\n', 'pdb2csv.py CANNOT BE USED!'
    elif not len(unk):
        pass
    return

screen_unknown_res(aa_seq)

## Compute KR-ratio
K=[]
for i in range(len(aa_seq)):
    if ' LYS ' in aa_seq[i]:
        K.append(aa_seq[i])

R=[]
for i in range(len(aa_seq)):
    if ' ARG ' in aa_seq[i]:
        R.append(aa_seq[i])

## Handle arithemtic pitfalls, i.e. 
## if no LYS or ARG residues in sequence, K:R = 0
## if no ARG residues, then define K:R to have value of 0.999

assert (len(K)==0 and len(R)==0) == False, "WARNING: No LYS or ARG \
                                              present in sequence"

if (len(K)==0 and len(R)==0):
    KR_ratio = float(0)
else:
    try:
        KR_ratio = (float(len(K)) / float(len(R)))
    except ZeroDivisionError:
        KR_ratio = float(.999)

## Compute DE-ratio
D=[]
for i in range(len(aa_seq)):
    if ' ASP ' in aa_seq[i]:
        D.append(aa_seq[i])

E=[]
for i in range(len(aa_seq)):
    if ' GLU ' in aa_seq[i]:
        E.append(aa_seq[i])

## Handle arithmetic pitfalls, i.e 
## if no ASP or GLU residues in sequence, D:E = 0
## if no GLU residues, then define D:E to have value of 0.999

assert (len(D)==0 and len(E)==0) == False, "WARNING: No ASP or GLU \
                                              present in sequence"
if (len(D)==0 and len(E)==0):
    DE_ratio = float(0)
else:
    try:
        DE_ratio = ( float(len(D)) / float(len(E)) )
    except ZeroDivisionError:
        DE_ratio = float(.999)
        

## Compute  electrostatic properties of sequence
Kfrac = ( float(len(K)) / float(len(aa_seq)) )
Rfrac = ( float(len(R)) / float(len(aa_seq)) )
Dfrac = ( float(len(D)) / float(len(aa_seq)) )
Efrac = ( float(len(E)) / float(len(aa_seq)) )
net_charge = Kfrac + Rfrac - Dfrac - Efrac
naa = len(aa_seq)


## Solubility sequence-level predictors (Price et al., 2011)
fracnumcharge = float(len(K) + len(R) + len(D) + len(E)) / float(naa)
fracabsnetcharge = float(abs(len(K) + len(R) - len(D) - len(E))) / float(naa)
fracnetcharge = float(len(K) + len(R) - len(D) - len(E)) / float(naa) # normalized net charge


## Convert 3-letter AA codes into 1-letter AA codes, e.g. 'ARG' => 'R'
## Find 3-letter codes in ATOM fields of PDB file using basic regular expression
import re
sequence=[]
for i in range(len(aa_seq)):
    pattern = r'\b[A-Z]{3}\b' # pattern = r'\s[A-Z]{3}\s'
    sequence.append(re.findall(pattern, aa_seq[i]))

for i in range(len(sequence)):
    sequence[i] = str(sequence[i])
    sequence[i] = sequence[i].strip("['']").strip("''")

assert len(sequence) == len(aa_seq), 'ERROR: Sequence length incorrect'
sequence = "".join(sequence)
assert len(sequence) == len(aa_seq) * 3, 'ERROR: Sequence length incorrect'

aacodes = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
           'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
           'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
           'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


def convert_sequence(sequence):
    """ Accepts a protein sequence in contiguous 3-letter code
        format and converts it to 1-letter code format, e.g.
        'LYSARGHISASPGLU' ==> 'KRHDE'
        INPUT: STRING, OUTPUT: STRING """
    if (len(sequence) % 3) != 0:
        raise ValueError("ERROR: Input sequence length must be a multiple of 3")
    else:
        pass
    abridged_seq = ""
    for i in range((len(sequence)/3)):
            abridged_seq += aacodes[sequence[(3*i):(3*i+3)]]
    assert len(abridged_seq) == len(aa_seq), "ERROR: Sequence Conversion Incorrect"
    return abridged_seq

final_seq = convert_sequence(sequence)


## Verify that amino acid frequences and proportions make sense
def validate_protein_seq(seq):
	""" Verify that AA sequence length is right """
	return len(seq) == (seq.count('A') + seq.count('V') + seq.count('I') + seq.count('M') + \
					    seq.count('C') + seq.count('F') + seq.count('W') + seq.count('R') + \
					    seq.count('S') + seq.count('T') + seq.count('N') + seq.count('Q') + \
					    seq.count('P') + seq.count('K') + seq.count('H') + seq.count('D') + \
					    seq.count('E') + seq.count('L') + seq.count('Y') + seq.count('G'))

assert validate_protein_seq(final_seq) == True, 'WARNING: \
AA frequencies do not add up to protein sequence length'


## Verify that protein sequence has 20 distinct residue types
def aa_comp(seq):
    """ Accepts a protein sequence and verifies it has 20 distinct
        AA types. If not, it reports how many types and returns the
        actual types """
    if len(set(seq))==20:
        print "Protein Sequence has all 20 Amino Acid Types"
    elif len(set(seq)) != 20:
        print "WARNING: protein sequence does not contain 20 distinct amino acids"
        print "sequence contains", len(set(seq)), "distinct AA types: ", '\n'
    return

aa_comp(final_seq)


# Compute sequence composition properties
Afrac = float( final_seq.count('A') ) / float(naa)
Cfrac = float( final_seq.count('C') ) / float(naa)
Ffrac = float( final_seq.count('F') ) / float(naa)
Gfrac = float( final_seq.count('G') ) / float(naa)
Hfrac = float( final_seq.count('H') ) / float(naa)
Ifrac = float( final_seq.count('I') ) / float(naa)
Lfrac = float( final_seq.count('L') ) / float(naa)
Mfrac = float( final_seq.count('M') ) / float(naa)
Nfrac = float( final_seq.count('N') ) / float(naa)
Pfrac = float( final_seq.count('P') ) / float(naa)
Qfrac = float( final_seq.count('Q') ) / float(naa)
Sfrac = float( final_seq.count('S') ) / float(naa)
Tfrac = float( final_seq.count('T') ) / float(naa)
Wfrac = float( final_seq.count('W') ) / float(naa)
Yfrac = float( final_seq.count('Y') ) / float(naa)
Vfrac = float( final_seq.count('V') ) / float(naa)

def validate_composition(seq):
    """ Verify that AA composition percentages sum up to a value that is
        very near 100 % """
    return float(.998) <= (Yfrac + Wfrac + Tfrac + Sfrac + Qfrac + Pfrac + \
                        Nfrac + Mfrac + Lfrac + Ifrac + Gfrac + Ffrac + \
                        Vfrac + Cfrac + Afrac + Efrac + Dfrac + Rfrac + \
                        Kfrac + Hfrac)

assert validate_composition(final_seq) == True, "WARNING: \
AA proportional compositions do not sum to 1.0 (less than 100%)"

## Search for low-complexity sequence regions
def find_low_complexity_regions(seq):
    #low_compl_reg=[]
    for i in set(final_seq):
        if i*4 in final_seq:
            pass
            #low_compl_reg.append(final_seq[i])

## Print protein name and primary sequence length
print argv[1], 'SUMMARY'
print '\n', 'Protein Name:', title
print '\n', 'Structure Determination Methodology:', expdta
print '\n', 'Protein Sequence:', '\n', final_seq
print '\n', 'naa', '=', naa

## Print KR-ratio and DE-ratio values
print '\n', 'Lysine Count', '=', final_seq.count('K')
print '\n', 'Arginine Count', '=', final_seq.count('R')
print '\n', 'K:R', '=', KR_ratio
print '\n', 'D:E', '=', DE_ratio

'''
## Print proportional composition of ionizable side chains
print '\n', 'Kfrac', '=', Kfrac
print '\n', 'Rfrac', '=', Rfrac
print '\n', 'Dfrac', '=', Dfrac
print '\n', 'Efrac', '=', Efrac
print '\n', 'Hfrac', '=', Hfrac

## Print proportional composition of all other side chains
print '\n', 'Afrac', '=', Afrac
print '\n', 'Cfrac', '=', Cfrac
print '\n', 'Vfrac', '=', Vfrac
print '\n', 'Ffrac', '=', Ffrac
print '\n', 'Gfrac', '=', Gfrac
print '\n', 'Ifrac', '=', Ifrac
print '\n', 'Lfrac', '=', Lfrac
print '\n', 'Mfrac', '=', Mfrac
print '\n', 'Nfrac', '=', Nfrac
print '\n', 'Pfrac', '=', Pfrac
print '\n', 'Qfrac', '=', Qfrac
print '\n', 'Sfrac', '=', Sfrac
print '\n', 'Tfrac', '=', Tfrac
print '\n', 'Wfrac', '=', Wfrac
print '\n', 'Yfrac', '=', Yfrac
'''

print '\n', 'net charge (Kfrac+Rfrac-Dfrac-Efrac)', '=', net_charge
print '\n', 'fraction of charged residues', '=', fracnumcharge
print '\n', 'fractional absolute net charge', '=', fracabsnetcharge
print '\n', 'fractional net charge', '=', fracnetcharge


## Write data out to file if user so desires 
import csv

def export_to_csv():
    """ Export data to csv file """
    name = argv[1][0:4]
    b = open(name + '_summary.csv', 'w')
    a = csv.writer(b)
    ## Specify output format
    column_titles = [[ 'PDB-ID', 'naa', 'K:R', 'D:E', 'Afrac', 'Cfrac', 'Dfrac', 'Efrac', 'Kfrac', 'Gfrac', \
                                'Hfrac', 'Ifrac',  'Kfrac', 'Lfrac', 'Mfrac', 'Nfrac', 'Pfrac', 'Qfrac', 'Rfrac', \
                                'Sfrac', 'Tfrac', 'Vfrac', 'Wfrac', 'Yfrac', 'net charge', \
                                'fracnumcharge', 'fracabsnetcharge', 'fracnetcharge' ]]
                                
    data = [[ name, naa, KR_ratio, DE_ratio, Afrac, Cfrac, Dfrac, Dfrac, Efrac, Gfrac, Hfrac, \
                    Ifrac, Kfrac, Lfrac, Mfrac, Nfrac, Pfrac, Qfrac, Rfrac, Sfrac, Tfrac, Vfrac, Wfrac, \
                    Yfrac, net_charge, fracnumcharge, fracabsnetcharge, fracnetcharge ]]
    ## Write to file 
    a.writerows(column_titles)
    a.writerows(data)
    b.close()
    return
    

def user_choice():
        """ Ask user if they want to export results to csv file """
        choice = raw_input("\nExport data to csv file? (Y/N): ")
        if choice in ("Y", "y"):
            export_to_csv()
        elif choice in("N", "n"):
            pass
        else:
            print "\nYour input was invalid, type either Y/y or N/n"
        return

user_choice()
