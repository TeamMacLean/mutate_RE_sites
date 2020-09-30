#!/usr/bin/env python

import os, sys
import argparse
from Bio import SeqIO
import time

Description="Program to recognise restriction enzyme sites and mutate or replace a nucleotide in the restriction site with another nucleotide such that it does not change the amino acid and the new codon has the highest score of occurrence as shown in the codon file"
usage="""
python {script} --input seq.fasta --re list_of_re_sites.txt --codon list_of_codons_with_scores.txt --output cleaned.fasta

""".format(script=sys.argv[0])

parser=argparse.ArgumentParser(description=Description, epilog=usage)
parser.add_argument("-i", "--input", action="store", dest="input", help="Fasta input file")
parser.add_argument("--re", action="store", dest="enzymes", help="Restriction enzymes in column 1 with restriction site in column 2")
parser.add_argument("--codon", action="store", dest="codonfile", help="A file with codon with occurrence score")
parser.add_argument("--out", action="store", dest="output", default="output.fasta", help="Output file name")

options=parser.parse_args()


if not options.codonfile:
    #print "Please provide a codon file showing the occurrences of each codon"
    exit(1)
elif not options.input:
    #print "Please provide an input fasta file with nucleotide sequences"
    exit(1)
elif not options.enzymes:
    #print "Please provide a list of restriction enzymes in column 1 and restriction sites in column 2"
    exit(1)
else:
    pass


def reverse_complement(ntseq):

    #reverse the ntseq
    ntseq=ntseq[::-1]
    from_this="actgnACTGN"
    to_this="tgacnTGACN"

    #complement the ntseq and return
    return ntseq.translate(str.maketrans(from_this, to_this))



#read the restriction enzyme sites

restriction_sites=[]
restriction_sites_lengths=[]
fh=open(options.enzymes)
for line in fh:
    line=line.rstrip()
    if line=="": continue
    linearray=line.split()
    restriction_sites.append(linearray[1].upper())
    restriction_sites.append(reverse_complement(linearray[1]).upper())
    restriction_sites_lengths.append(len(linearray[1]))

fh.close()
#print(restriction_sites)
#print(restriction_sites_lengths)

min_cut_site_length=min(restriction_sites_lengths)
max_cut_site_length=max(restriction_sites_lengths)



#read the codon file first
aminoacids=dict()
codons=dict()
codon_nucleotides={}
codon_scores={}
codonfile=open(options.codonfile)
for line in codonfile:
    line=line.rstrip()
    if line=="": continue

    linearray=line.split()
    aa=linearray[0]
    codon=linearray[1]
    value=float(linearray[2])

    if aa in aminoacids.keys():
        aminoacids[aa].append(codon)
    else:
        aminoacids[aa]=[codon]

    if aa in codons.keys():
        codons[aa].update({codon:value})
    else:
        codons[aa]={codon:value}

    codon_nucleotides[codon]=aa
    codon_scores[codon]=value

codonfile.close()

#print codons
#print codon_nucleotides

def get_codon_with_highest_score(aminoacid):
    score=0
    aminoacid_codons=aminoacids[aminoacid]
    for codon in aminoacid_codons:
        if codon_scores[codon] >= score:
            score=codon_scores[codon]

    return score

def get_score_for_codon(codon):
    return codon_scores[codon]

def get_codons_for_aminoacid(aminoacid):
    return aminoacids[aminoacid]

def get_aminoacid_for_codon(codon):
    for aminoacid in aminoacids.keys():
        if codon in aminoacids[aminoacid]:
            return aminoacid
    return None

def change_cut_site(subseq):

    if subseq in restriction_sites:
        pass

    #check where is the alternate codon with highest score

    high_scoring_codons={}
    #for origcodon in (subseq[:3], subseq[3:6], subseq[6:]):
    for origcodon in map(''.join, zip(*[iter(subseq)]*3)):
        aminoacid=codon_nucleotides[origcodon]  # gets aa for that codon
        codons_for_aminoacid=get_codons_for_aminoacid(aminoacid)
        #print("codons for aminoacid ", aminoacid, " are ", codons_for_aminoacid)
        #codons_for_aminoacid.delete(codon) # remove the existing codon from the list
        highest_score=0  #initiate the highest score as 0
        for aa_codon in codons_for_aminoacid:
            if aa_codon == origcodon:    # if the codon is same as original codon, don't add
                continue
            else:
                if aa_codon[0:2] == origcodon[0:2]:    #every amino acid has the first two Ns same
                    codon_score=get_score_for_codon(aa_codon)
                    if codon_score >= highest_score:     #if existing codon has the highest score, it will also be added
                        highest_score=codon_score
                        high_scoring_codons[aa_codon]={"score":highest_score, "original_codon":origcodon}


    #print("Highest scoring codons ", high_scoring_codons)
    #lets find the highest scoring codon and its score
    highest_scoring_codon=""
    highest_score=0
    original_codon_at_highest_score=""
    for codon in high_scoring_codons.keys():
        if codon == high_scoring_codons[codon]["original_codon"]:   #if original codon has the highest score, just skip it
            continue

        if high_scoring_codons[codon]["score"] > highest_score:
            highest_score=high_scoring_codons[codon]["score"]
            highest_scoring_codon=codon
            original_codon_at_highest_score=high_scoring_codons[codon]["original_codon"]


    #print ("highest scoring codon ", highest_scoring_codon, "highest score ", highest_score, "original codon ", original_codon_at_highest_score)
    #now lets reconstruct the sequence with highest scoring codon
    ntseq=""
    has_changed=False
    #for codon in list_of_codons_in_frame:
    for codon in map(''.join, zip(*[iter(subseq)]*3)):
        #print("codon in frame",codon)
        if codon == original_codon_at_highest_score and has_changed==False:
            #print("Changing codon ", codon, " to ", highest_scoring_codon)
            ntseq+=highest_scoring_codon
            has_changed=True
        else:
            ntseq+=codon

    # now add trailing sequence

    #ntseq=ntseq+bases_left_at_right

    #print("Returning after re mutation", ntseq)
    return ntseq

    #highest_scoring_codon=sorted(high_scoring_codons["score"], key=high_scoring_codons["score"].__getitem__, reverse=True)[0]
    #highest_score=sorted(high_scoring_codons.values())[0]

def is_re_at_the_end(sequence):
    for seq in restriction_sites:
        if sequence.endswith(seq):
            return True

    return False

def is_re_at_the_start(sequence):
    for seq in restriction_sites:
        if sequence.startswith(seq):
            return True
    return False

def look_in_re_sites(sequence):
    for seq in restriction_sites:
        if sequence ==  seq:
            return True, seq
        elif seq in sequence:
            return True, seq
    return False, None

def re_site_start(resite,sequence):
    return sequence.find(resite)


def find_re_and_replace_codon(ntseq):

    #loop through the RE sites and find the restriction sites
    ntseq_length = len(ntseq)
    new_ntseq = ""

    position=0
    while position <= ntseq_length - 1:
        subseq=ntseq[position:position + 9]
        found, re_site = look_in_re_sites(subseq)

        if found == True:
            # print("subseq ", subseq)
            re_start_point=re_site_start(re_site, subseq) + 1
            re_end_point = re_start_point + len(re_site) -1
            # print("Found restriction site ", re_site)
            # print("RE start point ", re_start_point, " and end point ", re_end_point)

            if re_start_point <=3:
                subseq_to_change = subseq[:6]
                after_change = change_cut_site(subseq_to_change)
                new_ntseq += subseq.replace(subseq[:6], after_change)
                # print("before change ", "subseq to change ", subseq_to_change, subseq, "After change ", after_change)
            else:
                subseq_to_change = subseq[3:]
                after_change = change_cut_site(subseq_to_change)
                new_ntseq += subseq.replace(subseq[3:], after_change)

                # print("before change ", subseq, "subseq to change ", subseq_to_change, "After change ", after_change)
            
            
            position=len(new_ntseq)

        else:
            new_ntseq+=ntseq[position:position+3]
            position +=3

    return new_ntseq


#now read fasta input sequence file
fastafile=open(options.input)
outfilehandler=open(options.output, "w")

##print restriction_sites
for rec in SeqIO.parse(fastafile, "fasta"):
    seqid=rec.id
    ntseq=str(rec.seq).upper()

    new_ntseq=find_re_and_replace_codon(ntseq)
    ##print "Before ", ntseq
    ##print ">" + seqid
    ##print new_ntseq
    outfilehandler.write(">" + seqid + "\n")
    outfilehandler.write(new_ntseq + "\n")

fastafile.close()
outfilehandler.close()
exit(0)
