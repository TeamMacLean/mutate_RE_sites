#!/usr/bin/env python

import os, sys
import argparse
from Bio import SeqIO
from operator import itemgetter


Description="Program to recognise restriction enzyme sites and mutate or replace a nucleotide in the restriction site with another nucleotide such that it does not change the amino acid and the new codon has the highest score of occurrence as shown in the codon file"
usage="""
python {script} --input seq.fasta --re list_of_re_sites.txt --codon list_of_codons_with_scores.txt --output cleaned.fasta

""".format(script=sys.argv[0])

parser=argparse.ArgumentParser(description=Description, epilog=usage)
parser.add_argument("-i", "--input", action="store", dest="input", help="Fasta input file")
parser.add_argument("--re", action="store", dest="enzymes", help="Restriction enzymes in column 1 with restriction site in column 2")
parser.add_argument("--codon", action="store", dest="codonfile", help="A file with codon with occurrence score")
parser.add_argument("--out", action="store", dest="output", default="output.fasta", help="Output file name")
parser.add_argument("--report", action="store", dest="report", default="report.txt", help="Report of changes")

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
    """
    Reverse complement a nucleotide sequence.
    
    Args:
        ntseq (str): Nucleotide sequence
        
    Returns:
        str: Reverse complement of the input sequence
    """
    #reverse the ntseq
    ntseq=ntseq[::-1]
    from_this="actgnACTGN"
    to_this="tgacnTGACN"

    #complement the ntseq and return
    return ntseq.translate(str.maketrans(from_this, to_this))


def get_restriction_site(re_sites):
    """
    Read the restriction enzyme sites from a file.
    
    Args:
        re_sites (str): Path to the file containing restriction enzyme sites
        
    Returns:
        tuple: A tuple containing two lists - restriction_sites and restriction_sites_lengths
    """
    #read the restriction enzyme sites

    restriction_sites=[]
    restriction_sites_lengths=[]
    fh=open(re_sites)
    for line in fh:
        line=line.rstrip()
        if line=="": continue
        linearray=line.split()
        restriction_sites.append(linearray[1].upper())
        restriction_sites.append(reverse_complement(linearray[1]).upper())
        restriction_sites_lengths.append(len(linearray[1]))

    fh.close()
    return restriction_sites, restriction_sites_lengths


def get_codon_scores(codonfile):
    """
    Read the amino acid codons and their occurrence scores from a file.
    
    Args:
        codonfile (str): Path to the file containing codons with occurrence scores
        
    Returns:
        dict: A dictionary with amino acids as keys and dictionaries of codons and their occurrence scores as values
    """
    # read the amino acid codons and their values
    aminoacids_codons_values = dict()
    with open(codonfile) as fcon:
        for line in fcon:
            line = line.rstrip()
            if line == "":
                continue
            else:
                linedata = line.split()
                aa = linedata[0]
                codon = linedata[1]
                value = float(linedata[2])
            
                if aa in aminoacids_codons_values.keys():
                    aminoacids_codons_values[aa].update({codon:value})
                else:
                    aminoacids_codons_values[aa] = {codon: value}
    

    return aminoacids_codons_values

def sort_codons_by_values(aa_codons):
    """
    Sort amino acid codons by occurrence score in descending order.
    
    Args:
        aa_codons (dict): A dictionary of codons and their occurrence scores for a specific amino acid
        
    Returns:
        list: A list of (codon, value) pairs sorted by occurrence score in descending order
    """
    return sorted(aa_codons.items(), key=itemgetter(1), reverse=True)

def get_codons_for_aminoacid(aminoacid):
    """
    Get a list of codons for a specific amino acid.
    
    Args:
        aminoacid (str): Amino acid
    
    Returns:
        list: A list of codons for the given amino acid
    """

    for aa in aminoacid_codons_values.keys():
        if aa == aminoacid:
            return list(aminoacid_codons_values.keys() )
    return None

def get_aminoacid_for_codon(codon):
    """
    Get the amino acid for a specific codon.
    
    Args:
        codon (str): Codon
        
    Returns:
        str: Amino acid corresponding to the given codon
    """

    for aa in aminoacid_codons_values.keys():
        if codon in aminoacid_codons_values[aa].keys():
            return aa
        else:
            continue
    
    return None

def change_cut_site(re_subseq, subseq):

    if re_subseq in restriction_sites:
        pass

    #check where is the alternate codon with highest score

    high_scoring_codons={}
    #for origcodon in (re_subseq[:3], re_subseq[3:6], re_subseq[6:]):
    codon_values_in_re_subseq = list()
    for origcodon in map(''.join, zip(*[iter(re_subseq)]*3)):
        aminoacid=get_aminoacid_for_codon(origcodon)  # gets aa for that codon
        codons_for_aminoacid=get_codons_for_aminoacid(aminoacid)
        
        # get the aminoacid codons by sorted values in descending order
        codon_values_in_re_subseq += sort_codons_by_values(aminoacid_codons_values[aminoacid])

    # sort the codons by values from all aminoacids in re_subseq
    sorted_codon_values_in_re_subseq = sorted(codon_values_in_re_subseq, key=itemgetter(1), reverse=True)

    #remove duplicates from the list, if any
    sorted_codon_values_in_re_subseq = list(dict.fromkeys(sorted_codon_values_in_re_subseq))
    #print(sorted_codon_values_in_re_subseq)
    ntseq = ""
    codon_replaced=False
    for codon in map(''.join, zip(*[iter(re_subseq)]*3)):
        aminoacid=get_aminoacid_for_codon(codon)  # gets aa for that codon
        for codon_in_order, value in sorted_codon_values_in_re_subseq:
            aa_for_codon_in_order = get_aminoacid_for_codon(codon_in_order)
            if codon == codon_in_order:
                continue
            elif aminoacid == aa_for_codon_in_order and codon_replaced==False:
                
                # test if replacing with this codon gives another RE site
                #print("replacing ", codon, " by ", codon_in_order, subseq.replace(codon, codon_in_order))
                found, re_site = look_in_re_sites(subseq.replace(codon, codon_in_order, 1) )
                if found == True:
                    #print("Found RE site again. Replacing with second best codon value")
                    continue
                else:
                    ntseq += codon_in_order
                    codon_replaced=True
                    break

            else:
                #print("Adding nonreplaced codon ", codon)
                ntseq += codon
                break

    return ntseq

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

def get_subseq_after_change(subseq, re_start_point, re_end_point):
    ''' replaces codon in the RE site '''
    if re_start_point <=3:
        subseq_to_change = subseq[:6]
        after_change = change_cut_site(subseq_to_change, subseq)
        # subseq_after_change = subseq.replace(subseq_to_change, after_change)
        
        new_subseq = subseq.replace(subseq_to_change, after_change)
        #print("before change ", "subseq to change ", subseq_to_change, subseq, "After change ", after_change, new_subseq)
    else:
        subseq_to_change = subseq[3:]
        after_change = change_cut_site(subseq_to_change, subseq)
        # subseq_after_change = subseq.replace(subseq_to_change, after_change)
        
        new_subseq = subseq.replace(subseq_to_change, after_change)

        #print("before change ", subseq, "subseq to change ", subseq_to_change, "After change ", after_change, new_subseq)

    return subseq_to_change, after_change, new_subseq

def find_re_and_replace_codon(ntseq):

    #loop through the RE sites and find the restriction sites
    ntseq_length = len(ntseq)
    new_ntseq = ""

    report = []
    position=0
    while position <= ntseq_length - 1:
        subseq=ntseq[position:position + 9]
        found, re_site = look_in_re_sites(subseq)
      
        if found == True:
            #print("subseq ", subseq)
            re_start_point=re_site_start(re_site, subseq) + 1
            re_end_point = re_start_point + len(re_site) -1
            #print("Found restriction site ", re_site)
            #print("RE start point ", re_start_point, " and end point ", re_end_point)

            subseq_to_change, subseq_after_change, new_subseq = get_subseq_after_change(subseq, re_start_point, re_end_point)
            report.append([subseq, new_subseq, position, position +  len(subseq)])
            
            new_ntseq += new_subseq
            position = len(new_ntseq)
        else:
            new_ntseq+=ntseq[position:position+3]
            position +=3
            

    return report, new_ntseq


restriction_sites, restriction_sites_lengths = get_restriction_site(options.enzymes)
aminoacid_codons_values = get_codon_scores(options.codonfile)


min_cut_site_length=min(restriction_sites_lengths)
max_cut_site_length=max(restriction_sites_lengths)

#now read fasta input sequence file
fastafile=open(options.input)
outfilehandler=open(options.output, "w")

reportout = open(options.report, 'w')
reportout.write("SEQID\tRESiteBefore\tRESiteAfter\tSTART\tEND\n")
##print restriction_sites
for rec in SeqIO.parse(fastafile, "fasta"):
    seqid=rec.id
    ntseq=str(rec.seq).upper()

    report, new_ntseq=find_re_and_replace_codon(ntseq)
    ##print "Before ", ntseq
    ##print ">" + seqid
    ##print new_ntseq
    outfilehandler.write(">" + seqid + "\n")
    outfilehandler.write(new_ntseq + "\n")

    for data in report:
        reportout.write(seqid + "\t" + data[0] + "\t" + data[1] + "\t" + str(data[2]) + "\t" + str(data[3]) + "\n")

fastafile.close()
outfilehandler.close()
exit(0)
