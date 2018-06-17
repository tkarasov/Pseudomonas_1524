from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio import SeqIO
from Bio.Seq import Seq
import os
from subprocess import Popen, PIPE
"""The goal of this script is to take a fasta record and identify all fourfold degenerate sites in the record. Then to output the alignments with only these fourfold sites"""



def altcodons(codon, table):
    """List codons that code for the same aminonacid / are also stop.

    @param codon
    @table code table id
    @return list of codons

    """
    tab = unambiguous_dna_by_id[table]

    if codon in tab.stop_codons:
        return tab.stop_codons

    try:
        aa = tab.forward_table[codon]
    except:
        return []

    return [k for (k, v) in tab.forward_table.iteritems()
            if v == aa and k[0] == codon[0] and k[1] == codon[1]]


def degeneration(codon, table):
    """Determine how many codons code for the same amino acid / are also stop

    @param codon the codon
    @param table code table id
    @param the number of codons also coding for the amino acid codon codes for

    """
    return len(altcodons(codon, table))


def is_x_degenerated(x, codon, table):
    """Determine if codon is x-fold degenerated.

    @param codon the codon
    @param table code table id
    @param true if x <= the degeneration of the codon

    """
    return (x <= len(altcodons(codon, table)))


def degenerated_subseq(seq, x, table):
    """Get a subsequence consisting of the x-fold degenerated codons only."""
    data = ""
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3].tostring()
        if is_x_degenerated(x, codon, table):
            data += codon
    return data

def list_degenerated(fasta_record, x, table):
    """first consider only codons that vary at the third site"""
    all_fastarec={}
    for rec in fasta_record:
        data = []
        for i in range(0, len(rec.seq), 3):
            codon = rec.seq[i:i + 3].tostring()
            if is_x_degenerated(x, codon, table):
                data.append(i+3)
        all_fastarec[rec.name]=data
    result = set(all_fastarec[all_fastarec.keys()[0]])
    for s in all_fastarec.keys()[1:]:
        result.intersection_update(all_fastarec[s])
    return(list(result))
    
def only_degenerate(all_concat, degenerated):
    """returns only those sites that are four-fold degenerate and found in all individuals"""
    for rec in all_concat:
        temp=''
        for number in (range(len(rec.seq)
)):
            if number in degenerated:
                temp=temp+rec.seq[number]
            if number not in degenerated:
                temp=temp+"-"
                #rec.seq[number]="A"
        rec.seq=Seq(temp)
        #seq_subset=('').join([rec[num] for num in degenerated])
        #rec.seq=Seq(seq_subset)
    return all_concat
    
    
#ClonalFrameML concat_107.newick /ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/clonalframe/all_concat.fasta concat_107_clonalframe_output


all_concat=list(SeqIO.parse(open("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/clonalframe/concat_107_no_recom.fasta"), "fasta"))
degenerated=list_degenerated(all_concat, 4, 1)
final=only_degenerate(all_concat, degenerated)

SeqIO.write(final, "/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/clonalframe/all_concat_fourfold_only.fasta", "fasta")

test=Popen(["/usr/bin/raxmlHPC-PTHREADS-SSE3",  "-s", "/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/clonalframe/all_concat_fourfold_only.fasta", "-n", "/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/clonalframe/all_concat_fourfold_only_tree.newick", "-m", "GTRGAMMA", "-p",  "344312987"])
#output = test.communicate()[0]


