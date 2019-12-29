import itertools
from itertools import chain
from itertools import combinations
import pandas as pd
import numpy as np
import networkx as nx
from decimal import Decimal
from decimal import *
import math
from math import log10
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from CAI import CAI
import heapq

def build_DiGraph(data_file):
# Generate a list of nodes and edges with which to generate a directed force-directed graph
    def hamming(var1, var2):
        return sum(aa1 != aa2 for aa1, aa2 in zip(var1, var2))

    g = nx.DiGraph()
    data = pd.read_csv(data_file)
    # Note that in the .csv file, the variants must be ordered by enrichment rate for the edge directions to be consistent
    var_list = data["variant"]
    enrich_list = data["Enrichment Rate"]
    k = 0
    for var in var_list:
        g.add_node(var, enrichment=enrich_list[k])
        k += 1
    for v1, v2 in combinations(var_list,2):
        if '*' not in str(v1) and '*' not in str(v2):
            if hamming(str(v1), str(v2)) == 1:
                # create edge between v1 and v2
                g.add_edge(str(v1), str(v2))
    return nx.write_graphml(g, "FDG_Directed.graphml")

def calc_CAI(read1_file, read2_file, population2, min_count = 5, top = 100, bottom = 100):
    # Iterates through reads for a given population and calculates the CAI for each unique read
    # Analysis is restricted to the top/bottom user-defined number of reads for a given population
    # The weights below are those for yeast in Sharp et al. These will be used in the CAI calculation
    weights = {"TTT": 0.113, "TTC": 1, "TTA": .117, "TTG": 1, "CTT": .006, "CTC": .003, "CTA": .039, "CTG": .003, "ATT": .823, "ATC": 1, "ATA": .003, "ATG": 1,
    "GTT": 1, "GTC": .831, "GTA": .002, "GTG": .018, "TCT": 1, "TCC": .693, "TCA": .036, "TCG": .005, "CCT": .047, "CCC": .009, "CCA": 1, "CCG": .002, "ACT": .921, "ACC": 1,
    "ACA": .012, "ACG": .006, "GCT": 1, "GCC": .316, "GCA": .015, "GCG": .001, "TAT": .071, "TAC": 1, "CAT": .245, "CAC": 1, "CAA": 1, "CAG": .007, "AAT": .053, "AAC": 1, "AAA": .135,
    "AAG": 1, "GAT": .554, "GAC":1, "GAA": 1, "GAG": .016, "TGT": 1, "TGC": .077, "TGG": 1, "CGT": .137, "CGC": .002, "CGA": .002, "CGG": .002, "AGT": .021, "AGC": .031, "AGA": 1,
    "AGG": .003, "GGT": 1, "GGC": .02, "GGA": .002, "GGG": .004}

    # rec2 should be a later screen step (e.g. rec1 = naive, rec2 = PS4)
    rec1_iterator = SeqIO.parse(read1_file,"fastq")
    rec2_iterator = SeqIO.parse(read2_file,"fastq")

    reads1 = {}
    for read in rec1_iterator:
        reads1[str(read.seq[13:19] + read.seq[39:42] + read.seq[48:54])] = reads1.get(str(read.seq[13:19] + read.seq[39:42] + read.seq[48:54]), 0) + 1

    reads2 = {}
    for read in rec2_iterator:
        reads2[str(read.seq[13:19] + read.seq[39:42] + read.seq[48:54])] = reads2.get(str(read.seq[13:19] + read.seq[39:42] + read.seq[48:54]), 0) + 1

    reads2_enrichments = {}
    for key in reads2.keys():
        if reads1.get(key) == None:
            pass
        elif Decimal(reads1[key]) >= min_count:
            reads2_enrichments[key] = Decimal(reads2[key])/Decimal(reads1[key])

    top_reads2 = heapq.nlargest(top, reads2_enrichments, key=reads2_enrichments.get)
    bottom_reads2 = heapq.nsmallest(bottom, reads2_enrichments, key=reads2_enrichments.get)

    with open(population2 + '_CAI_v_enrich.csv', 'w+') as fout:
        fout.write('Top ' + str(top) + ' reads' + '\n')
        fout.write('DNA Sequence' + ',' + 'AA Sequence' + ',' + 'Enrichment Rate' + ',' + 'CAI' + '\n')
        for var in top_reads2:
            cai = CAI(var, weights=weights)
            translation = Seq(var, generic_dna).translate()
            enrichment = reads2_enrichments[var]
            list = var + ',' + str(translation) + ',' + str(enrichment) + ',' + str(cai) + '\n'
            fout.write(list)
        fout.write('Bottom ' + str(bottom) + ' reads' + '\n')
        fout.write('DNA Sequence' + ',' + 'AA Sequence' + ',' + 'Enrichment Rate' + ',' + 'CAI' + '\n')
        for var in bottom_reads2:
            cai = CAI(var, weights=weights)
            translation = Seq(var, generic_dna).translate()
            enrichment = reads2_enrichments[var]
            list = var + ',' + str(translation) + ',' + str(enrichment) + ',' + str(cai) + '\n'
            fout.write(list)
    return

def calc_enrich(file1, file2, population1, population2, min_count):
    aa = 'AILVMFWYNCQSTDERHKGP*'
    population1_counts = {}
    population2_counts = {}
    for var in map(''.join, itertools.product(aa, repeat=5)):
        population1_counts[var] = 0
        population2_counts[var] = 0

    with open(file1) as fin:
        for line in fin:
            read = str(line).strip()
            population1_counts[read] += 1

    with open(file2) as fin:
        for line in fin:
            read = str(line).strip()
            population2_counts[read] += 1

    global population2_enrichment
    population2_enrichment = {}
    for key in population1_counts.keys():
        if Decimal(population1_counts[key]) >= min_count:
            population2_enrichment[key] = Decimal(population2_counts[key])/Decimal(population1_counts[key])
    return

def print_enrich(file1, file2, population1, population2, min_count = 15):
    from analysis import calc_enrich

    calc_enrich(file1, file2, population1, population2, min_count)

    with open(population2 + '_over_' + population1 + '_enrichment.csv', 'w+') as fout:
        fout.write('variant' + ',' + 'enrichment rate' + '\n')
        for key in population2_enrichment.keys():
            fout.write(key + ',' + str(population2_enrichment[key]) + '\n')
    return

def count_miscalls(file1, file2, population):
    countfor = 0
    miscallfor = 0
    corrfor = 0
    wt = 'TCCTCAGCT'
    w = range(8,13)
    record1_iterator = SeqIO.parse(file1,"fastq")
    for rec1 in record1_iterator:
        countfor += 1
        for k in chain(w):
            if rec1.seq[k] != wt[k-8]:
                miscallfor += 1
            elif rec1.seq[k] == wt[k-8]:
                corrfor += 1

    x = range(19,23)
    record1_iterator = SeqIO.parse(file1,"fastq")
    for rec1 in record1_iterator:
        for l in chain(x):
            if rec1.seq[l] != wt[l-14]:
                miscallfor += 1
            elif rec1.seq[l] == wt[l-14]:
                corrfor += 1

    countrev = 0
    miscallrev = 0
    corrrev = 0
    wt = 'TGATGGGCAG'
    y = range(8,12)
    record2_iterator = SeqIO.parse(file2,"fastq")
    for rec2 in record2_iterator:
        countrev += 1
        for m in chain(y):
            if rec2.seq[m] != wt[m-8]:
                miscallrev += 1
            elif rec2.seq[m] == wt[m-8]:
                corrrev += 1

    z = range(18,24)
    record2_iterator = SeqIO.parse(file2,"fastq")
    for rec2 in record2_iterator:
        for n in chain(z):
            if rec2.seq[n] != wt[n-14]:
                miscallrev += 1
            elif rec2.seq[n] == wt[n-14]:
                corrrev += 1

    with open("miscalls_" + population + ".csv","w+") as fout:
        fout.write("Forward Reads" + "\n")
        fout.write("No. Miscalls" + "," + "No. Correct" + "\n")
        fout.write(str(miscallfor) + "," + str(corrfor) + "\n")
        fout.write("Reverse Reads" + "\n")
        fout.write("No. Miscalls" + "," + "No. Correct" + "\n")
        fout.write(str(miscallrev) + "," + str(corrrev) + "\n")
    return

def count_Stops(residue_file, population):
    total_reads = 0
    stop_count = 0
    stop_dic = {}
    with open(residue_file) as fin:
        for line in fin:
            read = str(line).strip()
            total_reads += 1
            if "*" in read:
                stop_count += 1
                stop_dic[read] = stop_dic.get(read, 0) + 1
    with open(population + " Stop Codon Count.csv","w+") as fout:
        fout.write(population + " Stop Codon Count" + "," + population + " Total Reads" + "," + "Fraction of Total Reads" + "\n")
        fout.write(str(stop_count) + "," + str(total_reads) + "," + str(stop_count/total_reads) + "\n")
        fout.write("Variant" + "," + "log2(Count)" + "," + "Library" + "\n")
        for key in stop_dic.keys():
            list = key + "," + str(np.log2(stop_dic[key])) + "," + population + "\n"
            fout.write(list)
    return

def enrich_histogram(file1, file2, population1, population2, min_count = 15):
    from analysis import calc_enrich

    calc_enrich(file1, file2, population1, population2, min_count)

    hist_dic = {}
    for value in population2_enrichment.values():
        hist_dic[str(format(value, '.2g'))] = hist_dic.get(str(format(value, '.2g')), 0) + 1

    with open(population2 + "_over_" + population1 + "_hist.csv","w+") as fout:
        fout.write("Enrichment Rate" + "," + "Log10(Count)" + "," + "Count" + "\n")
        for key in hist_dic.keys():
            list = key + "," + str(log10(hist_dic[key])) + "," + str(hist_dic[key]) + "\n"
            fout.write(list)
    return

def single_enrich(file1, file2, min_count = 15, top = 25):
    aa = 'AILVMFWYNCQSTDERHKGP*'
    population1_counts = {}
    population2_counts = {}
    for var in map(''.join, itertools.product(aa, repeat=5)):
        population1_counts[var] = 0
        population2_counts[var] = 0

    with open(file1) as fin:
        for line in fin:
            read = str(line).strip()
            population1_counts[read] += 1

    with open(file2) as fin:
        for line in fin:
            read = str(line).strip()
            population2_counts[read] += 1

    population2_enrichment = {}
    for key in population1_counts.keys():
        if Decimal(population1_counts[key]) >= min_count:
            population2_enrichment[key] = Decimal(population2_counts[key])/Decimal(population1_counts[key])

    population2_enrichment_2 = {}
    for key in population2_counts.keys():
        if Decimal(population1_counts[key]) > 0:
            population2_enrichment_2[key] = Decimal(population2_counts[key])/Decimal(population1_counts[key])

    top_pop2 = heapq.nlargest(top, population2_enrichment, key=population2_enrichment.get)

    with open("single_enrich.csv","w+") as fout:
        fout.write("variant" + "," + "Enrichment Rate" + "," + "Single Residue Enrichment Rate" + "\n")
        k = 1
        for var in top_pop2:
            fout.write(var + "," + str(population2_enrichment[var]) + ",")
            if var[0] != "T":
                fout.write("T88" + var[0] + "," + str(population2_enrichment_2[var[0] + "QWLH"]) + "\n")
            if var[1] != "Q":
                fout.write("," + "," + "Q89" + var[1] + "," + str(population2_enrichment_2["T" + var[1] + "WLH"]) + "\n")
            if var[2] != "W":
                fout.write("," + "," + "W246" + var[2] + "," + str(population2_enrichment_2["TQ" + var[2] + "LH"]) + "\n")
            if var[3] != "L":
                fout.write("," + "," + "L249" + var[3] + "," + str(population2_enrichment_2["TQW" + var[3] + "H"]) + "\n")
            if var[4] != "H":
                fout.write("," + "," + "H250" + var[4] + "," + str(population2_enrichment_2["TQWL" + var[4]]) + "\n")
    return

def variant_enrich(variant, file1, file2, population1, min_count = 15):
    aa = 'AILVMFWYNCQSTDERHKGP*'
    population1_counts = {}
    population2_counts = {}
    for var in map(''.join, itertools.product(aa, repeat=5)):
        population1_counts[var] = 0
        population2_counts[var] = 0

    with open(file1) as fin:
        for line in fin:
            read = str(line).strip()
            population1_counts[read] += 1

    with open(file2) as fin:
        for line in fin:
            read = str(line).strip()
            population2_counts[read] += 1

    enrichment_rate = Decimal(population2_counts[variant])/Decimal(population1_counts[variant])
    print('The enrichment rate of variant {} is {:.5f}'.format(variant, enrichment_rate))

    if population1_counts[variant] < min_count:
        print('Note, {} is present in the {} library {} times'.format(variant, population1, population1_counts[variant]))

    return

def qScore_analysis(file1, file2, population):
    getcontext().prec = 5
    countfor = 0
    x = range(13,19)
    errorfor = []
    record1_iterator = SeqIO.parse(file1, "fastq")
    for rec1 in record1_iterator:
        countfor += 1
        for l in chain(x):
            score = rec1.letter_annotations["phred_quality"][l]
            errorfor.append(10**(score/-10))

    avgerrorfor = Decimal(sum(errorfor))/Decimal(len(errorfor))
    avgqfor = Decimal(-10*math.log10(avgerrorfor))

    # Repeating the analysis for the reverse reads
    countrev = 0
    y = range(12,18)
    errorrev = []
    record2_iterator = SeqIO.parse(file2, "fastq")
    for rec2 in record2_iterator:
        countrev += 1
        for m in chain(y):
            score = rec2.letter_annotations["phred_quality"][m]
            errorrev.append(10**(score/-10))

    z = range(24,27)
    record2_iterator = SeqIO.parse(file2,"fastq")
    for rec2 in record2_iterator:
        for n in chain(z):
            score = rec2.letter_annotations["phred_quality"][n]
            errorrev.append(10**(score/-10))

    avgerrorrev = Decimal(sum(errorrev))/Decimal(len(errorrev))
    avgqrev = Decimal(-10*math.log10(avgerrorrev))

    with open("qScore_" + population + ".csv","w+") as fout:
        fout.write("Forward Reads" + "\n")
        fout.write("Average Q Score" + "," + "Error Rate" + "," + "# Reads Analyzed" + "\n")
        fout.write(str(avgqfor) + "," + str(avgerrorfor) + "," + str(countfor) + "\n")
        fout.write("Reverse Reads" + "\n")
        fout.write("Average Q Score" + "," + "Error Rate" + "," + "# Reads Analyzed" + "\n")
        fout.write(str(avgqrev) + "," + str(avgerrorrev) + "," + str(countrev))
    return

def qScore_flanking(file1, file2, population):
    getcontext().prec = 5
    countfor = 0
    errorfor = []
    w = range(8,13)
    record1_iterator = SeqIO.parse(file1,"fastq")
    for rec1 in record1_iterator:
        countfor += 1
        for k in chain(w):
            score = rec1.letter_annotations["phred_quality"][k]
            errorfor.append(10**(score/-10))

    x = range(19,23)
    record1_iterator = SeqIO.parse(file1,"fastq")
    for rec1 in record1_iterator:
        for l in chain(x):
            score = rec1.letter_annotations["phred_quality"][l]
            errorfor.append(10**(score/-10))

    avgerrorfor = Decimal(sum(errorfor))/Decimal(len(errorfor))
    avgqfor = Decimal(-10*math.log10(avgerrorfor))

    # Repeat for reverse read

    countrev = 0
    errorrev = []
    y = range(8,12)
    record2_iterator = SeqIO.parse(file2,"fastq")
    for rec2 in record2_iterator:
        countrev += 1
        for m in chain(y):
            score = rec2.letter_annotations["phred_quality"][m]
            errorrev.append(10**(score/-10))

    z = range(18,24)
    record2_iterator = SeqIO.parse(file2,"fastq")
    for rec2 in record2_iterator:
        for n in chain(z):
            score = rec2.letter_annotations["phred_quality"][n]
            errorrev.append(10**(score/-10))

    avgerrorrev = Decimal(sum(errorrev))/Decimal(len(errorrev))
    avgqrev = Decimal(-10*math.log10(avgerrorrev))

    with open("qScore_flanking_" + population + ".csv","w+") as fout:
        fout.write("Forward Reads" + "\n")
        fout.write("Average Q Score" + "," + "Error Rate" + "," + "# Reads Analyzed" + "\n")
        fout.write(str(avgqfor) + "," + str(avgerrorfor) + "," + str(countfor) + "\n")
        fout.write("Reverse Reads" + "\n")
        fout.write("Average Q Score" + "," + "Error Rate" + "," + "# Reads Analyzed" + "\n")
        fout.write(str(avgqrev) + "," + str(avgerrorrev) + "," + str(countrev))
    return
