import fileinput
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable

def append_reads(read1_file, read2_file, population):
    reads2 = SeqIO.index(read2_file, 'fastq')
    rec1_iterator = SeqIO.parse(read1_file, 'fastq')
    with open(population + '_reads.fastq', 'w+') as fout:
        for rec1 in rec1_iterator:
            if rec1.id in reads2:
                next_read = SeqRecord(rec1.seq + reads2[rec1.id].seq.reverse_complement(), id=rec1.id, letter_annotations={"phred_quality":rec1.letter_annotations["phred_quality"] + reads2[rec1.id].letter_annotations["phred_quality"][::-1]})
                SeqIO.write(next_read, fout, 'fastq')
    return

def merge_files(population, read, *files):
    with open(population + '_' + read + '.fastq', 'w+') as fout:
        fin = fileinput.input(files=files)
        for line in fin:
            fout.write(line)
        fin.close()
    return

def translate_reads(read_file, population):
    record_iterator = SeqIO.parse(read_file, 'fastq')
    with open(population + '_residues.fastq', 'w+') as fout:
        for rec in record_iterator:
            translation = rec.seq[13:19].translate() + rec.seq[39:42].translate() + rec.seq[48:54].translate()
            line = str(translation)
            fout.write(line + '\n')
    return

def trim_qc(read1_file, read2_file, population, r1_start = 'CGTCCTGG', r1_end = 'CCATCTTCA', r1_len = 32, r2_start = 'GCAGTTGA', r2_end = 'GCAGAGG', r2_len = 34, min_Qscore = 20):
    # Trims Illumina R1 and R2 reads given the length and sequences for the 5'- and 3'-end of the desired trimmed reads.
    # Filters Illumina R1 and R2 files for a given population with a minimum Q score
    rec1_iterator = SeqIO.parse(read1_file, "fastq")
    rec2_iterator = SeqIO.parse(read2_file, "fastq")
    with open(population + '_R1_trim_qc.fastq', 'w+') as fout:
        for rec1 in rec1_iterator:
            if rec1.seq.count(r1_start) == 1:
                if rec1.seq.count(r1_end) == 1:
                    search_start = rec1.seq.find(r1_start)
                    search_end = rec1.seq.find(r1_end)
                    if len(rec1.seq[search_start:search_end+len(r1_end)]) == r1_len:
                        if min(rec1.letter_annotations['phred_quality'][search_start:search_end+len(r1_end)]) >= min_Qscore:
                            SeqIO.write(rec1[search_start:search_end+len(r1_end)], fout, 'fastq')

    with open(population + '_R2_trim_qc.fastq', 'w+') as fout:
        for rec2 in rec2_iterator:
            if rec2.seq.count(r2_start) == 1:
                if rec2.seq.count(r2_end) == 1:
                    search_start = rec2.seq.find(r2_start[:7])
                    search_end = rec2.seq.find(r2_end)
                    if len(rec2.seq[search_start:search_end+len(r2_end)]) == r2_len:
                        if min(rec2.letter_annotations['phred_quality'][search_start:search_end+len(r2_end)]) >= min_Qscore:
                            SeqIO.write(rec2[search_start:search_end+len(r2_end)], fout, 'fastq')
    return
