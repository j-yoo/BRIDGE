# Bridging non-overlapping reads illuminates high-order epistasis between distal protein sites in a GPCR

## Authors

Justin I Yoo, Patrick S. Daugherty, Michelle A. O'Malley*

## Installation and dependencies

BRIDGE analysis source code was tested using python v3.7.3 on Windows 10 (on PC with 2.4 GHz Intel Core i7 and 8GB RAM). The latest version of python can be downloaded [here.](https://www.python.org/downloads/)
The source code requires the following packages:

* NetworkX
* pandas
* numpy
* matplotlib
* seaborn
* Biopython
* CAI

### Installing packages

A convenient method to install python and most of the required packages is through [Anaconda](https://www.anaconda.com/distribution/), which will download the latest version of python and many data science packages in a single installation. If Anaconda is downloaded, then the only packages that will need to be downloaded are 'Biopython' and 'CAI'.

These packages can be downloaded using 'pip' through the Anaconda prompt:

```
pip install 'package name'
```

For example:

```
pip install biopython
```
```
pip install CAI
```

The expected time required to install python and all packages on a 'normal' desktop computer is under 30 minutes.

### Running the code

Running the code below will generate the data presented in the manuscript.
In the "Example data set" section, instructions are given to run the code on a set of provided test files.
Note that the source code and file inputs have to be in the same folder for the code to run properly.

### Expected run time

The expected run time to process and analyze the manuscript data on a 'normal' desktop computer is 2 - 3 hours.

## Implementing code to process raw read files

The following subsections describe code written to process raw .fastq read files. The relevant code is provided in the 'Manuscript Data' folder. The code is executed directly using the python interpreter.

### Merge read files from different flow cell lanes

Forward (R1) and reverse (R2) read files for each library (e.g. naive, PS1, etc.) were produced from 4 flow cell lanes. In the process_reads module, the `merge_files()` function is used to merge all data into single R1 and R2 files for easier downstream processing.

```
merge_files(population, read, *files)
```

Here, the `population` argument takes the library name (e.g. naive, PS1, etc.), the `read` argument takes text identifying the file as R1 or R2, and the `*files` argument takes an arbitrary number of **unzipped** .fastq files that will be merged into a single file. To generate the merged .fastq files used for this study, run the code for all libraries as shown in the example below for the naive library:

```
import process_reads

process_reads.merge_files('naive', 'R1', 'Naive_R1_1.fastq', 'Naive_R1_2.fastq', 'Naive_R1_3.fastq', 'Naive_R1_4.fastq')
```

The output file will be named population_read.fastq (e.g. naive_R1.fastq).

### Trimming and quality-filtering raw reads

The next step comprises trimming and filtering paired-end reads by their quality scores using the `trim_qc()` function in the process_reads module.

```
trim_qc(read1_file, read2_file, population, r1_start = 'CGTCCTGG', r1_end = 'CCATCTTCA', r1_len = 32, r2_start = 'GCAGTTGA', r2_end = 'GCAGAGG', r2_len = 34, min_Qscore = 20)
```

Here, the `read1_file` and `read2_file` arguments take the respective R1 and R2 .fastq files generated using the `merge_files()` function for a given library (e.g. naive_R1.fastq and naive_R2.fastq). The `population` argument is used to name the output file as 'population_R1_trim_qc.fastq' and 'population_R2_trim_qc.fastq'. The `r1_start` and `r1_end` arguments take the sequences (5' to 3') corresponding to the 5' and 3' ends, respectively, of the trimmed forward read. The flanking sequences are discarded. The `r1_len` argument takes the anticipated length of the trimmed sequence in order to remove reads with insertions or deletions altering the reading frame. Similarly, the `r2_start`, `r2_end`, and `r2_len` arguments take the reverse read 5'-end sequence, 3'-end sequence, and trimmed sequence length. The `min_Qscore` argument sets the quality score threshold used to retain reads of sufficient quality. If no value is given, 20 is used by default. To generate the data presented in the manuscript, run the code for all merged R1 and R2 read files as shown in the example below for the naive library:

```
process_reads.trim_qc('naive_R1.fastq', 'naive_R2.fastq', 'naive')
```

There will be two files generated using this function: population_R1_trim_qc.fastq and population_R2_trim_qc.fastq.

### Append reads with matching flow cell (x,y) coordinates

Trimmed, quality-filtered read pairs are matched by (x,y) coordinates encoded in their sequence id and appended using the `append_reads()` function in the process_reads module. The reverse read is reverse-complemented in this step.

```
append_reads(read1_file, read2_file, population)
```

Here, `read1_file` and `read2_file` are the trimmed and quality-filtered R1 and R2 files for a given population. To generate the data presented in the manuscript, run the code for all trimmed/filtered R1 and R2 read files as shown in the example below for the naive library:

```
process_reads.append_reads('naive_R1_trim_qc.fastq', 'naive_R2_trim_qc.fastq', 'naive')
```

The output file will be named 'population_reads.fastq'.

### Translate mutated codons to amino acid sequence

The mutated codons are sliced (i.e. extracted from the sequence) and translated to their corresponding amino acids using the `translate_reads()` function in the process_reads module.

```
translate_reads(read_file, population)
```

Here, the `read_file` argument takes the .fastq file generated from the `append_reads()` function. To generate the data presented in the manuscript, run the code for all appended read files as shown in the example below for the naive library:

```
process_reads.translate_reads('naive_reads.fastq', 'naive')
```

The output file will be named 'population_residues.fastq'.

## Analysis of processed reads

### Calculate enrichment rates

The enrichment rates of variants across two given populations is calculated using the `print_enrich()` function in the analysis module. This function can be used to generate data to populate Table 1 and to generate data used to create Figure 3a.

```
print_enrich(file1, file2, population1, population2, min_count = 15)
```

Here, the enrichment rate is calculated as the abundance of a given variant in `population2` (e.g. Post Sort 4) divided by its abundance in `population1` (e.g. naive). The `file1` and `file2` arguments take their respective population_residues.fastq files, which must be located in the same folder as the code. The `min_count` argument sets the minimum count threshold in population1 to retain a variant in the analysis. If no value is given, 15 is used by default. To generate the data presented in the manuscript, run the code with the appropriate residue files as shown in the example below for the PS4 and naive libraries:

```
import analysis

analysis.print_enrich('naive_residues.fastq', 'PS4_residues.fastq', 'naive', 'PS4')
```

The output file will be named 'population2_over_population1_enrichment.csv'. This file will contain two columns:

* a unique variant's sliced amino acid sequence
* its respective enrichment rate.

### Generate histograms of PS4/naive and PS1/naive enrichment rates

Data used to create Figure 3a were obtained using the `enrich_histogram()` function in the analysis module.

```
enrich_histogram(file1, file2, population1, population2, min_count = 15)
```

Similar to `print_enrich()`, the enrichment rate is calculated as the abundance of a given variant in `population2` (e.g. Post Sort 4) divided by its abundance in `population1` (e.g. naive). The `file1` and `file2` arguments take their respective population_residues.fastq files, which must be located in the same folder as the code. The `min_count` argument sets the minimum count threshold in population1 to retain a variant in the analysis. If no value is given, 15 is used by default. To generate the data presented in Figure 3a, run the code with PS4/naive and PS1/naive residue files as shown in the example below for the PS4 and naive libraries:

```
import analysis

analysis.enrich_histogram('naive_residues.fastq', 'PS4_residues.fastq', 'naive', 'PS4')
```

The output file will be named 'population2_over_population1_hist.csv'. This file will contain the following information:

* an enrichment rate
* log10(count) of variants with the respective enrichment rate
* the count of variants with the respective enrichment rate

To generate Figure 3a, please see section 'Generating figures', subsection 'Figure 3a'.

### Calculate enrichment rate of specific variant

The enrichment rate of a specific variant can be calculated using `variant_enrich()`. This function can be used to generate data to populate Supplementary Table 2.

```
variant_enrich(variant, file1, file2, population1, min_count = 15)
```

Here, the `variant` argument takes the variant of interest. The `file1` and `file2` arguments take residue.fastq files (e.g. naive_residues.fastq and PS4_residues.fastq files, respectively).  The `population1` argument takes the name of the library comprising file1. The `min_count` argument sets the minimum count threshold in population1 to retain a variant in the analysis. If no value is given, 15 is used by default. If the variant of interest is not in population1 with abundance greater than or equal to min_count, the function will print:

'Note, 'variant' is present in the 'population1' library # times'

As an example, the PS4/naive enrichment rate of variant TAWLH is calculated using the following code:

```
analysis.variant_enrich('TAWLH', 'naive_residues.fastq', 'PS4_residues.fastq', 'naive')
```

### Count stop codons within each population

The number of variants containing stop codons is calculated using the `count_Stops()` function in the analysis module. The data generated from this function are used to create Supplementary Fig. 6.

```
count_Stops(residue_file, population)
```

Here, the `residue_file` argument takes the 'population_residues.fastq' file for a given population and generates a .csv file containing the following information:

* total number of variants in the population containing stop codons
* total number of reads analyzed in the population
* fraction of variants that contain at least 1 stop codon
* a list of unique variants containing at least 1 stop codon and their abundance

To generate the data presented in the manuscript, run the code with each of the residues files as shown in the example below for the naive library:

```
analysis.count_Stops('naive_residues.fastq', 'naive')
```

### Calculate codon adaptation index (CAI)

The codon adaptation index (CAI) is calculated using the `calc_CAI()` function in the analysis module. The data generated from this function are used to create Supplementary Fig. 7, which plots CAI against enrichment rate for the 100 most enriched and 100 most depleted variants.

```
calc_CAI(read1_file, read2_file, population2, min_count = 5, top = 100, bottom = 100)
```

Here, the `read1_file` and `read2_file` arguments take the population_reads.fastq files (generated from `append_reads()`) corresponding to the populations from which enrichment rates will be calculated (e.g. read2_file corresponds to PS4 and read1_file corresponds to naive). The `min_count` argument sets the minimum DNA read count threshold in population1 to retain a variant in the analysis. If no value is given, 5 is used by default. Note this value is lower than the default min_count value used in `print_enrich()` to account for synonymous codons encoding for the same variant. The `top` and `bottom` arguments take the number of variants with the highest (top) or lowest (bottom) enrichment rates to include in the analysis. If no values are given, 100 is used by default for both arguments.

The output file will be named 'population2_CAI_v_enrich.csv' and will contain the following information:

* Variant's sliced DNA and corresponding amino acid sequence
* Variant's enrichment rates
* Variant's CAI

To generate the data presented in the manuscript, run the code with the appropriate residues files as shown in the example below for the naive and PS4 libraries:

```
analysis.calc_CAI('naive_reads.fastq', 'PS4_reads.fastq', 'PS4')
```

### Calculate single substitution enrichment rates

The single substitution enrichment rates used to generate Figure 3c can be calculated using `single_enrich()` in the analysis module.

```
single_enrich(file1, file2, min_count = 15, top = 25)
```

Here, the `file1` and `file2` arguments take the naive_residues.fastq and PS4_residues.fastq files, respectively. The `min_count` argument sets the minimum variant count in file1 to retain a variant in the analysis. If no value is given, 15 is used by default. The `top` argument takes the number of variants with the greatest enrichment rates to include in the analysis. If no value is given, 25 is used by default. The output file will be named 'single_enrich.csv' and will contain the following information:

* Variant's sliced amino acid sequence
* Variant's enrichment rates
* Enrichment rate of each single substitution comprising the respective variant

To generate the data presented in the manuscript, run the code with the appropriate residue files as shown in the example below for the PS4 and naive libraries:

```
analysis.single_enrich('naive_residues.fastq', 'PS4_residues.fastq')
```

### Calculate most highly-enriched variants with proximal single, double, or triple substitutions

The 20 most highly-enriched variants containing proximal single (e.g. XQWLH), double (e.g. XXWLH), or triple (i.e. TQXXX) substitutions (Figure 3b) can be calculated by executing proximal_subs.py in the interpreter. Note 'naive_residues.fastq' and 'PS4_residues.fastq' files must be in the same folder as the proximal_subs.py file. In python v3, this can be done by entering the following text in the interpreter:

```
exec(open('proximal_sub.py').read())
```

The output file will be named 'PS4_proximal_subs.csv' and will contain the following information:

* variant amino acid sequence
* variant enrichment rate
* Population

This file is used by the fig3b.py file described in the 'Generating figures' section to generate Figure 3b.

### Determine average quality score of mutated bases

The average Q scores for the mutated bases (Supplementary Table 3) were calculated using the `qScore_analysis()` function in the analysis module.

```
qScore_analysis(file1, file2, population)
```

Here, the `file1` and `file2` arguments take the R1 and R2 files, respectively, output from the trim_qc() function. The `population` argument takes the population name. To generate the data presented in the manuscript, run the code for all trim_qc read files as shown in the example below for the naive library:

```
analysis.qScore_analysis('naive_R1_trim_qc.fastq', 'naive_R2_trim_qc.fastq', 'naive')
```

The output file will be named 'qScore_population.csv' and will contain the following information for the forward and reverse reads:

* average Q score for the mutated bases
* average error rate for the mutated bases
* total number of reads analyzed

### Determine average quality score of constant, flanking bases

The average Q scores for the constant, flanking bases (Supplementary Table 4) were calculated using the `qScore_flanking()` function in the analysis module.

```
qScore_flanking(file1, file2, population)
```

Here, the `file1` and `file2` arguments take the R1 and R2 files, respectively, output from the trim_qc() function. The `population` argument takes the population name. To generate the data presented in the manuscript, run the code for all trim_qc read files as shown in the example below for the naive library:

```
analysis.qScore_flanking('naive_R1_trim_qc.fastq', 'naive_R2_trim_qc.fastq', 'naive')
```

The output file will be named 'qScore_flanking_population.csv' and will contain the following information for the forward and reverse reads:

* average Q score for the constant bases
* average error rate for the constant bases
* total number of reads analyzed

### Count miscalls in constant, flanking bases

The observed error rates for the constant, flanking bases (Supplementary Table 4) were calculated by counting the number of miscalls at those positions using the `count_miscalls()` function in the analysis module.

```
count_miscalls(file1, file2, population)
```

Here, the `file1` and `file2` arguments take the R1 and R2 files, respectively, output from the trim_qc() function. The `population` argument takes the population name. To generate the data presented in the manuscript, run the code for all trim_qc read files as shown in the example below for the naive library:

```
analysis.count_miscalls('naive_R1_trim_qc.fastq', 'naive_R2_trim_qc.fastq', 'naive')
```

The output file will be named 'miscalls_population.csv' and will contain the following information for the forward and reverse reads:

* number of miscalls
* number of correct calls

### Generate nodes and edges for force-directed graph

A directed .graphml file containing nodes and edges used to generate Figure 5 and Supplementary Fig. 10 was generated using the `build_DiGraph()` function in the analysis module. Here, directionality is imposed such that smaller nodes (i.e. variants with lower enrichment rates) point towards larger nodes (i.e. variants with greater enrichment rates).

```
build_DiGraph(data_file)
```

Here, the `data_file` argument takes an enrichment .csv file. To obtain data presented in the manuscript, this .csv file comprised PS4/naive enrichment rates generated using `analysis.print_enrich('naive_residues.fastq', 'PS4_residues.fastq', 'naive', 'PS4')`. To generate the data presented in the manuscript, run the code using the PS4_over_naive_enrichment.csv file as input as shown in the example below for the naive library:

```
build_DiGraph('PS4_over_naive_enrichment.csv')
```

The output file will be named 'FDG_Directed.graphml' and can be analyzed further using [Gephi](https://gephi.org/) to generate the exact figures presented in Figure 5 and Supplementary Fig. 10.

## Generating figures

The following steps provide instructions for generating figures that were created using code and processed NGS data.

### Figure 3a

1. Use the `enrich_histogram()` function to generate PS4_over_naive_hist.csv and PS1_over_naive_hist.csv files.
2. Move the .csv files into the 'figures' folder
3. Ensure the path directory is set to the figures folder. If the following code is run, the returned file path should match the 'figures' folder in which the 'fig3a.py' code is located.

```
import os
os.getcwd()
```

4. Execute fig3a.py directly in the python interpreter as shown below:

```
exec(open('fig3a.py').read())
```

### Figure 3b

Note that the generation of this figure requires some manual data processing.

1. Run the 'proximal_subs.py' file in the interpreter (as described in the appropriate subsection above).
2. Move the .csv file into the 'figures' folder
3. Rank variants in each population by their enrichment rate
4. Delete variants outside of the top 20 in each population
5. From the 'PS4_over_naive_enrichment.csv' file generated using the `print_enrich()` function, copy and paste the top 20 variants (by enrichment rate) with mutations in distal positions into the 'PS4_proximal_subs.csv' file. The population label for these variants will be 'XXXXX'
6. Ensure the path directory is set to the figures folder.
7. Execute fig3b.py directly in the python interpreter.

### Supplementary Figure 6b

Note that the generation of this figure requires some manual data processing.

1. Run the `count_Stops()` function from the analysis module for each population (as described in the appropriate subsection above).
2. Move all .csv files into the 'figures' folder
3. In one .csv file, delete the first two rows so that cell A1 contains the label, "Variant"
4. From the remaining files, copy all of the data below the "Variant", "log2(Count)", and "Library" headings to the first .csv file.
5. Ensure the path directory is set to the figures folder.
6. Rename this file "stop_box"
7. Execute figS6b.py directly in the python interpreter.

## Application of BRIDGE to library containing three distal loci

The following steps provide instructions for processing an A<sub>2</sub>aR library containing mutations at positions T88<sup>3.36</sup> and Q89<sup>3.37</sup> (locus 1), L249<sup>6.51</sup> and H250<sup>6.52</sup> (locus 2), and S281<sup>7.46</sup> (locus 3) (see Supplementary Figure 14).

### Merge read files from different flow cell lanes

If necessary, follow the instructions provided previously to merge forward (R1) and reverse (R2) read files from separate lanes into single R1 and R2 files.

### Trimming and quality-filtering raw PS1_reads

In order to properly trim and quality-filter raw reads, the `trim_qc()` function is modified as shown below.

```
def trim_qc(read1_file, read2_file, population, r1_start = 'CGTCCTGG', r1_end = 'CCATCTTCA', r1_len = 32, r2_start = 'GCAGTTGA', r2_end = 'GCAGAGG', r2_len = 34, r3_start = 'TCACAAC', r3_end = 'ATTGGTG', r3_len = 17, min_Qscore = 20):
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
                            locus2_read = rec2[search_start:search_end+len(r2_end)]
            if rec2.seq.count(r3_start) == 1:
                if rec2.seq.count(r3_end) == 1:
                    search_start = rec2.seq.find(r3_start)
                    search_end = rec2.seq.find(r3_end)
                    if len(rec2.seq[search_start:search_end+len(r3_end)]) == r3_len:
                        if min(rec2.letter_annotations['phred_quality'][search_start:search_end+len(r3_end)]) >= min_Qscore:
                            locus3_read = rec2[search_start:search_end+len(r3_end)]
            SeqIO.write(locus3_read+locus2_read, fout, 'fastq')
    return
```

### Append reads with matching flow cell (x,y) coordinates

The code utilized to match non-overlapping paired-end reads (i.e. `append_reads()`) does not have to be modified as it relies solely on flow cell coordinates encoded within the FASTQ sequence identifier and is sequence-agnostic. Accordingly, previous instructions can be used to append the .fastq files generated from the modified `trim_qc()` function.

### Translate mutated codons to amino acid sequence

To translate mutated codons into their corresponding amino acids, modify the `translate_reads()` function as shown below:

```
def translate_reads(read_file, population):
    record_iterator = SeqIO.parse(read_file, 'fastq')
    with open(population + '_residues.fastq', 'w+') as fout:
        for rec in record_iterator:
            translation = rec.seq[13:19].translate() + rec.seq[48:54].translate() + rec.seq[73:76].translate()
            line = str(translation)
            fout.write(line + '\n')
    return
```
