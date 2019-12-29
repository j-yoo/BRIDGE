import itertools
from decimal import Decimal
from math import log10

noT = 'AILVMFWYNCQSDERHKGP*'
noQ = 'AILVMFWYNCSTDERHKGP*'
noW = 'AILVMFYNCQSTDERHKGP*'
noL = 'AIVMFWYNCQSTDERHKGP*'
noH = 'AILVMFWYNCQSTDERKGP*'
naive_XQWLH_Counts = {}
naive_TXWLH_Counts = {}
naive_TQXLH_Counts = {}
naive_TQWXH_Counts = {}
naive_TQWLX_Counts = {}
PS4_XQWLH_Counts = {}
PS4_TXWLH_Counts = {}
PS4_TQXLH_Counts = {}
PS4_TQWXH_Counts = {}
PS4_TQWLX_Counts = {}
naive_XXWLH_Counts = {}
naive_TQXXH_Counts = {}
naive_TQWXX_Counts = {}
naive_TQXLX_Counts = {}
PS4_XXWLH_Counts = {}
PS4_TQXXH_Counts = {}
PS4_TQWXX_Counts = {}
PS4_TQXLX_Counts = {}
naiveTripCounts = {}
PS4TripCounts = {}
for var1, var2, var3, var4, var5 in map(''.join, itertools.product(noT, noQ, noW, noL, noH, repeat=1)):
    naive_XQWLH_Counts[var1 + "QWLH"] = 0
    naive_TXWLH_Counts["T" + var2 + "WLH"] = 0
    naive_TQXLH_Counts["TQ" + var3 + "LH"] = 0
    naive_TQWXH_Counts["TQW" + var4 + "H"] = 0
    naive_TQWLX_Counts["TQWL" + var5] = 0
    PS4_XQWLH_Counts[var1 + "QWLH"] = 0
    PS4_TXWLH_Counts["T" + var2 + "WLH"] = 0
    PS4_TQXLH_Counts["TQ" + var3 + "LH"] = 0
    PS4_TQWXH_Counts["TQW" + var4 + "H"] = 0
    PS4_TQWLX_Counts["TQWL" + var5] = 0
    naive_XXWLH_Counts[var1 + var2 + "WLH"] = 0
    naive_TQXXH_Counts["TQ" + var3 + var4 + "H"] = 0
    naive_TQWXX_Counts["TQW" + var4 + var5] = 0
    naive_TQXLX_Counts["TQ" + var3 + "L" + var5] = 0
    PS4_XXWLH_Counts[var1 + var2 + "WLH"] = 0
    PS4_TQXXH_Counts["TQ" + var3 + var4 + "H"] = 0
    PS4_TQWXX_Counts["TQW" + var4 + var5] = 0
    PS4_TQXLX_Counts["TQ" + var3 + "L" + var5] = 0
    naiveTripCounts["TQ" + var1 + var2 + var3] = 0
    PS4TripCounts["TQ" + var1 + var2 + var3] = 0

with open("naive_residues.fastq") as fin:
    for line in fin:
        read = str(line).strip()
        if read in naive_XQWLH_Counts:
            naive_XQWLH_Counts[read] += 1
        elif read in naive_TXWLH_Counts:
            naive_TXWLH_Counts[read] += 1
        elif read in naive_TQXLH_Counts:
            naive_TQXLH_Counts[read] += 1
        elif read in naive_TQWXH_Counts:
            naive_TQWXH_Counts[read] += 1
        elif read in naive_TQWLX_Counts:
            naive_TQWLX_Counts[read] += 1
        elif read in naive_XXWLH_Counts:
            naive_XXWLH_Counts[read] += 1
        elif read in naive_TQXXH_Counts:
            naive_TQXXH_Counts[read] += 1
        elif read in naive_TQWXX_Counts:
            naive_TQWXX_Counts[read] += 1
        elif read in naive_TQXLX_Counts:
            naive_TQXLX_Counts[read] += 1
        elif read in naiveTripCounts:
            naiveTripCounts[read] += 1

with open("PS4_residues.fastq") as fin:
    for line in fin:
        read = str(line).strip()
        if read in PS4_XQWLH_Counts:
            PS4_XQWLH_Counts[read] += 1
        elif read in PS4_TXWLH_Counts:
            PS4_TXWLH_Counts[read] += 1
        elif read in PS4_TQXLH_Counts:
            PS4_TQXLH_Counts[read] += 1
        elif read in PS4_TQWXH_Counts:
            PS4_TQWXH_Counts[read] += 1
        elif read in PS4_TQWLX_Counts:
            PS4_TQWLX_Counts[read] += 1
        elif read in PS4_XXWLH_Counts:
            PS4_XXWLH_Counts[read] += 1
        elif read in PS4_TQXXH_Counts:
            PS4_TQXXH_Counts[read] += 1
        elif read in PS4_TQWXX_Counts:
            PS4_TQWXX_Counts[read] += 1
        elif read in PS4_TQXLX_Counts:
            PS4_TQXLX_Counts[read] += 1
        elif read in PS4TripCounts:
            PS4TripCounts[read] += 1

XQWLH_Enrichments = {}
TXWLH_Enrichments = {}
TQXLH_Enrichments = {}
TQWXH_Enrichments = {}
TQWLX_Enrichments = {}
XXWLH_Enrichments = {}
TQXXH_Enrichments = {}
TQWXX_Enrichments = {}
TQXLX_Enrichments = {}
PS4_TripEnrichments = {}
for key in PS4_XQWLH_Counts.keys():
    if Decimal(naive_XQWLH_Counts[key]) >= 15:
        XQWLH_Enrichments[key] = Decimal(PS4_XQWLH_Counts[key])/Decimal(naive_XQWLH_Counts[key])
for key in PS4_TXWLH_Counts.keys():
    if Decimal(naive_TXWLH_Counts[key]) >= 15:
        TXWLH_Enrichments[key] = Decimal(PS4_TXWLH_Counts[key])/Decimal(naive_TXWLH_Counts[key])
for key in PS4_TQXLH_Counts.keys():
    if Decimal(naive_TQXLH_Counts[key]) >= 15:
        TQXLH_Enrichments[key] = Decimal(PS4_TQXLH_Counts[key])/Decimal(naive_TQXLH_Counts[key])
for key in PS4_TQWXH_Counts.keys():
    if Decimal(naive_TQWXH_Counts[key]) >= 15:
        TQWXH_Enrichments[key] = Decimal(PS4_TQWXH_Counts[key])/Decimal(naive_TQWXH_Counts[key])
for key in PS4_TQWLX_Counts.keys():
    if Decimal(naive_TQWLX_Counts[key]) >= 15:
        TQWLX_Enrichments[key] = Decimal(PS4_TQWLX_Counts[key])/Decimal(naive_TQWLX_Counts[key])
for key in PS4_XXWLH_Counts.keys():
    if Decimal(naive_XXWLH_Counts[key]) >= 15:
        XXWLH_Enrichments[key] = Decimal(PS4_XXWLH_Counts[key])/Decimal(naive_XXWLH_Counts[key])
for key in PS4_TQXXH_Counts.keys():
    if Decimal(naive_TQXXH_Counts[key]) >= 15:
        TQXXH_Enrichments[key] = Decimal(PS4_TQXXH_Counts[key])/Decimal(naive_TQXXH_Counts[key])
for key in PS4_TQWXX_Counts.keys():
    if Decimal(naive_TQWXX_Counts[key]) >= 15:
        TQWXX_Enrichments[key] = Decimal(PS4_TQWXX_Counts[key])/Decimal(naive_TQWXX_Counts[key])
for key in PS4_TQXLX_Counts.keys():
    if Decimal(naive_TQXLX_Counts[key]) >= 15:
        TQXLX_Enrichments[key] = Decimal(PS4_TQXLX_Counts[key])/Decimal(naive_TQXLX_Counts[key])
for key in PS4TripCounts.keys():
    if Decimal(naiveTripCounts[key]) >= 15:
        PS4_TripEnrichments[key] = Decimal(PS4TripCounts[key])/Decimal(naiveTripCounts[key])

with open("PS4_proximal_subs.csv","w+") as fout:
    fout.write("variant" + "," + "Enrichment Rate" + "," + "Population" + "\n")
    for key in XQWLH_Enrichments.keys():
        list = key + "," + str(XQWLH_Enrichments[key]) + "," + "XQWLH/TXWLH/TQXLH/TQWXH/TQWLX" + "\n"
        fout.write(list)
    for key in TXWLH_Enrichments.keys():
        list = key + "," + str(TXWLH_Enrichments[key]) + "," + "XQWLH/TXWLH/TQXLH/TQWXH/TQWLX" + "\n"
        fout.write(list)
    for key in TQXLH_Enrichments.keys():
        list = key + "," + str(TQXLH_Enrichments[key]) + "," + "XQWLH/TXWLH/TQXLH/TQWXH/TQWLX" + "\n"
        fout.write(list)
    for key in TQWXH_Enrichments.keys():
        list = key + "," + str(TQWXH_Enrichments[key]) + "," + "XQWLH/TXWLH/TQXLH/TQWXH/TQWLX" + "\n"
        fout.write(list)
    for key in TQWLX_Enrichments.keys():
        list = key + "," + str(TQWLX_Enrichments[key]) + "," + "XQWLH/TXWLH/TQXLH/TQWXH/TQWLX" + "\n"
        fout.write(list)
    for key in XXWLH_Enrichments.keys():
        list = key + "," + str(XXWLH_Enrichments[key]) + "," + "XXWLH/TQXXH/TQWXX/TQXLX" + "\n"
        fout.write(list)
    for key in TQXXH_Enrichments.keys():
        list = key + "," + str(TQXXH_Enrichments[key]) + "," + "XXWLH/TQXXH/TQWXX/TQXLX" + "\n"
        fout.write(list)
    for key in TQWXX_Enrichments.keys():
        list = key + "," + str(TQWXX_Enrichments[key]) + "," + "XXWLH/TQXXH/TQWXX/TQXLX" + "\n"
        fout.write(list)
    for key in TQXLX_Enrichments.keys():
        list = key + "," + str(TQXLX_Enrichments[key]) + "," + "XXWLH/TQXXH/TQWXX/TQXLX" + "\n"
        fout.write(list)
    for key in PS4_TripEnrichments.keys():
        list = key + "," + str(PS4_TripEnrichments[key]) + "," + "TQXXX" + "\n"
        fout.write(list)
