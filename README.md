# TargetSeq-QC
This is a perl script which can generate some statistics on the Targeted re-sequencing bam file.

## $0 -bed Amplicon.bedtools.bed -bam Sample_RMS222_A_H3FHYAFXX.trim.bam -out Sample_RMS222_A_H3FHYAFXX.QC.txt -int 10 -int 100 -int 500 -int 1000 -sum 1
#Options
###-h, -help, --help Print this message.
**-bed**    Bed file containing the locations where statistics to be generated. (Required)
**-bam**    Bam file on which the you whould like to generate the statistics, should be indexed. (Required)
**-sum**    1 if you want aggregate entry for the whole bed file. (Default does not aggregate.)
Overlapping regions in bed file will be merged to avoid overcounting.
**-int**    %of Based covered with min intX. Can be specified multiple times. (Required)
**-out**    Output file name. (Required)
