# tngs

The tngs project aims to provide a collection of tools supporting work with microbial genomics. This is the perl edition which is designed to provide a collection of useful command line tools in a single script with few dependencies. 

Current tools:

fileinfo  (basic qc info from fastq or fasta WGS-assembly)

downsample (downsamles fastq files)

filtercontigs (remove short and/or low-coverage contigs)

kmeroverlapp (counts how many k-mers present in a reference genome are present in an assembly)

kmerhisto (makes histogram of k-mers in fastq read files, estimates genome size, measures gc-bias)

