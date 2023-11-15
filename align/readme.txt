align_reads.py

Very basic aligner that generates a k-mer hash list for candidate genome positions and performs edit-dist against each position, returning the best position for each read.

Takes as input a fastq file or a simple text file, along with a referenge genome.  The genome is indexed when loaded, so no need to create an index file up front.  The reads are mapped against each candidate position in the genome and edit distance is used to determine the best position.  The full length of the read is used during edit distance, so this is not a great seed & extend, more of a seed & align algorithm.  It runs very fast on 10,000 reads, taking only a couple seconds, including generation of the k-mer hash list for the genome.

For reverse complement read alignments, a reverse complement genome is created, rather than reverse-comp the reads.  This preserves the error profile of the reads for better mapping.  Reversing the reads would result in higher initial errors when generating the k-mer for the read, wich would result in lower reverse complement alignments, thus a copy of the genome is generated, reverve complemented, and a second k-mer list generated.  This uses 2x the RAM, but ensures accurate alignments.

to run as a demo, simply run the align_reads.py with no arguments.  The default is to load up the included PhiX174 genome and the simulated 10,000 reads generated form PhiX with some random base substitutions.

