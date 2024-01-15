Structural Variation Detection

basic example that detects insertions & deletions in long reads
The example loads up the ecoli genome, adds "markers" at 4bp pattern locations,
creates a fake long read with an insertion, and another with a deletion,
then analyzes the read's anchor patterns to determine any change with respect
to the reference.

The thought is that a technology could nick a genome at specific locations, using a pattern, I
use 'cgat' in this case, and replace with a marker, like a fluorescent label for example,
then detect the positions of those markers as they are observed by a sensor.  So in a
system like a microscope with an avalanche photo diode camera system and a ploymerase
attached to a slide, it would pull the single strand along as it incorporates, and the
marker could be observed over time, and it's position on the DNA strand determined, based
on the incorporation rate.  So multiple markers can be aligned (anchored) to exstablish the
rate, and used to determine insertions or deletions using both ends of the read.

