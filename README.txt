=== HapMaker Documentation ===
23 May 2013 Nozomu Okuda

-- Introduction --

HapMaker is a DNA sequence evolution simulator.  Its purpose is to generate DNA
sequences which differ from a given DNA sequence, especially with regards to
haplotypic variance.  Hopefully, HapMaker will be useful for development of
diploid (and maybe even polyploid) genome assemblers.

-- Dependencies --

In order to run HapMaker, a computer must have Perl and the BioPerl library
installed.  HapMaker was developed using Perl v5.16.0 and BioPerl v1.6.1 in a
Linux (Red Hat Enterprise Linux 6.2) environment.

-- Usage --

HapMaker's usage is as follows:

    perl HapMaker.pl [ref] [no-change] [het] [var] [loc]

The [ref] input parameter is the path to a fasta format DNA sequence file
which will be used as the reference sequence from which HapMaker creates the
variant haplotype.  Currently, HapMaker can handle only one contiguous DNA
sequence as a reference for producing the variant haplotype.  Thus, HapMaker
will consider only the first sequence given in [ref].

The [no-change] input parameter is the path to a no-change file.  In the
context of HapMaker, a no-change file marks which regions of the reference
sequence will not be considered for mutation in producing the variant
haplotype sequence.  The no-change file follows the following format:

    * Each region not to be changed is separated by a newline
    * The base number (the number is 1-based) where the region of no-change starts
        is first on the line
    * A hyphen (-) immediately follows the start base
    * The last base of the no-change region immediately follows the hyphen
    * The start base must be less than or equal to the last base
    * If every base is to be considered for mutation, the no-change file
        should have nothing in it except for a newline character

The [het] input parameter is a decimal number representing the relative amount
of difference between the reference sequence and variant sequence to be
produced.  It should be a decimal number between 0 (exclusive) and 1
(exclusive), representing the percent difference between the reference and
variant sequence.  Behavior is undefined when het is set to a percentage higher
than the percentage of bases available for mutation on the reference sequence.
HapMaker calculates percent difference between the two sequences by the number
of bases changed divided by the total number of bases in the reference
sequence.  HapMaker considers single nucletoide polymorphism events as having
changed one base, and indel events as having changed a number of bases equal to
the number of bases inserted/deleted in that event.  As such mutation events
are based on Perl's standard random number generator, there is no guarantee
that the heterozygosity between the reference and variant sequences will be
exactly as specified by [het], but HapMaker will print out the percent
difference it calculates, according to the description above, when it has
completed generating a variant sequence.

The [var] input parameter is a string used to name the variant sequence both
in the header of the output fasta file as well as to name the output fasta
file itself.  Thus, [var].fa will be an output file in the working directory
where HapMaker was run.

The [loc] input parameter is a string used to name the loc file.  The loc file
will be an output file in the working directory where HapMaker was run.  A loc
file explicitly states the mapping of a base in the variant sequence to a base
in the reference sequence.  As in the no-change file, the bases are counted on
a 1-based system.  The mapping is simply the line number in the file (which
represents the base number in the variant sequence) to whatever is printed on
that line (which corresponds to the reference sequence).  A number in the loc
file represents the base number of the reference sequence to which the base
number (as specified by the line number) of the variant sequence corresponds.
In cases of a single nucleotide polymorphism event, the line will also contain
the following notation:

    [reference sequence's base]->[variant sequence's base]

A hyphen in the loc file represents a base not present in the reference
sequence.  In other words, the base of the variant sequence at the location
specified by the line number on which a hyphen occurs was inserted during an
insertion event.  When numbers from successive lines are different by more
than one, a deletion event occurred there.
