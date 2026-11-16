use strict;
use warnings;

# change the fastp merged tag to BAM format tag fm:Z: for passing through minimap2 as fastq
# count reads and bases to assess what fastp filtering and merging did to prepare_fastq.pl output
# final output tag format in the FASTQ name line is:
#   QNAME [any previous tags] fm:Z:nRead1,nRead2
# where:
#     fm indicates "fastp merge"
#     nRead1, the number of read1 bases present in the final merged read
#     nRead2, the number of read2 bases present in the final merged read
# the tag is not added to unmerged reads and will be missing from the eventual BAM file

# initialize reporting
our $action = "adjust_merge_tags";
my ($nInputReads, $nMergedReads, $nUnmergedReads, $nUnmergedEvents, $nPairedEvents) = (0) x 10;
my @baseCounts = (0, 0, 0);

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
resetCountFile();

# constants
use constant {
    READ1   => 1,
    READ2   => 2,
    MERGED  => 3,
};

# run the interleaved pairs
my $lineN = 0;
my ($prevQName, $readN) = ("");
while(my $line = <STDIN>){
    if($lineN % 4 == 0){ # the QNAME line
        $nInputReads++;
        chomp $line;
        my @f = split(" ", $line); 

        # QNAME [tags] merged_150_15 means that 150bp are from read1, and 15bp are from read2
        if(@f > 1 and $f[-1] =~ m/merged_(\d+)_(\d+)/){ 
            $line = join(" ", @f[0..$#f-1], "fm:Z:$1,$2")."\n";
            $nMergedReads++;

        # a single or unmerged read
        } else {
            $line .= "\n";
            $nUnmergedReads++;
            $prevQName eq $f[0] or $nUnmergedEvents++;
        }

        $readN = $prevQName eq $f[0] ? READ2 : READ1;
        $readN == READ2 and $nPairedEvents++;
        $prevQName = $f[0];
    } elsif($lineN % 4 == 1){ # the sequence line
        $baseCounts[$readN] += length($line) - 1;
    }
    print $line;
    $lineN++;
}

# print summary information
my $nInputEvents    = $nMergedReads + $nUnmergedEvents;
my $nOrphanedReads  = $nUnmergedReads - 2 * $nPairedEvents;
my $percentMerged   = $nMergedReads   / $nInputEvents * 100;
my $percentOrphaned = $nOrphanedReads / $nInputEvents * 100;
printCount(commify($nInputReads),        'nInputReads',       'input reads');
printCount(commify($nInputEvents),       'nInputEvents',      'input events');
printCount(commify($nMergedReads),       'nMergedReads',      'input merged reads == events with merged read pairs');
printCount(commify($nUnmergedReads),     'nUnmergedReads',    'input unmerged reads (single, paired, or orphaned reads)');
printCount(commify($nUnmergedEvents),    'nUnmergedEvents',   'input events with unmerged reads (single reads, unmerged read pairs, or orphaned reads)');
printCount(commify($nPairedEvents),      'nPairedEvents',     'input events with paired reads == read2 count');
printCount(commify($nOrphanedReads),     'nOrphanedReads',    'input read1 without a read2 partner = events that lost (or never had) a read partner');
printCount(roundCount2($percentMerged),  'percentMerged',     'percent of reads pairs merged');
printCount(roundCount2($percentOrphaned),'percentOrphaned',   'percent of reads pairs orphaned');
printCount(commify($baseCounts[READ1]),  'baseCounts[READ1]', 'read1 bases in input events (single, first paired, merged, or orphaned)');
printCount(commify($baseCounts[READ2]),  'baseCounts[READ2]', 'read2 bases in input events (second paired)');


# from https://github.com/OpenGene/fastp?tab=readme-ov-file#merge-paired-end-reads

# In the output file, a tag like merged_xxx_yyy will be added to each read name 
# to indicate that how many base pairs are from read1 and from read2, respectively. 
# For example,  
#     @NB551106:9:H5Y5GBGX2:1:22306:18653:13119 1:N:0:GATCAG merged_150_15 
# means that 
#     150bp are from read1, and 15bp are from read2
# fastp prefers the bases in read1 since they usually have higher quality than read2.
