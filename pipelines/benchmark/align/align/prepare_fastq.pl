use strict;
use warnings;

# action:
#     prepare FASTQ for alignment from input read files of different types
#     enforce PLATFORM_MIN_INSERT_SIZE
#     trim fixed-length short-read platform reads to $FIXED_READ_LENGTH
#        typically, this is done to trim the low-quality 3' "extra base"
#     propagate required tags from unaligned bam through FASTQ to aligned BAM
#     tabulate observed insert size distributions, i.e., non-fixed read lengths, before any downstream filtering or projection
# expects:
#     source $MODULES_DIR/align/set_read_file_vars.sh
#     source $MODULES_DIR/library/set_library_vars.sh
#     $FIXED_READ_LENGTH when $IS_FIXED_READ_LENGTH==TRUE
# output: 
#     FASTQ on stdout with tags as QNAME extensions
#         is interleaved FASTQ for paired reads

# initialize reporting
our $action = "prepare_fastq";
my ($nInputReads, $nInputEvents, $nOutputEvents, $baseCounts) = (0) x 10;
my %insertSizeCounts;

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
resetCountFile();

# environment variables
fillEnvVar(\my $READ_PAIR_TYPE,      'READ_PAIR_TYPE');
fillEnvVar(\my $READ_LENGTH_TYPE,    'READ_LENGTH_TYPE');
fillEnvVar(\my $IS_FIXED_READ_LENGTH,'IS_FIXED_READ_LENGTH');
fillEnvVar(\my $FIXED_READ_LENGTH,   'FIXED_READ_LENGTH');
fillEnvVar(\my $READ_FILE_TYPE,      'READ_FILE_TYPE');
fillEnvVar(\my $READ_1_FILES,        'READ_1_FILES');
fillEnvVar(\my $READ_2_FILES,        'READ_2_FILES');
fillEnvVar(\my $PLATFORM_MIN_INSERT_SIZE,     'PLATFORM_MIN_INSERT_SIZE');
fillEnvVar(\my $UNFILTERED_INSERT_SIZES_FILE, 'UNFILTERED_INSERT_SIZES_FILE');
$READ_1_FILES =~ s/[\n|\r]/ /g;
$READ_2_FILES =~ s/[\n|\r]/ /g;
my $isPairedReads = ($READ_PAIR_TYPE eq "paired") ? 1 : 0;
my $sizePlotBinSize = $READ_LENGTH_TYPE eq "short" ? 10 : 250; # bp

# constants
use constant {
    _QNAME  => 0,
    _TAGS   => 1,
    _SEQ    => 2,
    _QUAL   => 3,
    _PASSED_SIZE_FILTER => 4,
};

# open and pass the file input handles
my $keepQnameTags;
if($READ_FILE_TYPE eq "unaligned.bam" or $READ_FILE_TYPE eq "bam"){
    if($isPairedReads){
        die "$action error: paired reads not yet supported for unaligned bam input\n";
    } else {
        $keepQnameTags = 1;
        foreach my $ubamFile(split(" ", $READ_1_FILES)){
            my $modTags    = "ML,MM"; # tags do not need to be present
            my $ontTags    = "ch,rn,fn,tl";
            my $pacBioTags = "ff,ec,dt,dd,sk"; # not ip,pw
            open my $inH, "-|", "samtools view $ubamFile | ".
                                "samtools fastq -T $modTags,$ontTags,$pacBioTags -  2>/dev/null" 
                or throwError("could not open $ubamFile: $!");
            runSingleReads($inH);
            close $inH;
        }
    }
} elsif($READ_FILE_TYPE eq "fastq.gz"){
    if($isPairedReads){
        open my $inH1, "-|", "zcat $READ_1_FILES" or throwError("could not open $READ_1_FILES: $!");
        open my $inH2, "-|", "zcat $READ_2_FILES" or throwError("could not open $READ_2_FILES: $!");
        runReadPairs($inH1, $inH2);
        close $inH1;
        close $inH2; 
    } else {
        open my $inH, "-|", "zcat $READ_1_FILES" or throwError("could not open $READ_1_FILES: $!");
        runSingleReads($inH);
        close $inH;
    }
} elsif($READ_FILE_TYPE eq "fastq"){
    if($isPairedReads){
        open my $inH1, "-|", "cat $READ_1_FILES" or throwError("could not open $READ_1_FILES: $!");
        open my $inH2, "-|", "cat $READ_2_FILES" or throwError("could not open $READ_2_FILES: $!");
        runReadPairs($inH1, $inH2);
        close $inH1;
        close $inH2; 
    } else {
        open my $inH, "-|", "cat $READ_1_FILES" or throwError("could not open $READ_1_FILES: $!");
        runSingleReads($inH);
        close $inH;
    }
} else {
    die "unrecognized READ_FILE_TYPE: $READ_FILE_TYPE\n"; 
}

# run single reads
# no need to multi-thread since the downstream alignment process is rate limiting
sub runSingleReads {
    my ($inH) = @_;
    while (my $read = getRead($inH)){
        $nInputReads++;
        $nInputEvents++; 
        my $len = length($$read[_SEQ]);
        $insertSizeCounts{int($len / $sizePlotBinSize) * $sizePlotBinSize}++;
        if($$read[_PASSED_SIZE_FILTER]){
            $nOutputEvents++;
            my $nameLine = $$read[_TAGS] ? "$$read[_QNAME] $$read[_TAGS]" : $$read[_QNAME];
            print join("\n", $nameLine, $$read[_SEQ], '+', $$read[_QUAL]), "\n";
            $baseCounts += $len;
        }
    }
}

# run paired reads
sub runReadPairs {
    my ($inH1, $inH2) = @_;
    while (my $read1 = getRead($inH1)){
        my $read2 = getRead($inH2);
        $nInputReads += 2;
        $nInputEvents++;
        if($$read1[_PASSED_SIZE_FILTER] and $$read2[_PASSED_SIZE_FILTER]){
            $nOutputEvents++;
            print join("\n", $$read1[_QNAME], $$read1[_SEQ], '+', $$read1[_QUAL]), "\n"; # print interleaved read pairs
            print join("\n", $$read1[_QNAME], $$read2[_SEQ], '+', $$read2[_QUAL]), "\n";
            $baseCounts += length($$read1[_SEQ]);
            $baseCounts += length($$read2[_SEQ]);
        }
    }
}

# parse a FASTQ set of 4 lines
sub getRead {
    my ($inH) = @_;
    my $name = <$inH>;
    $name or return; # EOF
    chomp $name;
    my ($qName, $tags);
    if($keepQnameTags){
        ($qName, $tags) = split(/\s/, $name, 2);
    } else {
        ($qName) = split(/\s/, $name); # discard any platform-specific non-name content on FASTQ name line
        $tags = "";
    }
    my $seq     = <$inH>;
    my $discard = <$inH>;
    my $qual    = <$inH>; 
    chomp $seq;
    chomp $qual;
    if($IS_FIXED_READ_LENGTH){
        $seq  = substr($seq,  0, $FIXED_READ_LENGTH);
        $qual = substr($qual, 0, $FIXED_READ_LENGTH);
    }
    return [
        $qName, 
        $tags, 
        $seq, 
        $qual, 
        length($seq) >= $PLATFORM_MIN_INSERT_SIZE
    ];
}

# print summary information
printCount(commify($nInputReads),   'nInputReads',   'input reads');
printCount(commify($nInputEvents),  'nInputEvents',  'input events, i.e., single reads or read pairs');
printCount(commify($nOutputEvents), 'nOutputEvents', "output events after size filtering >= $PLATFORM_MIN_INSERT_SIZE bp");
printCount(commify($baseCounts),    'baseCounts',    'read bases in output events');

# print insert size counts
unless($isPairedReads){
    open my $outH, '>', "$UNFILTERED_INSERT_SIZES_FILE" or throwError("could not open insert sizes file: $!");
    my @insertSizes = sort { $a <=> $b } keys %insertSizeCounts;
    print $outH join("\n", map { 
        join("\t", 
            "unfiltered",
            $_, 
            $insertSizeCounts{$_}
        )
    } @insertSizes), "\n";
    close $outH;
}
