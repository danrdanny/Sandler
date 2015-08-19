#!/usr/bin/perl

use strict;
use Getopt::Std;

# User-defined variables
my $minReadLength 	= 100;
my $minReadQuality 	= 60;
my $minBaseQuality 	= 20;
my $eventFrequency	= 30; # discard events seen in more than this percentage of reads within a 1kb window
my $minVCFScore		= 200;

## Command-line options
my %opts;
getopts('a:b:c:r:s:h', \%opts); # Values in %opts

### Usage Statement
if ($opts{'h'} || !$opts{'a'} || !$opts{'b'} || !$opts{'s'}) {
	print "
	Usage statement.

	Usage: 
		Typical:
		\$ Sandler.pl -a parentA.vcf -b parentB.vcf -s pooledBAM.bam

		With multiple bam files (quoted, with a comma, and no spaces between file names):
		\$ Sandler.pl -a parentA.tsv -b parentB.tsv -s \"pooledBAM1.bam,pooledBAM2.bam,etc\"

		With repeatmasker file:
		\$ Sandler.pl -a parentA.vcf -b parentB.vcf -s pooledBAM.bam -r repeats.txt

		Only looking at chrX and chr3R:
		\$ Sandler.pl -a parentA.vcf -b parentB.vcf -s pooledBAM.bam -r repeats.txt -c \"chrX,chr3R\"

	Required flags:

	-a	Parent A VCF or tsv file of SNPs. VCF file must be v4.0 or greater. If supplying
		a tsv file the format should be <chr>\\t<ID>\\t<base>, lines that begin with \#
		are ignored.

	-b	Same as -a, for parent B.

	-s	bam files to be analyzed. 

	Optional flags:

	-c	List of chromosomes to check, may or may not include ranges. If passing more than one
		chromosome with a comma you MUST quote the line. Examples are:

			chr21				Only check chr21
			\"chrX,chr18\"			Checks chrX and chr18
			\"chrX:1000000-5000000,chr18\"	Checks chrX positions 1-5mm and chr18

	-h	This helpful help.

	-r	Repeatmasker file - highly reccomended to cut down on the number of false positives.
		The file must be in the format (example from ce10.fa.out):

   			   SW  perc perc perc  query      position in query
			score  div. del. ins.  sequence    begin     end
			  508   0.0  0.0  0.0  chrI            1     432 
 			 1226  10.0  0.0  0.0  chrI          566     595

		Files can typically be found here (last accessed 2015-08):
			http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html



	\n";
	exit 0;
}

#------------------------------------#
#  Subroutine to check memory usage  #
#------------------------------------#
sub memUsage {
	my $gb;
	if (-e "/proc/$$/status") {
		open DATA, "< /proc/$$/status" or die "Unable to read /proc/$$/status: $!\n";
		local $/;
		<DATA> =~ m/^VmSize:\s+(.*?)$/m;
		my($kb) = $1 =~ /(\d+)/;
		$gb = sprintf("%0.1f", $kb / 1000 / 1000);
		$gb .= " Gigabytes";
	} else {
		$gb = "unable to check on this OS\n";
	}
	return $gb;
}

#-----------------------#
#  Parse bam file list  #
#-----------------------#
my(@alignmentFiles,$alignmentFileCount); 
foreach my $file (split /\,/, $opts{'s'}) {
	$file =~ s/\s+//g;
	next unless $file =~ /[a-z0-9]/i;

	die "Input file $file not found. Exiting.\n" if !-e $file;

	push(@alignmentFiles,$file);
	++$alignmentFileCount;
}
print "[".localtime(time)."] $alignmentFileCount BAM file(s) found.\n";

my %chrSizes;
my $bamHeader = `samtools view -H $alignmentFiles[0]`;
foreach my $chunk (split /\n/, $bamHeader) {
	my($chr) = $chunk =~ /SN\:(\w+)/;
	my($len) = $chunk =~ /LN\:(\d+)/;

	$chrSizes{$chr} = $len;
}

#---------------------------------------------------------------------------------#
#  Let's not waste anyone's time, check to make sure essential files are present  #
#---------------------------------------------------------------------------------#
die "Parent A vcf/tsv file $opts{'a'} not found. Exiting.\n" if !-e $opts{'a'};
die "Parent B vcf/tsv file $opts{'b'} found. Exiting.\n" if !-e $opts{'b'};
print "[".localtime(time)."] Parent A and B vcf/tsv files found.\n";

#-------------------------#
#  Parse chromosome list  #
#-------------------------#
my(%chromosomes,$chromosomeCount);
if ($opts{'c'}) {
	#print "[".localtime(time)."] -c flag given, splitting custom chromosome list.\n";
	foreach my $chr (split /\,/, $opts{'c'}) {
		$chr =~ s/\s+//g;
		my($lowRange,$highRange);
		my $range = 1;
		if ($chr =~ /\:/) {
			($chr,$lowRange,$highRange) = $chr =~ /(\w+)\:(\d+)\-(\d+)/;
			$range = "$lowRange-$highRange";
		}

		next unless $chr =~ /[a-z0-9]/i;
		$chromosomes{$chr} = $range;
		++$chromosomeCount;
	}
	print "[".localtime(time)."] $chromosomeCount individual chromosome(s) passed in. Summary:\n\n";
	print "\t\t\t\t\tChr\tRange (if given)\n";

	foreach my $chr (sort keys %chromosomes) {
		my $range = $chromosomes{$chr};
		   $range = undef if $range == 1;
		print "\t\t\t\t\t$chr\t$range\n";
	}
	print "\n";
} else {
	print "[".localtime(time)."] No specific chromosomes given, will search all in parental files.\n";
}

#--------------------------#
#  Open repeatmasker file  #
#--------------------------#
my %repeats;
if ($opts{'r'}) {
	print "[".localtime(time)."] -r flag given, opening repeatmasker file $opts{'r'}\n";

	my $countMaskedBases = 0;
	open INF,"$opts{'r'}" or die "can't open file $opts{'r'}: $!";
	while (<INF>) {
		$_ =~ s/^\s+//;
		my @F = split /\s+/, $_;

		my($chr,$lowBase,$highBase) = ($F[5],$F[6],$F[7]);
		$chr =~ s/chr//;

		next unless $chr =~ /[0-9]/;

		# if specific chromsomes are targeted, we can add only those to the hash
		# this saves some serious memory
		next if $chromosomeCount > 0 && !$chromosomes{$chr};
	
		foreach ($lowBase..$highBase) {
			$repeats{$chr}{$_} = 1;
			++$countMaskedBases;
		}
	}
	print "[".localtime(time)."] repeatmasker file read, $countMaskedBases bases masked.\n";
	print "[".localtime(time)."] WARN: 0 bases masked!\n" if $countMaskedBases == 0;
} else {
	print "[".localtime(time)."] No repeatmasker file given - careful, this could get messy!\n";
}

print "[".localtime(time)."] Memory usage: ".memUsage()."\n";

#----------------------------------#
#  Open parental vcf or tsv files  #
#----------------------------------#

my(%parentASNPs,%parentBSNPs,%snpPos,$parentCount);
foreach my $file ($opts{'a'}, $opts{'b'}) {
	$parentCount++;
	my $parent = "ParentA";
	   $parent = "ParentB" if $parentCount == 2;

	my($skipScore,$skipHet,$skipIndel,$skipRepeat,$skipAlt) = (0,0,0,0,0);

	open INF,"$file" or die "Can't open file $file: $!";
	while (<INF>) {
		my(@F) = split /\t/, $_;

		my($chr,$id,$base) = ($F[1],$F[2],$F[3]);

		if ($F[5] > 0) { # probably a vcf file
			$chr  = $F[0]; #CHROM
			$id   = $F[1]; #POS
			$base = $F[4]; #ALT


			next if $chromosomeCount > 0 && !$chromosomes{$chr};

			++$skipScore if $F[5] < $minVCFScore;
			next if $F[5] < $minVCFScore;

			++$skipHet if $_ =~ /0\/1/; # only keep homozygous reads
			next if $_ =~ /0\/1/;

			++$skipIndel if $_ =~ /INDEL/;
			next if $_ =~ /INDEL/;

			++$skipAlt if $base =~ /\,/;
			next if $base =~ /\,/;
		}

		$skipRepeat++ if $repeats{$chr}{$id};
		next if $repeats{$chr}{$id};

		next if $chromosomeCount > 0 && !$chromosomes{$chr};

		next unless $chr;
		next unless $id;
		next unless $base;

		$snpPos{$chr}{$id} = $parent;
		$parentASNPs{$chr}{$id} = $base if $parentCount == 1;
		$parentBSNPs{$chr}{$id} = $base if $parentCount == 2;
	}
	close INF;

	printf "[".localtime(time)."] \n";
	printf "[".localtime(time)."] $parent $file opened, summary:\n";
	printf "[".localtime(time)."] %30s %8d\n", "Skip d/t repeat", $skipRepeat;
	printf "[".localtime(time)."] %30s %8d\n", "Skip d/t score < $minVCFScore", $skipScore;
	printf "[".localtime(time)."] %30s %8d\n", "Skip d/t het SNP", $skipHet;
	printf "[".localtime(time)."] %30s %8d\n", "Skip d/t two alt alleles", $skipAlt;
	#printf "[".localtime(time)."] %30s %8d\n", "Skip d/t depth <= $minParentalDepth", $skipD    epth;
}

#-----------------------------------#
#  Get rid of shared parental SNPs  #
#-----------------------------------#
my(%count,$skipSharedSNP);
foreach my $chr (sort keys %snpPos) {
	foreach my $id (sort keys %{$snpPos{$chr}}) {
		if ($parentASNPs{$chr}{$id} eq $parentBSNPs{$chr}{$id}) {
			$snpPos{$chr}{$id} = undef;
			$parentASNPs{$chr}{$id} = undef;
			$parentBSNPs{$chr}{$id} = undef;
			++$skipSharedSNP;
		} else {
			$count{$chr}++;
		}
	}
}

printf "[".localtime(time)."] Compared ParentA to ParentB. Skipped $skipSharedSNP d/t shared SNP.\n";
printf "[".localtime(time)."] \n";
printf "[".localtime(time)."] SNPs per chromosome:\n";
foreach my $chr (sort keys %count) {
	print "\t\t\t\t$chr\t$count{$chr}\n";
}
printf "[".localtime(time)."] \n";

#-----------------------------------------------------------------------------------------#
#  Calculate gaps between SNPs - letting us know where it is feasible to identify events  #
#  Also build feeder file for coverage heatmap						  #
#-----------------------------------------------------------------------------------------#
my(%gapSummary,%graphGaps);
foreach my $chr (sort keys %snpPos) {
	my($lastParent,$lastID);
	foreach my $id (sort {$a<=>$b} keys %{$snpPos{$chr}}) {
		my $currParent = $snpPos{$chr}{$id};

		if ($lastParent && $lastParent ne $currParent) {
			my $gap = $id - $lastID;

			$gapSummary{$chr}{'1-250'}++      if $gap <= 250;
			$gapSummary{$chr}{'251-500'}++    if $gap > 250 && $gap <= 500;
			$gapSummary{$chr}{'501-1000'}++   if $gap > 500 && $gap <= 1000;
			$gapSummary{$chr}{'1001-2000'}++  if $gap > 1000 && $gap <= 2000;
			$gapSummary{$chr}{'2001-5000'}++  if $gap > 2000 && $gap <= 5000;
			$gapSummary{$chr}{'5001-10000'}++ if $gap > 5000 && $gap <= 10000;
			$gapSummary{$chr}{'>10000'}++ 	  if $gap > 10000;

			push(@{$graphGaps{$chr}},$id) if $gap <= 1000;
		}

		$lastParent = $currParent;
		$lastID = $id;
	}
}

foreach my $chr (sort keys %gapSummary) {
	foreach my $range ("1-250","251-500","501-1000","1001-2000","2001-5000","5001-10000",">10000") {
		printf "%5s\t%15s\t%8s\n", $chr, $range, $gapSummary{$chr}{$range};
	}
	print "\n";
}

# Build heatmap file
print "[".localtime(time)."] Making heatmap feeder file for gaps < 1kb.\n";
my $output = "Chr\tPosition\tPercent\n";
foreach my $chr (sort keys %gapSummary) {
	print "[".localtime(time)."]   Chr: $chr\n";
	my $maxPos 	= $chrSizes{$chr};
	my $currMin 	= 0;
	my $step 	= 5000; # in kb
	my $currMax 	= 5000;
	while ($currMax < $maxPos) {
		my $count = 0;
		foreach my $id (@{$graphGaps{$chr}}) {
			next unless $id >= $currMin && $id <= $currMax;
			++$count;
		}
		my $percent = sprintf("%0.3f", ( ($count / $step) * 100));
		$output .= "$chr\t$currMax\t$percent\n";
			
		$currMin = $currMax + 1;
		$currMax += $step;
	}
}
open OUTF,">./coverageHeatmapData.tsv";
print OUTF $output;
close OUTF;

print "[".localtime(time)."] Memory usage: ".memUsage()."\n";

#-------------------------------------------------------#
#  Open each bamfile and look for evidence of exchange  #
#-------------------------------------------------------#
foreach my $chr (sort keys %gapSummary) {
  foreach my $bamFile (@alignmentFiles) {
    my $count = 0;
    my $readID;
    my($firstStartID,$firstEndID,$firstSequence,$secondStartID,$secondEndID,$secondSequence);

    print "Samtools view\n";
    `samtools view -b $bamFile $chr > Sandler.temp.bam`; # we select only the chr we're interested in
    print "Samtools sort\n";
    `samtools sort -n Sandler.temp.bam > Sandler.temp.sorted.bam`;
    print "Samtools view\n";
    `samtools view Sandler.temp.sorted.bam > Sandler.temp.sam`;
    #`rm -f Sandler.temp.sorted.bam Sandler.temp.bam`;
		
    open INF,"./Sandler.temp.sam" or die "can't open file Sandler.temp.sam: $!";
    while (<INF>) {
	my @F = split /\t/, $_;

	if ($count == 0) {
		$firstSequence = undef;
		$secondSequence = undef;
		$firstEndID = undef;
		$firstStartID = undef;
		$secondEndID = undef;
		$secondStartID = undef;
	} # end of if count == 0
	
	if ($count == 1) {
		my $skip = 1;
		my($lowRange,$highRange) = $chromosomes{$chr} =~ /(\d+)-(\d+)/;
		if ($lowRange > 0 && $highRange > 0) {
			$skip = 0 if $firstStartID >= $lowRange && $firstStartID <= $highRange;
		} else {
			$skip = 0;
		}

		$count = 0;
		if ($readID eq $F[0] && $skip == 0) { # the read IDs need to be the same
			my $ogSequence = $F[9];
			my $seq;
			my $splitInfo = $F[5];

			foreach (split /([0-9]+[A-Z])/, $splitInfo) {
				my($len,$action) = $_ =~ /([0-9]+)([A-Z])/;
				next unless $len =~ /[0-9]/;

				if ($action =~ /M/) {
					$ogSequence =~ /^(.{$len})/;
					$seq .= $1;
					$ogSequence =~ s/^.{$len}(.+)/$1/;
				} elsif ($action =~ /D/) {
					foreach (1..$len) {
						$seq .= ".";
					}
				} elsif ($action =~ /I/) {
					$ogSequence =~ s/^.{$len}(.+)/$1/;
				} elsif ($action =~ /S/) {
					$ogSequence =~ s/^.{$len}(.+)/$1/;
				} else {
					print "Unknown action: $action: $len\t$readID\n";
				}
			}

			next if !$seq;
			$secondSequence = $seq;
			my $length = length($secondSequence);
			$secondStartID = $F[3];
			$secondEndID = $F[3] + $length - 1;

			#print "$secondStartID\t$secondEndID\n$secondSequence\n\n";

			my($parentACount,$parentBCount);
			foreach ($firstStartID..$firstEndID) {
				next unless $snpPos{$chr}{$_};
				#print "Snp at $_: $snpPos{$_}, $parents{$_}\n";
				$parentACount++ if $snpPos{$chr}{$_} =~ /parentA/;
				$parentBCount++ if $snpPos{$chr}{$_} =~ /parentB/;
			}
			foreach ($secondStartID..$secondEndID) {
				next unless $snpPos{$chr}{$_};
				#print "Snp at $_: $snpPos{$_}, $parents{$_}\n";
				$parentACount++ if $snpPos{$chr}{$_} =~ /parentA/;
				$parentBCount++ if $snpPos{$chr}{$_} =~ /parentB/;
			}

			#print "$firstStartID, $firstEndID\n$firstSequence\n";

			my($parentAMatch,$parentBMatch,$output);
			my $tempLoc = $firstStartID;
			foreach (split //, $firstSequence) {
				my $match;
				$match = "parentA" if $parentASNPs{$chr}{$tempLoc} eq $_ && $snpPos{$chr}{$tempLoc} eq "parentA";
				$match = "parentB" if $parentBSNPs{$chr}{$tempLoc} eq $_ && $snpPos{$chr}{$tempLoc} eq "parentB";
	
				$parentAMatch++ if $match eq "parentA";
				$parentBMatch++ if $match eq "parentB";
				$output .=  "$match\t$tempLoc\t" if $match;
				$tempLoc++;
			}

			#print "$secondStartID\t$secondEndID\n$secondSequence\n\n";

			my $tempLoc = $secondStartID;
			foreach (split //, $secondSequence) {
				#$strandSNPs{$tempLoc} = $_;
				#print "$tempLoc: $_\n";

				my $match;
				$match = "parentA" if $parentASNPs{$chr}{$tempLoc} eq $_ && $snpPos{$chr}{$tempLoc} eq "parentA";
				$match = "parentB" if $parentBSNPs{$chr}{$tempLoc} eq $_ && $snpPos{$chr}{$tempLoc} eq "parentB";
	
				$parentAMatch++ if $match eq "parentA";
				$parentBMatch++ if $match eq "parentB";
				$output .=  "$match\t$tempLoc\t" if $match;
				$tempLoc++;
				#sleep 1;
			}

			#print "csMatch: $csMatch\n";
			#print "wMatch: $wMatch\n";

			if ($parentACount > 0 && $parentBCount > 0 && $parentAMatch > 0 && $parentBMatch > 0) {
				my($max,$min);
				foreach ($firstStartID,$firstEndID,$secondStartID,$secondEndID) {
					$max = $_ if $_ > $max || !$min;
				}
				foreach ($firstStartID,$firstEndID,$secondStartID,$secondEndID) {
					$min = $_ if $_ < $min || !$min;
				}

				my $fragment = $max - $min;
				#next if $fragment > 2000; # not sure why I did that

				#print "$firstStartID, $firstEndID\n$firstSequence\n";
				#print "$secondStartID\t$secondEndID\n$secondSequence\n\n";
				print "$readID\t$min\t$max\t$fragment\t\t";
				print "$firstStartID\t$firstEndID\t$secondStartID\t$secondEndID\t$parentACount\t$parentAMatch\t$parentBCount\t$parentBMatch\t\t$output\n";
			}
		}

		$count = 0;
	} else { # hits this when count == 0
		$readID = $F[0];
		my $ogSequence = $F[9];
		my $seq;
		my $splitInfo = $F[5];

		foreach (split /([0-9]+[A-Z])/, $splitInfo) {
			my($len,$action) = $_ =~ /([0-9]+)([A-Z])/;
			next unless $len =~ /[0-9]/;

				if ($action =~ /M/) {
				$ogSequence =~ /^(.{$len})/;
				$seq .= $1;
				$ogSequence =~ s/^.{$len}(.+)/$1/;
			} elsif ($action =~ /D/) {
				foreach (1..$len) {
					$seq .= ".";
				}
			} elsif ($action =~ /I/) {
				$ogSequence =~ s/^.{$len}(.+)/$1/;
			} elsif ($action =~ /S/) {
				$ogSequence =~ s/^.{$len}(.+)/$1/;
			} else {
				print "Unknown action: $action: $len\t$readID\n";
			}

			#print "$action: $len\n";
			#print "new: $seq\n";
			#print "old: $ogSequence\n";
			#print "\n";
		}

		if ($seq) {
			$firstSequence = $seq;
			my $length = length($firstSequence);
			$firstStartID = $F[3];
			$firstEndID = $F[3] + $length - 1;

			#print "$firstStartID\t$firstEndID\n$firstSequence\n";
			#sleep 3;

			$count++;
		}
	} # end of if count == 1
      } # end of open INF
    #`rm -f ./Sandler.temp.sam`;
  } # end of bamFile loop
} # end of chromosome loop
