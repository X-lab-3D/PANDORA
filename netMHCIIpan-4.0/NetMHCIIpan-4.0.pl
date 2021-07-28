#! /usr/bin/perl -w

# Author: Edita Karosiene, 07-01-2013
#         Massimo Andreatta, 20-12,2014
# 	      Kamilla Kjaergaard Jensen, 2017
#	      Birkir Reynisson, 2019
# update 20-07-16: removing affinity-based thresholds to define weak and strong binders.
#                   the new defaults for strong and weak binders are ranks 2% and 10% 
# Updated 10-10-19: Option added for BA and EL output neuron read-outs. Also prediction with context encoding

use strict;

use Getopt::Long;
use Env;
use Sys::Hostname;

#set environment variable for $NETMHCIIpan
my $NETMHCIIpan = $ENV{'NETMHCIIpan'};
# determine platform
my $UNIX = `uname -s`;
my $AR = `uname -m`;
chomp $UNIX;
chomp $AR;
my $platform_name = "$UNIX" . "_" . "$AR";
# get environment variable for $PLATFORM
my $PLATFORM = $ENV{'NetMHCIIpanPLAT'}; 

#Specify R version (works with R-2.14)
#my $R = qx(which R-2.14);
my $R = "/tools/bin/R-2.14";

# specify version of the method
my $version = "4.0";

# variables used to save options
my ($BA, $context, $rdir, $inptype, $length, $exclude_offset, $tdir, $choose_chain, $chain_a, $chain_b, $alleles, $thr_rank_S, $thr_rank_W, $filter, $rank_f, $file, $dirty, $hlaseq, $hlaseq_a, $sort, $hist, $histhr, $xls, $xlsfile, $unique, $w, $list, $verbose, $help, $pboth);

# associating options values with variables
&GetOptions (
'BA'=>\$BA,
'context'=>\$context,
'rdir:s' => \$rdir,
'tdir:s' => \$tdir,
'ex_offset' => \$exclude_offset,
'inptype:i' => \$inptype,
'length:s' => \$length,
'rankS=s' => \$thr_rank_S,
'rankW=s' => \$thr_rank_W,
'choose' => \$choose_chain,
'cha=s' => \$chain_a,
'chb=s' => \$chain_b,
'filter' => \$filter,
'rankF=s' => \$rank_f,
'xlsfile:s' => \$xlsfile,
'a:s' => \$alleles,
'f:s'=>   \$file,
'hlaseq:s' => \$hlaseq,
'hlaseqA:s' => \$hlaseq_a,
's' => \$sort,
'hist' => \$hist,
'histhr=s' => \$histhr,
'u' => \$unique,
'list' => \$list,
'xls' => \$xls,
'w' => \$w,
'dirty' => \$dirty,
'v' => \$verbose,
'h' => \$help
);

# Program directory
$rdir = $NETMHCIIpan; 

# Define session Id
my $sId = $$;

### Scripts, files and data
#my $nnalign_pan_player = "$PLATFORM/bin/nnalign_gaps_pan_play-2.0b_PA";
my $nnalign_pan_player = "$PLATFORM/bin/nnalign_gaps_pan_play_MA_wgt_MN_v2_two_outputs_context_allelelist_pboth";
my $mhcfsa2psseq = "$PLATFORM/bin/mhcfsa2psseq";
my $make_P1_histograms = "$PLATFORM/bin/make_P1_hist_frames.R";
my $convert_HLA_list = "$rdir/data/convert_pseudo.dat";
my $pseudoseq_file = "$rdir/data/pseudosequence.2016.all.X.dat";
my $blosum_freq = "$rdir/data/blosum62.freq_rownorm";
my $allelelist = "$rdir/data/allelelist.txt";
#my $synlists = "$rdir/data/synlists";
my $synlists = "$rdir/data/synlists.bin";

#my $wwwdir = "/usr/opt/www/pub/CBS/services/NetMHCIIpan-3.2/tmp";
my $NetMHCIIpanWWWDIR = $ENV{'NetMHCIIpanWWWDIR'};
my $NetMHCIIpanWWWPATH = $ENV{'NetMHCIIpanWWWPATH'};

my $hist_limit = 20; #max number of histograms 

##############################################################################################################
# Define default values for optional variables and check if provided options and files are in the right format
##############################################################################################################

# When the help option is defined
if (defined $help) {
	&usage();
	exit;
}
#When the list option was defined
if (defined $list) {
	system ("cat $convert_HLA_list | cut -f2 -d ' '");
	exit;
}

unless (defined $file) {
	if (defined $w) {
		print "ERROR: no input sequences\n";
		exit();
	}
	&usage();
	exit();
}

# temporary directory
if (!defined $tdir) {
	$tdir = $ENV{'TMPDIR'} . "/tmp_$sId";
}
elsif (defined $tdir and $tdir eq ""){
	print "ERROR: insufficient arguments for option -tdir\nCheck ussage of the program using -h option\n";
	exit;
}

####################DEFINING SYNLISTS
my ($synlistNew);
if (!defined $context) {
	$synlistNew = "$rdir/data/synlist_nocontext.bin";
}
else {
	$synlistNew = "$rdir/data/synlist_context.bin";
}

#############BA logic

my $nout = 1;#nout is always 1

if (defined $BA){#pboth is only active if BA is selected
	$pboth = 1;
}

$exclude_offset = 1; #Offset is not defined in this version, should always be turned off

# input file format
my $input_format = "";
if (!defined $inptype) {
	$inptype = 0;
	$input_format = "FASTA";	
}
elsif ($inptype == 0) {
	$input_format = "FASTA";
}
elsif ($inptype == 1) {
	$input_format = "PEPTIDE";
}
unless ($inptype == 0 or $inptype == 1) {
	print "ERROR: input type should be 1 for peptide format and 0 for FASTA format\n";
	exit;
}

#peptide length
if (!defined $length and $inptype == 0) {
	$length = 15;
}
elsif (defined $length and $length eq "") {
	print "ERROR: insufficient arguments for option -length\nCheck ussage of the program using -h option\n";
	exit;
}
elsif (defined $length and $inptype == 0 and $length ne "") {
	if ($length < 9 ) {
		print "ERROR: peptide length should be 9 or more amino acids\n";
		exit;
	}	
}

# rank threshold for weak binders (default 10%)
if (!defined $thr_rank_W) {
	$thr_rank_W = 10;
}
unless (&isAnumber($thr_rank_W)) {
	print "ERROR: threshold rank value for weak binders should be a number\n";
	exit;
}

# rank threshold for strong binders (default 2%)
if (!defined $thr_rank_S) {
	$thr_rank_S = 2;
}
unless (&isAnumber($thr_rank_S)) {
	print "ERROR: threshold rank value for strong binders should be a number\n";
	exit;
}

# Checking if hla sequence was defined
my $hlaseq_option=0;
if (defined $hlaseq and $hlaseq ne "") {
	$hlaseq_option = 1;
	$histhr=100;
	my $first = 0;
	open (IN, "<", $hlaseq) or die "Can not open the file with MHC sequence $hlaseq $!\n";
	while (defined (my $line =<IN>)) {
		chomp $line;
		if ($line =~ m/^>/) {
			$first++;
		}
		elsif ($line !~ m/^>/ and $line !~ m/^[a-zA-Z]+$/) {
			print "ERROR: provided MHC sequence of beta chain is in a wrong format\n";
			exit;
		}
	}
	close IN;
	if ($first == 0) { 
		print "ERROR: provided MHC sequence of beta chain is not in FASTA format\n";
		exit;
	}
	# Checking if hla sequence was defined for alpha chain
	if (defined $hlaseq_a and $hlaseq_a ne "") {
		$hlaseq_option = 2;
		my $first = 0;
		open (IN, "<", $hlaseq_a) or die "Can not open the file with MHC sequence $hlaseq_a $!\n";
		while (defined (my $line =<IN>)) {
			chomp $line;
			if ($line =~ m/^>/) {
				$first++;
			}
			elsif ($line !~ m/^>/ and $line !~ m/^[a-zA-Z]+$/) {
				print "ERROR: provided MHC sequence of alpha chain is in a wrong format\n";
				exit;
			}
		}
		close IN;
		if ($first == 0) { 
			print "ERROR: provided MHC sequence of alpha chain is not in FASTA format\n";
			exit;
		}
	}
	elsif (defined $hlaseq_a and $hlaseq_a eq "") {
		print "ERROR: insufficient arguments for option -hlaseqA\nCheck usage of the program using -h option\n";
		exit;
	}	
}
elsif (defined $hlaseq_a and !defined $hlaseq) {
	print "ERROR: FASTA files for both MHC sequence chains should be specified, please enter MHC beta chain sequence in FASTA format\n";
	exit;
}		
elsif (defined $hlaseq and $hlaseq eq "") {
	print "ERROR: insufficient arguments for option -hlaseq\nCheck ussage of the program using -h option\n";
	exit;
}

# rank threshold for filtering output
unless (&isAnumber($rank_f)) {
	print "ERROR: threshold rank value for filtering should be a number\n";
	exit;
}

# Give the default values for rank and affinity

if (!defined $rank_f) {
	$rank_f = 10;
}

# cutoff value for core histograms
if (!defined $histhr) {
	$histhr = 100;
}
unless (&isAnumber($histhr)) {
	print "ERROR: threshold on %Rank for histograms should be a number\n";
	exit;
}

# File name for xls output
my $xlsfilename;
if (!defined $xlsfile) {
	$xlsfilename = "NetMHCIIpan_out.xls";
}
elsif (defined $xlsfile and $xlsfile eq "") {
	print "ERROR: insufficient arguments for option -xlsfile\nCheck ussage of the program using -h option\n";
	exit;
}
else {
	$xlsfilename = $xlsfile;
}

# Store allele list (one allele or a comma separated file with many alleles)
my @alleles_input = ();

if (! defined $alleles and ! defined $choose_chain and ! defined $hlaseq) {
	$alleles = ("DRB1_0101");
}
elsif (! defined $alleles and defined $choose_chain) {
	if (defined $chain_a and defined $chain_b) {
		$alleles = ("HLA-" . $chain_a . "-" . $chain_b);
	}
	elsif ((defined $chain_a) and (!defined $chain_b)) {
		print "ERROR: beta chain is not defined\n";
		exit;
	}
	elsif ((!defined $chain_a) and (defined $chain_b)) {
		print "ERROR: alpha chain is not defined\n";
		exit;
	}
	else {
		print "ERROR: alpha and beta chains are not defined \n";
		exit;
	}
}
elsif (defined $alleles and defined $choose_chain) {
	print "ERROR: you should only choose one option from '-a' and '-choose'\n";
	exit;
}
	
elsif (defined $alleles and $alleles eq "") {
	print "ERROR: insufficient arguments for option -a\nCheck ussage of the program using -h option\n";
	exit;
}

# converting allele names to the new nomenclature
my @alleles = ();
my %allele_names = ();

open(IN, "<", $convert_HLA_list) or die "Can't open the file for converting allele names $convert_HLA_list $! \n";
while (defined (my $line = <IN>)) {
	chomp $line;
	my @tmp = split (" ", $line);
	my $old_name = $tmp[1];
	my $new_name = $tmp[0];
	$allele_names{$old_name} = $new_name; ##3.2 is trained on simpler nomenclature
}
close IN;

# splitting the list of alleles separated by commas into array of alleles
if (! defined $hlaseq) {
	@alleles_input = split (",", $alleles);
	foreach my $a (@alleles_input) {				
		unless (exists $allele_names{$a}) {
			print "ERROR: Could not find allele $a - refer to the list of available molecules\n";
			exit;
		}	
		push @alleles, $a;
	}
}
else {
	@alleles = ("USER_DEF");
}

my $tot_molecules = $#alleles+1;

# input file	
if (! defined $file and ! defined $help) {
	print "ERROR: no input data\n";
	print "Usage: ./NetMHCIIpan-4.0.pl [-h] [args] -f [fastafile/peptidefile]\n";
	exit;
}
elsif (defined $file and $file eq ""){
	print "ERROR: insufficient arguments for option -f\nCheck usage of the program using -h option\n";
	exit;
}
elsif (defined $file) {
	unless (-e $file) {
		print "ERROR: input file $file doesn't exist\n";
		exit;
	}
}

# Testing if the input file is in FASTA format when no type specified and in peptide format when peptide type is specified
my %peplist;
my $flag_expected_aff = 0;

if ($inptype == 0) {
	open (IN, "<", $file) or die "Can not open input file $! \n";
	my $first = 0;
	while ($first == 0 and defined (my $input_line =<IN>)) {
		if ($input_line !~ m/^>/) {
			print "ERROR: Input file is not in FASTA format\n";
			exit;
		}
		$first = 1;
	}
	close IN;
}
# If input is in PEPTIDE format
elsif ($inptype == 1) {
	open (IN, "<", $file) or die "Can not open input file $! \n";
	while (defined (my $input_line =<IN>)) {
		chomp $input_line;
		next if length($input_line)<=1;
		
		my $pep = "";
		if ($input_line =~ m/(\S+)\s+(\S+)/) {
			$pep = $1;
			my $ann = $2;
			unless (isAnumber($ann)) {
				print "ERROR: Wrong format. The annotation on the second column must be numeric\n";
				exit;
			}
			$peplist{$pep} = $ann;
			$flag_expected_aff = 1;	
		} elsif ($input_line =~ m/(\S+)/) {
			$pep = $1;
			if ($flag_expected_aff == 1) {
				print "ERROR: Inconsistent format. Each line in the file should have the same format, either:\n1. Two columns with peptide and a numerical annotation or\n2. A single column of peptide sequences\n";
				exit;
			}	
			$flag_expected_aff = 0;
			$peplist{$pep} = "-99.999";
		}
		if ($pep !~ m/^[A-Za-z]+$/) {
				print "ERROR: Wrong format in input file - unknown charachers at line:\n$input_line\n";
				exit;
		}		
	}
	close IN;
}

my %synlist;
$synlist{"DR"} = "$synlists/DR_synlist";
$synlist{"DP"} = "$synlists/DP_synlist";
$synlist{"DQ"} = "$synlists/DQ_synlist";
$synlist{"H2"} = "$synlists/H-2_synlist";
$synlist{"ALL"} = "$synlists/ALL_synlist";

# creating temporary directory
mkdir $tdir, 0777;
unless (-d $tdir) {
	print "ERROR: directory $tdir doesn't exist\n";
	exit;
}

# Header flag
my $print_count = 0;

# Putting all the pseudosequences into the hash so that we have a look up table for each allele (used for ranks and etc.)
my %pseudoseqs = ();
open (IN, "<", $pseudoseq_file) or die "Can not open the file $pseudoseq_file $! \n";
while (defined (my $line =<IN>)) {
	chomp $line;
	my @tmp = split (" ", $line);
	my $allele = $tmp[0];
	my $seq = $tmp[1];
	$pseudoseqs{$allele} = $seq;
}
close IN;

# Defining hashes to be used for saving the results for making xls file
my %pos = ();
my %pep = ();
my %prot_id = ();
my %log = ();
my %rel = ();
my %nm = ();
my %bl = ();
my %rank = ();
my %targ = ();
my %histograms = ();
my $nhist=0;
#BA hashes for xls file
my %log_BA = ();
my %nm_BA = ();
my %rank_BA = ();


############################################################
### If MHC sequence has been uploaded - hlaseq option(s) ###
############################################################
if ($hlaseq_option != 0) {
	my $allele = "USER_DEF";
	my $rank_flag = 0;
	my $is_dpdq = 0;
	my $RESULT = "";
	
	# Getting pseudosequence for running the predictions
	my $pseudosequence = "";
	my $pseudoseqA = "";
	my $pseudoseqB = "";
	if ($hlaseq_option == 1) {
		$pseudoseqB = `cat $hlaseq | $mhcfsa2psseq -p $rdir/data/full_beta.positions -m $rdir/data/BLOSUM30 -gf 14 -r $rdir/data/DRB10101.fsa -n -- | grep -v '#' | cut -f2 -d ' '`;
		$pseudoseqA = `cat $rdir/data/DRA.pseudo`;
		chomp $pseudoseqA;
		chomp $pseudoseqB;		
	}
	elsif ($hlaseq_option == 2) {
		$pseudoseqA = `cat $hlaseq_a | $mhcfsa2psseq -p $rdir/data/full_alpha.positions -m $rdir/data/BLOSUM30 -gf 14 -r $rdir/data/DRA10101.fsa -n -- | grep -v '#' | cut -f2 -d ' '`;
		$pseudoseqB = `cat $hlaseq | $mhcfsa2psseq -p $rdir/data/full_beta.positions -m $rdir/data/BLOSUM30 -gf 14 -r $rdir/data/DRB10101.fsa -n -- | grep -v '#' | cut -f2 -d ' '`;
		chomp $pseudoseqA;
		chomp $pseudoseqB;		
	}
	
	if ($pseudoseqA eq "" or $pseudoseqB eq "") {
		print "ERROR: an error occured when getting pseudosequence from MHC fasta sequence\n";
		exit;
	}
	else {
		$pseudosequence = $pseudoseqA . $pseudoseqB;
	}
	if (length($pseudosequence) != 34) {
		print "ERROR: pseudosequence was not obtained correctly: $pseudosequence\n";
		exit;
	}
	
	my $new_pseudosequence = &changeGap($pseudosequence);	
	
	open (OUT, ">", "$tdir/pseudosequence_$$") or die "Can not open the file for writing '$tdir/pseudosequence_$$' $! \n";
	print OUT "USER_DEF $new_pseudosequence\n";
	close OUT;
	
	open (OUT, ">", "$tdir/allelelist_$$") or die "Can not open the file for writing '$tdir/allelelist_$$' $! \n";
	print OUT "USER_DEF USER_DEF\n";
	close OUT;
	##Use ALL synlist for user-defined pseudosequences
	my $synlist_file = $synlist{"ALL"};
	
	# Split FASTA files
	if ($inptype == 0) {
		my $FASTA_files = &SplitFiles($file);
		
		foreach my $fasta_file (@{$FASTA_files}) {
					
			# Preparing data to submit for predictions
			my ($identity, $input_file) = &PrepareInput($allele, $inptype, $fasta_file, $tdir);
	
			# Running the prediction script
			my $output_file = "$tdir/output_$$.out";
			my $cmd = "cat  $input_file | $nnalign_pan_player ";
			$cmd .= "-offset " unless defined($exclude_offset);
			$cmd .= "-context " unless !defined($context);
			$cmd .= "-pboth " unless !defined($pboth);
			#$cmd .= "-aX -afs -encmhc -elpfr 0 -eplen -1 -fl 3 -l 9 -nout $nout -classI 13,14,15,16,17,18,19,20,21 -mhc $pseudoseq_file -blf $blosum_freq -allelelist $allelelist -sbin $synlistNew -- | grep -v '#' > $output_file";
			$cmd .= "-aX -afs -encmhc -elpfr 0 -eplen -1 -fl 3 -l 9 -nout $nout -classI 13,14,15,16,17,18,19,20,21 -mhc $tdir/pseudosequence_$$ -blf $blosum_freq -allelelist $tdir/allelelist_$$ -sbin $synlistNew -- | grep -v '#' > $output_file";
			##system ($cmd);		
			doit( $cmd);
		
			# Accessing the results and preparing them to print as an output
			$RESULT = &AccessResults($allele, $output_file, $identity, $rdir, $thr_rank_S, $thr_rank_W, $rank_flag, 0);
			
			# Saving the results into hashes for creation of xls file
			if (defined $xls) {
				&populateXLShash($RESULT,$allele);
			}
	
			# Filtering the results if the option for printing only the best binding core was specified
			if (defined $unique) {
				$RESULT = &FindBestCore($RESULT);
			}
			
			# Modifying the results if the filtering or the sorting was used
			my ($FINAL_RESULT, $count_strong, $count_weak) = &ModifyResults($RESULT, $rank_f, $tdir, $sort, $rank_flag, $is_dpdq);			
		
			if (defined($hist) and defined($w)) {
				open (IN, "<", "$tdir/final_results.out") or die "Can not open the file $tdir/final_results.out for reading $!\n";

				while (defined (my $line = <IN>)) {										
						my @fields = split (" ", $line);
						my $all = $fields[1];
						my $pep = $fields[2];
						my $nm = $fields[8];

						if ($nhist<$hist_limit) {   #there is no %rank for custom MHCs
							$histograms{"$all-$pep"} = &make_histogram($all,$pep);
							$nhist++;
						}
				}				
				close IN;
			}
		
			# Printing the results
			&Print($allele, $is_dpdq, $FINAL_RESULT, $count_strong, $count_weak, $input_format, $pseudosequence, $thr_rank_S, $thr_rank_W, $rank_f, $rank_flag);	
		}
	}	
	else { 
		# Preparing data to submit for predictions
		my ($identity, $input_file) = &PrepareInput($allele, $inptype, $file, $tdir);
	
		# Running the prediction script
		
		my $output_file = "$tdir/output_$$.out";
		
		my $cmd = "cat  $input_file | $nnalign_pan_player ";
		$cmd .= "-offset " unless defined($exclude_offset);
		$cmd .= "-context " unless !defined($context);
		$cmd .= "-pboth " unless !defined($pboth);
		#$cmd .= "-aX -afs -encmhc -elpfr 0 -eplen -1 -fl 3 -l 9 -nout $nout -classI 13,14,15,16,17,18,19,20,21 -mhc $pseudoseq_file -blf $blosum_freq -allelelist $allelelist -sbin $synlistNew -- | grep -v '#' > $output_file";
		$cmd .= "-aX -afs -encmhc -elpfr 0 -eplen -1 -fl 3 -l 9 -nout $nout -classI 13,14,15,16,17,18,19,20,21 -mhc $tdir/pseudosequence_$$ -blf $blosum_freq -allelelist $tdir/allelelist_$$ -sbin $synlistNew -- | grep -v '#' > $output_file";
		##system ($cmd);	
		doit( $cmd);
		
		unlink "$tdir/pseudosequence_$$";
		unlink "$tdir/allelelist_$$";
	
		# Accessing the results and preparing them to print as an output
		$RESULT = &AccessResults($allele, $output_file, $identity, $rdir, $thr_rank_S, $thr_rank_W, $rank_flag, 0);
	
		# Saving the results into hashes for creation of xls file
		if (defined $xls) {
			&populateXLShash($RESULT,$allele);
		}
		
		# Modifying the results if the filtering or the sorting was used
		my ($FINAL_RESULT, $count_strong, $count_weak) = &ModifyResults($RESULT, $rank_f, $tdir, $sort, $rank_flag, $is_dpdq);

		if (defined($hist) and defined($w)) {
			open (IN, "<", "$tdir/final_results.out") or die "Can not open the file $tdir/final_results.out for reading $!\n";
			while (defined (my $line = <IN>)) {										
					my @fields = split (" ", $line);
					my $all = $fields[1];
					my $pep = $fields[2];
					my $nm = $fields[8];

					
#					if ($nm < $histhr) {  #only if predicted affinity is smaller than threshold
					if ($nhist<$hist_limit) {   #there is no %rank for custom MHCs
						$histograms{"$all-$pep"} = &make_histogram($all,$pep);
						$nhist++;
					}
					
			}				
			close IN;
		}		

		# Printing the results
		&Print($allele, $FINAL_RESULT, $count_strong, $count_weak, $input_format, $pseudosequence, $thr_rank_S, $thr_rank_W, $rank_f, $rank_flag);
	}
	unlink "$tdir/pseudosequence_$$";
	unlink "$tdir/allelelist_$$";
}	
	
##########################################
### If allele or allele list was given ###
##########################################
if ($hlaseq_option == 0){
	foreach my $allele (@alleles) {
		my $rank_flag = 1;
		my $RESULT = "";
		my $is_dpdq=0;
		
		my $synlist_file = $synlist{"DR"};
		if ($allele =~ m/DP/) {
			$synlist_file = $synlist{"DP"};
			$is_dpdq=1;
		} elsif ($allele =~ m/DQ/) {
			$synlist_file = $synlist{"DQ"};
			$is_dpdq=1;
		} elsif ($allele =~ m/H-2/) {
			$synlist_file = $synlist{"H2"};
		}
		# Split FASTA files
		if ($inptype == 0) {
			my $FASTA_files = &SplitFiles($file);
		
			foreach my $fasta_file (@{$FASTA_files}) {	
				# Preparing the data to submit for predictions
				my ($identity, $input_file) = &PrepareInput($allele, $inptype, $fasta_file, $tdir);
			
				# Running the prediction script
				my $output_file = "$tdir/output_$$.out";
				
				my $cmd = "cat  $input_file | $nnalign_pan_player ";
				$cmd .= "-offset " unless defined($exclude_offset);
				$cmd .= "-context " unless !defined($context);
				$cmd .= "-pboth " unless !defined($pboth);
				$cmd .= "-aX -afs -encmhc -elpfr 0 -eplen -1 -fl 3 -l 9 -nout $nout -classI 13,14,15,16,17,18,19,20,21 -mhc $pseudoseq_file -blf $blosum_freq -allelelist $allelelist -sbin $synlistNew -- | grep -v '#' > $output_file";
				##system ($cmd);	
				doit( $cmd);
		
				# Accessing the results and preparing them from printing 
				$RESULT = &AccessResults($allele, $output_file, $identity, $rdir, $thr_rank_S, $thr_rank_W, $rank_flag, \%pseudoseqs);
			
				# Saving the results into hashes for creation of xls file
				if (defined $xls) {
					&populateXLShash($RESULT,$allele);
				}
	
				# Filtering the results if the option for printing only the best binding core was specified
				if (defined $unique) {
					$RESULT = &FindBestCore($RESULT);
				}		
				
				# Modifying the results if the filtering or the sorting was used
				my ($FINAL_RESULT, $count_strong, $count_weak) = &ModifyResults($RESULT, $rank_f, $tdir, $sort, $rank_flag, $is_dpdq);
		
				if (defined($hist) and defined($w)) {
					open (IN, "<", "$tdir/final_results.out") or die "Can not open the file $tdir/final_results.out for reading $!\n";

					while (defined (my $line = <IN>)) {										
						my @fields = split (" ", $line);
						my $all = $fields[1];
						my $pep = $fields[2];
						my $nm = $fields[8];
						my $rk = $fields[9];
						next if $rk eq 'NA';
						
#						if ($nm < $histhr) {  #only if predicted affinity is smaller than threshold
						if ($rk <= $histhr) {   #only if predicted rank is smaller than threshold
							last if $nhist>=$hist_limit;
							$histograms{"$all-$pep"} = &make_histogram($all,$pep);
							$nhist++;
						}
					}				
					close IN;
				}			
		
				# Printing the results
				&Print($allele, $is_dpdq, $FINAL_RESULT, $count_strong, $count_weak, $input_format, 0, $thr_rank_S, $thr_rank_W, $rank_f, $rank_flag);
			}
		}
		else {
			# Preparing the data to submit for predictions
			my ($identity, $input_file) = &PrepareInput($allele, $inptype, $file, $tdir);
		
			# Running the prediction script
			my $output_file = "$tdir/output_$$.out";
			
			my $cmd = "cat  $input_file | $nnalign_pan_player ";
			$cmd .= "-offset " unless defined($exclude_offset);
			$cmd .= "-context " unless !defined($context);
			$cmd .= "-pboth " unless !defined($pboth);
			$cmd .= "-aX -afs -encmhc -elpfr 0 -eplen -1 -fl 3 -l 9 -nout $nout -classI 13,14,15,16,17,18,19,20,21 -mhc $pseudoseq_file -blf $blosum_freq -allelelist $allelelist -sbin $synlistNew -- | grep -v '#' > $output_file";
			##system ($cmd);	
			doit( $cmd);
		
			# Accesing the results and preparing them to print as an output
			$RESULT = &AccessResults($allele, $output_file, $identity, $rdir, $thr_rank_S, $thr_rank_W, $rank_flag, \%pseudoseqs);
			
			# Saving the results into hashes for creation of xls file
			if (defined $xls) {
				&populateXLShash($RESULT,$allele);
			}
		
			# Modifying the results if the filtering or the sorting was used
			my ($FINAL_RESULT, $count_strong, $count_weak) = &ModifyResults($RESULT, $rank_f, $tdir, $sort, $rank_flag, $is_dpdq);			
		
			if (defined($hist) and defined($w)) {
				open (IN, "<", "$tdir/final_results.out") or die "Can not open the file $tdir/final_results.out for reading $!\n";

				while (defined (my $line = <IN>)) {										
						my @fields = split (" ", $line);
						my $all = $fields[1];
						my $pep = $fields[2];
						my $nm = $fields[8];
						my $rk = $fields[9];
						next if $rk eq 'NA';
						
#						if ($nm < $histhr) {  #only if predicted affinity is smaller than threshold
						if ($rk <= $histhr) {   #only if predicted rank is smaller than threshold
							last if $nhist>=$hist_limit;
							$histograms{"$all-$pep"} = &make_histogram($all,$pep);
							$nhist++;
						}

				}				
				close IN;
			}			
			# Printing the results
			&Print($allele, $is_dpdq, $FINAL_RESULT, $count_strong, $count_weak, $input_format, 0, $thr_rank_S, $thr_rank_W, $rank_f, $rank_flag);
		}
	}
}

#######################################
### Preparing XLS file if specified ###
#######################################
# Acessing and printing  the results into the tab separated file for .xls output if the option "xls" option was specified
if (defined $xls) {
	foreach my $allele1 (@alleles) {
		foreach my $allele2 (@alleles) {
			if ($#{$pos{$allele1}} != $#{$pos{$allele2}}) {
				print "ERROR occured when creating Excel file: $allele1 and $allele2 resulted into different number of peptides\n";
				exit;
			}
		}
	}		
	my $first = "FALSE";
	my $file_name;
	if (defined $w) {
		$file_name = "$NetMHCIIpanWWWDIR/$$"."_NetMHCIIpan.xls";
	}
	else {
		$file_name = $xlsfilename;
	}
	open (OUT, ">", $file_name) or die "Can not open the file $!\n";
	foreach my $query_allele (@alleles) {
		if ($first eq "FALSE") {
			print OUT "\t\t\t\t$query_allele";
			$first = "TRUE";
		}
		else {
			if (defined $BA){
				print OUT "\t\t\t\t\t$query_allele";
			}
			else{
				print OUT "\t\t$query_allele";
			}

		}
	}
	$first = "FALSE";
	foreach my $query_allele (@alleles) {
		if ($first eq "FALSE") {
			if (!exists $rank{$alleles[0]}) {
				if (defined $BA){
					print OUT "\nPos\tPeptide\tID\tTarget\tScore\tScore_BA\tnM";
				}
				else{
					print OUT "\nPos\tPeptide\tID\tTarget\tScore";	
				}
			}
			else {
				if (defined $BA){
					print OUT "\nPos\tPeptide\tID\tTarget\tScore\tRank\tScore_BA\tnM\tRank_BA";
				}
				else{
					print OUT "\nPos\tPeptide\tID\tTarget\tScore\tRank";
				}
			}			
			$first = "TRUE";
		}
		else {
			if (defined $BA){
				print OUT "\tScore\tRank\tScore_BA\tnM\tRank_BA";
			}
			else{
				print OUT "\tScore\tRank";
			}			
		}
	}
	print OUT "\tAve\tNB\n";

	for (my $i = 0; $i <= $#{$pos{$alleles[0]}}; $i++) {
		my $nb = $bl{$alleles[0]}[$i];
		my $log_sum = $log{$alleles[0]}[$i];
		if (!exists $rank{$alleles[0]}) {
			if (defined $BA){
				print OUT "$pos{$alleles[0]}[$i]\t$pep{$alleles[0]}[$i]\t$prot_id{$alleles[0]}[$i]\t$targ{$alleles[0]}[$i]\t$log{$alleles[0]}[$i]\tlog_BA{$alleles[0]}[$i]\tnm_BA{$alleles[0]}[$i]\t";
			}
			else{
				print OUT "$pos{$alleles[0]}[$i]\t$pep{$alleles[0]}[$i]\t$prot_id{$alleles[0]}[$i]\t$targ{$alleles[0]}[$i]\t$log{$alleles[0]}[$i]\t";
			}			
		}
		else {
			if (defined $BA){
				print OUT "$pos{$alleles[0]}[$i]\t$pep{$alleles[0]}[$i]\t$prot_id{$alleles[0]}[$i]\t$targ{$alleles[0]}[$i]\t$log{$alleles[0]}[$i]\t$rank{$alleles[0]}[$i]\t$log_BA{$alleles[0]}[$i]\t$nm_BA{$alleles[0]}[$i]\t$rank_BA{$alleles[0]}[$i]\t";
			}
			else{
				print OUT "$pos{$alleles[0]}[$i]\t$pep{$alleles[0]}[$i]\t$prot_id{$alleles[0]}[$i]\t$targ{$alleles[0]}[$i]\t$log{$alleles[0]}[$i]\t$rank{$alleles[0]}[$i]\t";
			}			
		}
		for (my $n = 1; $n <= $#alleles; $n++){
			$nb += $bl{$alleles[$n]}[$i];
			$log_sum += $log{$alleles[$n]}[$i];
			if (defined $BA){
				print OUT "$log{$alleles[$n]}[$i]\t$rank{$alleles[$n]}[$i]\t$log_BA{$alleles[$n]}[$i]\t$nm_BA{$alleles[$n]}[$i]\t$rank_BA{$alleles[$n]}[$i]\t";
			}
			else{
				print OUT "$log{$alleles[$n]}[$i]\t$rank{$alleles[$n]}[$i]\t";
			}			
		}
		my $avg = $log_sum/scalar(@alleles);
		$avg = sprintf("%.4f", $avg);
		print OUT "$avg\t$nb\n";
	}
	close OUT;
	my $short_name = "$NetMHCIIpanWWWPATH/$$"."_NetMHCIIpan.xls";
	if (defined $w) {
		print "Link to output xls file <a href='$short_name'>NetMHCIIpan_out.xls</a>\n";
	}	
}
#################### Finished creating .xls file #################################################

# Deleting temporary directory if the dirty mode was not chosen
if (!defined $dirty) {
	system ("rm -r $tdir");
}
else {
	print "\n\nTemporary files have been saved here: $tdir\n";
}
###################################
###   S U B R O U T I N E S	###
###################################

## system calls
sub doit {

	
	if ( defined($verbose) ) {
		print "# doit: $_[0]\n";
	}

	system( $_[0] );
	
}

## check for integer
sub isInt {
    my $test=shift;
    
    if ($test =~ m/^\d+$/ && $test>=0) {
	return 1; }
    else {
	return 0; }
}

## check format of a number
sub isAnumber {
    my $test = shift;

    eval {
        local $SIG{__WARN__} = sub {die $_[0]};
        $test += 0;
    };
    if ($@) {
	return 0;}
    else {
	return 1;} 
}

## if there is a gap in a pseudosequence, change the position of it
sub changeGap {
    my $seq = shift;
    my $new_seq = "";
    my $no_gaps = 0;
    for (my $i = 0; $i < length($seq); $i++) {
    	my $symbol = substr ($seq, $i, 1);
	if ($symbol eq "-") {
		$no_gaps = 1;
	}
	
    }
    if ($no_gaps == 1) {
    	for (my $i = 0; $i < length($seq); $i++) {
    		my $symbol = substr ($seq, $i, 1);
		if ($i < 6 or $i > 6 ){
			if ($symbol ne "-") {
				$new_seq .= $symbol;
			}
		}
		elsif ($i == 6) {
			$new_seq .= "$symbol-";			
		}
	}
    }
    else {	
    	$new_seq = $seq;
}	
    return $new_seq;
}

sub SplitFiles {
	my ($file) = @_;
	my $flag = 0;
	my @fasta_files_array = ();
	my $fasta_file;
	open (IN, "<", $file) or die "Can not open input file $file for reading $! \n";
	my $unique_id_count = 0;
	while (defined (my $line = <IN>)) {
		chomp $line;
		
		if ($line =~ m/^>(\S+)\s*/) {## MN has changed  m/^>(\S+)\s+/ to  m/^>(\S+)\s*/	04092013
			$unique_id_count++;
			if ($flag == 0) {
				$fasta_file = $tdir ."/". $1 . "_" . $unique_id_count;
				open (OUT, ">", $fasta_file) or die "Can not open input file $fasta_file for writing $! \n";
				print OUT "$line\n";
				$flag++;
			}
			else {
				close OUT;
				push (@fasta_files_array, $fasta_file);
				$fasta_file = $tdir . "/". $1 . "_" . $unique_id_count;
				open (OUT, ">", $fasta_file) or die "Can not open input file $fasta_file for writing $! \n";
				print OUT "$line\n";	
			}
						
		}
		if ($line !~ m/^>/) {
			print OUT "$line\n";
		}		
	}
	close OUT;
	push (@fasta_files_array, $fasta_file);
	return (\@fasta_files_array);		
}
sub PrepareInput {
	my ($allele, $inptype, $file, $tdir) = @_;
	my $identity = "Sequence";
	# If the input was in PEPTIDE format
	if ($inptype == 1) {
		open (OUT, ">", "$tdir/input_file.dat") or die "Can not open the file '$tdir/input_file.dat' for writing $!\n";
		open (IN, "<", $file) or die "Can not open the file $file for reading $! \n";
		my @tmp = ();
		my $pep = 0;
		my $aff = 0;
		my $con = 0;
		while (defined (my $line = <IN>)) {
			chomp $line;
			if ($line =~ m/(\S+)\s+(\S+)/ and defined $context) {
				@tmp = split (" ", $line);
				$pep = $tmp[0];
				$aff = $tmp[1];
				$con = $tmp[2];
				if (length($pep) < 9 ) {
					print "ERROR: peptide length should be 9 or more amino acids: $line\n";
					exit;
				}
				print OUT "$pep\t$aff\t$allele\t$con\n";				
			}
			elsif ($line =~ m/(\S+)\s+(\S+)/) {
				my $pep = $1;
				my $aff = $2;
				if (length($pep) < 9 ) {
					print "ERROR: peptide length should be 9 or more amino acids: $line\n";
					exit;
				}
				print OUT "$pep\t$aff\t$allele\n";				
			}
			elsif ($line =~ m/(\S+)/) {
				my $pep = $1;
				if (length($pep) < 9 ) {
					print "ERROR: peptide length should be 9 or more amino acids\n";
					exit;
				}
				print OUT "$pep\t-99.999\t$allele\n";
			}
			
		}
		close IN;
		close OUT;
	}	
	elsif ($inptype == 0 and !defined $context) {
		my $seq = "";
		open (IN, "<", $file) or die "Can not open input file $file for reading $! \n";
		while (defined (my $line = <IN>)) {
			chomp $line;
			if ($line =~ m/^>(\S+)/) {
				$identity = $1;
			} else {
			   $seq .= $line;
			}
		}
		close IN;
		open (OUT, ">", "$tdir/input_file.dat") or die "Can not open the file '$tdir/input_file.dat' for writing $!\n";
		for (my $i=0; $i < length($seq) - $length + 1; $i++) {
			my $pep = substr ($seq, $i, $length);
			print OUT "$pep\t-99.999\t$allele\n";
		}
		close OUT;
	}
	elsif ($inptype == 0 and defined $context) {
		my $seq = "";
		open (IN, "<", $file) or die "Can not open input file $file for reading $! \n";
		while (defined (my $line = <IN>)) {
			chomp $line;
			if ($line =~ m/^>(\S+)/) {
				$identity = $1;
			} else {
			   $seq .= $line;
			}
		}
		close IN;
		my $pfr_lgt = 3;
		open (OUT, ">", "$tdir/input_file.dat") or die "Can not open the file '$tdir/input_file.dat' for writing $!\n";
		for (my $i=0; $i < length($seq) - $length + 1; $i++) {
			my $pep = substr ($seq, $i, $length);

	        my $nstart = $i-$pfr_lgt > 0 ? $i-$pfr_lgt : 0;
    	    my $nterm_L = $i>=$pfr_lgt ? $pfr_lgt : $i;
    
	        my $cstart = $i+$length;   
    	    my $cterm_L = length($seq)-$pfr_lgt>$cstart ? $pfr_lgt : length($seq)-$cstart;

	        my $nterm = substr($seq, $nstart, $nterm_L);
    	    my $cterm = substr($seq, $cstart, $cterm_L);
            my $nterm_lig = substr($pep,0,$pfr_lgt);
            my $cterm_lig = substr($pep,$length-$pfr_lgt,$pfr_lgt);
        	##Add gaps to complete the window
	        #$nterm = 'X' x ($pfr_lgt-$nterm_L) . $nterm;
    	    #$cterm = $cterm . 'X' x ($pfr_lgt-$cterm_L);    
    	    #Patch to account for C-terminal ligands
	        $nterm = 'X' x ($pfr_lgt-$nterm_L) . $nterm;
    	    $cterm = $cterm . 'X' x ($pfr_lgt-$cterm_L);        	    
        
        	my $context_out = $nterm . $nterm_lig . $cterm_lig . $cterm;
			print OUT "$pep\t-99.999\t$allele\t$context_out\n";
		}
		close OUT;
	}	

	return ($identity, "$tdir/input_file.dat");
}

sub getRank {
	my ($Score, $thresholdFile,$allele_sequence) = @_;
	my @RANKS = ();
	my @SCORES = ();
	my $rank=100.0;
	open (IN1, "<", "$thresholdFile") or die "Can not open the file $rdir/data/thresholds/$allele_sequence.thr_cmb $! \n";
	while (defined (my $line =<IN1>)) {				
		chomp $line;
		next if $line=~m/^#/;
		my @tmp = split (" ", $line);
		push (@RANKS, $tmp[0]);
		push (@SCORES, $tmp[2]);
	}
	close IN1;
	my $flag = 0;
	for (my $i = 0; $i <= $#RANKS; $i++) {
		#print($Score);
		if ($Score >= $SCORES[$i] and $flag == 0) {
			$flag = 1;
			#$rank = $RANKS[$i];
			#interpolation of ranks: y = y_1 + (y_2-y_1)*/((x-x_1)/(x_2-x_1))
			#print "Rank: $RANKS[$i]";
			#print "Score: $SCORES[$i]";
			if ($i==0){#If score is larger than largest score in rank table, exception
				$rank=0.0;
			}
			else{
				$rank = $RANKS[$i-1]+($RANKS[$i]-$RANKS[$i-1])*(($Score-$SCORES[$i-1])/($SCORES[$i]-$SCORES[$i-1]));
			}
		}
		if ($i == $#RANKS and $Score < $SCORES[$i]) {
			$rank = $RANKS[$#RANKS];
		}
	}
	#print("$rank\n");
	$rank = sprintf("%6.2f", $rank);
	return ($rank);

}

sub AccessResults {
	my ($allele, $output_file, $identity, $rdir, $thr_rank_S, $thr_rank_W , $rank_flag, $pseudoseqs) = @_;
	my $RESULT = "";
	my $pos=1;
	my $rank = "NA";
	open (IN, "<", $output_file) or die "Can not open the file $output_file for reading: $!";
	while (defined (my $line = <IN>) ) {
		chomp $line;
		if ($line =~ m/Error\. Cannot find MHC name/) {
			print "ERROR: Could not find MHC name, see the list of possible molecules\n";
			exit;
		}
		if ($line =~ m/Error\. Wrong line format/) {
			print "ERROR: Input contains lines in a wrong format. Check if you have empty lines\n";
			exit;
		}
		if ($line =~ m/Error\. Not elements in trainlist/) {
			print "ERROR: Wrong input format.\n";
			exit;
		}
		my @tmp = split (" ", $line);
		my $core = $tmp[0];
		my $offset = $tmp[1];
		my $score = $tmp[3];
		my $peptide = $tmp[4];
		my $reliability = $tmp[10];
		my $output_allele = $tmp[11];
		my $score_BA = 0.0;
		if (defined $BA){
			if ( defined $context ) {
				$score_BA = $tmp[15];
			}	
			else {
				$score_BA = $tmp[14];
			}
		}
		
		# Finding the position of the core in the peptide (unnecessary? use $offset?)
		my $length = length($peptide);
		
		my $pos_core = $offset;
		
		my $expected = exists($peplist{$peptide}) ? $peplist{$peptide} : -99.999;
		$expected = sprintf("%7.3f", $expected);
			
		# Calculating affinity from the log score
		##my $aff_EL = "NA";
		my $aff_BA = 0.0;
		if (defined $BA){
			$aff_BA = exp((1-$score_BA)*log(50000));	
		}
		
		$score = sprintf("%13.9f", $score);
		my $rank_EL = 0.0;
		my $rank_BA = 0.0;
		# Finding the rank if allele name was given and not MHC sequence
		if ($rank_flag == 1) {
			
			#my @RANKS = ();
			#my @SCORES = ();
			
			my $allele_sequence = $pseudoseqs{$allele};

			######Thresold logic - Select the correct threshold file based on the input parameters nout and context
			#my $threshold_file;
			my $threshold_file_EL;
			my $threshold_file_BA;
			if (!defined $context){
				$threshold_file_BA = "$rdir/data/netMHCIIpan-4.0_thresholds/$allele_sequence.thr.ba";
				$threshold_file_EL = "$rdir/data/netMHCIIpan-4.0_thresholds/$allele_sequence.thr.el";
			}
			else {
				$threshold_file_BA = "$rdir/data/netMHCIIpan-4.0_thresholds/$allele_sequence.thr.bacont";
				$threshold_file_EL = "$rdir/data/netMHCIIpan-4.0_thresholds/$allele_sequence.thr.elcont";				
			}
			
			$rank_EL = getRank($score,$threshold_file_EL,$allele_sequence);
			if (defined $BA){
				$rank_BA = getRank($score_BA,$threshold_file_BA,$allele_sequence);
			}
			



		}
		else {
			$rank_EL = -99;
			if (defined $BA){
				$rank_BA = -99;
			}
			
		}
	
		## Finding the level of binding
		my $level ="";
		if ($rank_flag == 1) {
			if ($rank_EL <= $thr_rank_S) {
				$level = "<=SB";
			}
			elsif ($rank_EL <= $thr_rank_W) {
				$level = "<=WB";
			}
		}
		
		#print("$level\n");

		# Saving the results into one variable
		if (defined $BA){
			$RESULT .= "$pos $output_allele $peptide $identity $pos_core $core $reliability $score $rank_EL $expected $score_BA $aff_BA $rank_BA $level\n";
	    	}
	    	else{
	    		$RESULT .= "$pos $output_allele $peptide $identity $pos_core $core $reliability $score $rank_EL $expected $level\n";
	    	}	
	    	$pos++;
	}	
	return $RESULT;
}

sub ModifyResults {
	my ($RESULT, $rank_f, $tdir, $sort, $rank_flag, $is_dpdq) = @_;
	my $RESULT_MOD = "";
	my $count_strong = 0;
	my $count_weak = 0;
	my $FINAL_RESULT = "";

	my @result_lines = split ("\n", $RESULT);
	foreach my $result_line (@result_lines) {
		my @scores = split (" ", $result_line);
		my $pos = sprintf("%4s", $scores[0]);
		my $output_allele = $is_dpdq ? sprintf("%23s", $scores[1]) : sprintf("%13s", $scores[1]);
		my $peptide = sprintf("%20s", $scores[2]);
		my $identity = sprintf("%15.15s", $scores[3]);
		my $pos_core = sprintf("%4s", $scores[4]);		
		my $core = sprintf("%11s", $scores[5]);
		my $rel = sprintf("%9.3f", $scores[6]);
		my $score = sprintf("%13.6f", $scores[7]);
		#my $aff = sprintf("%4s", $scores[8]);
		#my $aff = 0;
		#my $aff = -99.9;
		#if ($nout==0){
		#	$aff = sprintf("%12.2f", $scores[8]);
		#}
		#elsif ($nout==1){
		#	$aff = $aff;
		#}
		 
		#my $rank = $scores[9] eq "NA" ? sprintf("%8s", $scores[9]) : sprintf("%8.2f", $scores[9]);
		my $rank = $scores[8] eq "NA" ? sprintf("%8s", $scores[8]) : sprintf("%8.2f", $scores[8]);
		my $expected = $scores[9] eq "-99.999" ? sprintf("%8s", "NA") : sprintf("%8.3f", $scores[9]);
		my $score_BA = 0.0;
		my $aff_BA = 0.0;
		my $rank_BA = 0.0;
		my $level = "";
		if (defined $BA){
			$score_BA = sprintf("%13.6f", $scores[10]);
			$aff_BA = sprintf("%13.2f", $scores[11]);
			$rank_BA = sprintf("%8.2f", $scores[12]);
			$level = defined $scores[13] ? $scores[13] : "";	
		
		}
		else{
			$level = defined $scores[10] ? $scores[10] : "";	
		}
		#my $level = defined $scores[11] ? $scores[11] : "";
		
		# Finding the number of strong and weak binders
		if ($level eq "<=SB") {
			$count_strong++;			
		}
		elsif ($level eq "<=WB") {
			$count_weak++;			
		}
		
		#some formatting
		#$score = sprintf("%13s", $score);		
		#$score_BA = sprintf("%13s", $score_BA);
		#$aff = sprintf("%13s", $aff);
		my $level_formatted =  sprintf("%6s", $level);
		
		#If the filter for filtering output was defined
		if (defined($filter) && $rank_flag==1) {
			if ($rank <= $rank_f) {

				if (defined $BA){
					$RESULT_MOD .= "$pos $output_allele $peptide $pos_core $core $rel $identity $score $rank $expected $score_BA $aff_BA $rank_BA $level_formatted\n";	
				}
				else{
					$RESULT_MOD .= "$pos $output_allele $peptide $pos_core $core $rel $identity $score $rank $expected $level_formatted\n";	
				}

				
			}
			if ($level eq "<=SB" and $rank > $rank_f) {
				$count_strong--;
			}
			if ($level eq "<=WB" and $rank > $rank_f) {
				$count_weak--;
			}	
		}
		else {
				if (defined $BA){
					$RESULT_MOD .= "$pos $output_allele $peptide $pos_core $core $rel $identity $score $rank $expected $score_BA $aff_BA $rank_BA $level_formatted\n";	
				}
				else{
					$RESULT_MOD .= "$pos $output_allele $peptide $pos_core $core $rel $identity $score $rank $expected $level_formatted\n";	
				}			
		}
	}
		
	# Printing result lines into the file
	open (OUT, ">", "$tdir/results.out") or die "Can not open the file $tdir/results.out for writing $!\n";
	print OUT $RESULT_MOD;
	close OUT;

	# If the sort option was specified, sort the results based on affinity	
	if (defined $sort)  {
		system("cat $tdir/results.out | sort -nrk8 > $tdir/final_results.out");
		open (IN, "<", "$tdir/final_results.out") or die "Can not open the file $tdir/final_results.out for reading $!\n";
		while (defined (my $line = <IN>)) {			
			$FINAL_RESULT .= $line;
		}			
		close IN;
	}
	else {
		$FINAL_RESULT = $RESULT_MOD;
		system("cp $tdir/results.out $tdir/final_results.out");
	}
	return ($FINAL_RESULT, $count_strong, $count_weak);
}

sub populateXLShash{
	my ($RESULT,$allele) = @_;
	my @lines = split ("\n", $RESULT);		
	foreach my $result_line (@lines) {
		chomp $result_line;
		my @tmp = split (" ", $result_line);			
		push (@{$pos{$allele}}, $tmp[0]);
		push (@{$pep{$allele}}, $tmp[2]);
		push (@{$prot_id{$allele}}, $tmp[3]);
		push (@{$rel{$allele}}, $tmp[6]);
		push (@{$log{$allele}}, $tmp[7]);
		push (@{$nm{$allele}}, "NA");
		push (@{$rank{$allele}}, $tmp[8]);
		if ($input_format eq 'PEPTIDE'){
			push (@{$targ{$allele}}, $tmp[9]);	
		}
		elsif ($input_format eq 'FASTA'){
			push (@{$targ{$allele}}, "NA");	
		}
		if ($tmp[-1] eq "<=WB" || $tmp[-1] eq "<=SB") {
			push (@{$bl{$allele}}, 1);
		}
		else {
			push (@{$bl{$allele}}, 0);
		}
		if (defined $BA){
			push (@{$log_BA{$allele}}, $tmp[10]);
			push (@{$nm_BA{$allele}}, $tmp[11]);
			push (@{$rank_BA{$allele}}, $tmp[12]);
		}
	}

}

sub Print {
	my ($allele, $is_dpdq, $FINAL_RESULT, $count_strong, $count_weak, $input_format, $pseudosequence, $thr_rank_S, $thr_rank_W, $rank_f, $rank_flag) = @_;
	# defining the names for the initial lines to print 
	my $pos_print = sprintf("%4s", "Pos");
	my $pos_core_print = sprintf("%4s", "Of");
	my $allele_print = $is_dpdq ? sprintf("%23s", "MHC") : sprintf("%13s", "MHC"); 
	my $peptide_print = sprintf("%20s", "Peptide");
	my $core_print = sprintf("%11s", "Core");
	my $reliability_print = sprintf("%9s", "Core_Rel");
	my $identity_print = sprintf("%15.15s", "Identity");
	#my $score_print = sprintf("%13s", "1-log50k(aff)");
	my $score_print = sprintf("%13s", "Score_EL");
	#my $affinity_print = sprintf ("%13s", "Affinity(nM)");
	my $rank_print = sprintf("%8s", "\%Rank_EL");
	my $score_print_BA = sprintf("%13s", "Score_BA");
	my $affinity_print_BA = sprintf ("%13s", "Affinity(nM)");
	my $rank_print_BA = sprintf("%8s", "\%Rank_BA");
	my $level_print = " BindLevel";
	my $expected_affinity_print = sprintf("%8s", "Exp_Bind");
	
	if ($print_count == 0) {
		print "# NetMHCIIpan version $version\n\n" , 
      		"# Input is in $input_format format\n\n";
		
		if ($inptype ==0) {
			print "# Peptide length $length\n\n";
		}
		if (defined $BA and !defined $context){
			print "# Prediction Mode: EL+BA\n\n"
		}
		if (!defined $BA and !defined $context){
			print "# Prediction Mode: EL\n\n"
		}
		if (defined $BA and defined $context){
			print "# Prediction Mode: EL+BA with Context\n\n"
		}
		if (!defined $BA and defined $context){
			print "# Prediction Mode: EL with Context\n\n"
		}
		if ($rank_flag == 0) {
			print "# User-defined MHC pseudo sequence $pseudosequence USER_DEF\n\n";
		}

		if ($rank_flag == 1) {
			print "# Threshold for Strong binding peptides (\%Rank)\t$thr_rank_S%\n",
			"# Threshold for Weak binding peptides (\%Rank)\t$thr_rank_W%\n";
		} else {
			print "# NB: Strong and Weak binders can not be defined with a user-defined MHC pseudo-sequence.\n"
		}
		if (defined($filter)) {
			print "\n# Threshold for filtering output (\%Rank)\t$rank_f%\n",
		}
      		
		$print_count++;
	}
		
	# Printing the results for each allele
	print "\n# Allele: $allele\n";	

	print '-' x 140; 
	if (defined $BA){
		print "\n$pos_print $allele_print $peptide_print $pos_core_print $core_print $reliability_print $identity_print $score_print $rank_print $expected_affinity_print $score_print_BA $affinity_print_BA $rank_print_BA $level_print\n";
	}
	else{
		print "\n$pos_print $allele_print $peptide_print $pos_core_print $core_print $reliability_print $identity_print $score_print $rank_print $expected_affinity_print $level_print\n";
	}
	print '-' x 140; 
	print "\n";

	open (IN, "<", "$tdir/final_results.out") or die "Can not open the file $tdir/final_results.out for reading $!\n";	
	while (defined (my $line = <IN>)) {	
			chomp $line;				
			my @fields = split (" ", $line);
			my $all = $fields[1];
			my $pep = $fields[2];
			if (defined $w && exists($histograms{"$all-$pep"})) {
				my $link = $histograms{"$all-$pep"};				
				$line .= "\t<a href=\"$link\" target='_blank'>Core_Histogram</a>";
			}
			print "$line\n";
	}
	close IN;
	
	print '-' x 140; 
	print "\nNumber of strong binders: $count_strong Number of weak binders: $count_weak\n";
	print '-' x 140; 
	print "\n";
}

sub FindBestCore {
	my $RESULTS_ALL= $_[0];
	my $RESULT = "";
	my @result_lines = split("\n", $RESULTS_ALL);
	my @uniq_lines;
	my $count_strong = 0;
	my $count_weak = 0;
	##first line	
	my @firstline = split(" ", $result_lines[0]);
	my $core_pos_prev = $firstline[4];
	my $core_prev = $firstline[5];
	my $score_best = $firstline[7];	
	my $best_line = $result_lines[0];
	##remaining lines
	for (my $n=1; $n<=$#result_lines; $n++) {
		
		my $line = $result_lines[$n];
		my @line = split( " ", $line);
		my $core_pos = $line[4];
		my $core = $line[5];
		my $score = $line[7];
	
		if (($core eq $core_prev) and ($core_pos == $core_pos_prev-1))  {
			if ($score > $score_best) {
			   $best_line = $line;
			   $score_best = $score;
			}
		} else {
		   push (@uniq_lines, $best_line);
		   $best_line = $line;
		   $score_best = $score;	
		}
		$core_pos_prev = $core_pos;
		$core_prev = $core;
	}
	#last line	
	push (@uniq_lines, $best_line);

	$RESULT = join("\n", @uniq_lines);
	$RESULT .= "\n";
	return ($RESULT);
}

sub make_histogram {
	my $allele = $_[0];
	my $allele_nn = exists($allele_names{$allele}) ? $allele_names{$allele} : $allele;
	my $peptide = $_[1];
	my $pdf = "$NetMHCIIpanWWWDIR/$$"."_$allele.$peptide.pdf";
	
	if (defined($hist) && (length($R)<1  || $R =~ /not found/)) {
		print "#Warning. R-2.14 does not seem to be installed on this system, graphics won't be visualized.";
		return "0";
	}	
	
	##choose synlist
	
	my $synlist_file = $synlist{"DR"};
	if ($allele =~ m/DP/) {
		$synlist_file = $synlist{"DP"};
	} elsif ($allele =~ m/DQ/) {
		$synlist_file = $synlist{"DQ"};
	} elsif ($allele =~ m/H-2/) {
		$synlist_file = $synlist{"H2"};
	}

	my $resfile = "$tdir/$allele.$peptide.printall.txt";
	my $histfile = "$tdir/$allele.$peptide.hist.txt";

	my $pseudo = $hlaseq_option == 0 ? $pseudoseq_file : "$tdir/pseudosequence_$$";
	#my $output_file = "$tdir/output_$$.out";
	my $cmd = "echo $peptide | $nnalign_pan_player ";
	$cmd .= "-offset " unless defined($exclude_offset);	
	$cmd .= "-context " unless !defined($context);
	$cmd .= "-pboth " unless !defined($pboth);
	$cmd .= "-printall -aX -afs -encmhc -elpfr 0 -eplen -1 -fl 3 -l 9 -nout $nout -classI 13,14,15,16,17,18,19,20,21 -mhc $pseudoseq_file -blf $blosum_freq -allelelist $allelelist -sbin  $synlistNew -- | grep -v '#' > $resfile";
	##system($cmd);
	doit($cmd);

	##read results
	my $prev=-99;
	my $nets=0;
	my $highest=-1;
	my $bs;
	my $bc;
	my %votes;
	my $minv=0;
	my $maxv=0;

	open (R,'<',$resfile) or die "Cannot open file $resfile\n";
	open (H, '>', $histfile) or die "Cannot create file $histfile\n";
	print H "Net\tcore1\tstart1\taff1\n";

	while (defined(my $l=<R>)) {
		chomp $l;
		my @e = split(' ',$l);
		next if $#e>6;
		
		if ($l=~ m/\S+\s+(-*\d+)\s+(-*\d+)\s+(\S+)\s+\S+\s+(-*\d+\S+)/)
		{
			my $seq = $1;
			my $start = $2;
			my $core = $3;
			my $aff = $4;
		
			if ($seq <= $prev) {
				print H "$nets\t$bc\t$bs\t$highest\n";
				$votes{$bs} = exists($votes{$bs}) ? $votes{$bs}+1 : 1;
				$minv=$bs if $bs<$minv;
				$maxv=$bs if $bs>$maxv;
				$highest = -1;
				$nets++;
			}
			
			if ($aff>$highest) {
				$highest=$aff;
				$bs = $start;
				$bc = $core;
			}
			
			$prev = $seq;		
        }
	}
	close R;

	#last network
	print H "$nets\t$bc\t$bs\t$highest\n";
	$votes{$bs} = exists($votes{$bs}) ? $votes{$bs}+1 : 1;
	$minv=$bs if $bs<$minv;
	$maxv=$bs if $bs>$maxv;
	$nets++;
	close H;
	
	$cmd = "cat $make_P1_histograms | $R --vanilla --args $histfile $peptide $allele_nn $pdf > $tdir/log.R";
	system($cmd);
	
	###HERE
#	$cmd = "cp $resfile $NetMHCIIpanWWWDIR; cp $histfile $NetMHCIIpanWWWDIR";
#	system($cmd); 
	
	my $link = "$NetMHCIIpanWWWPATH/$$"."_$allele.$peptide.pdf";
	return ($link);		
}

# program usage	
sub usage {
	print "\nUsage: ./NetMHCIIpan-4.0.pl [-h] [args] -f [fastafile/peptidefile]\n";
	print "Command line options:\n\n";
	printf ("%-16s\t%-32s\t%-12s\n",  "PARAMETER", "DEFAULT VALUE", "DESCRIPTION");
	printf ("%-16s\t%-33s\t%-12s\n",  "[-BA]", "0", "Include BA predictions, default is EL only");
	#printf ("%-16s\t%-33s\t%-12s\n",  "[-nout int]", "1", "Result from output neuron of network: 1 for EL, 0 for BA. Should not be used with -BA.");
	printf ("%-16s\t%-33s\t%-12s\n",  "[-context]", "0", "Predict with context encoding");
	printf ("%-16s\t%-33s\t%-12s\n",  "[-rdir dirname]", "$NETMHCIIpan", "Home directory for NetMHCIIpan");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-tdir dirname]", "$TMPDIR/tmp_\$\$", "Temporary directory");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-a name]", "DRB1_0101", "HLA allele");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-choose]", " ", "Choose alpha and beta chains separately");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-cha name]", "", "Alpha chain name");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-chb name]", "", "Beta chain name");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-rankS float]", "2", "Threshold for strong binders (\%Rank)");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-rankW float]", "10", "Threshold for weak binders (\%Rank)");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-filter]", "0", "Toggle filtering of output");	
	printf ("%-16s\t%-32s\t%-12s\n",  "[-rankF float]", "10", "Threshold for filtering output (\%Rank), if -filter option in on");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-hlaseq filename]", " ", "File with full length MHC beta chain sequence (used alone for HLA-DR)");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-hlaseqA filename]", " ", "File with full length MHC alpha chain sequence (used with -hlaseq option)");	
	printf ("%-16s\t%-32s\t%-12s\n",  "[-inptype int]", "0", "Input type [0] FASTA [1] Peptide");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-length int_array]", "15", "Peptide length. Necessary for FASTA input only.");	
	printf ("%-16s\t%-32s\t%-12s\n",  "[-s]", "0", "Sort output on descending affinity");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-u]", "0", "Print unique binding core only");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-ex_offset]", "0", "Exclude offset correction");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-f filename]", " ", "File with the input data");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-xls]", "0", "Save output into xls file");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-xlsfile filename]", "NetMHCIIpan_out.xls", "File name for xls output");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-dirty]", "0", "Dirty mode, leave tmp dir+files");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-w]", "0", "w option for webface");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-list]", "0", "Print the list of possible alleles and exit");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-v]", "0", "Verbose mode" );
	printf ("%-16s\t%-32s\t%-12s\n",  "[-h]", "0", "Print this message and exit");
}
