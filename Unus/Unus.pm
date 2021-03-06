package Unus::Unus;

use strict;
use Getopt::Long;
use File::Basename;
use Pod::Usage;
use Log::Log4perl qw(:easy);
use Switch;
use Parallel::ForkManager;

sub new {
	my ($class,@args) = @_;
	my $self = {
		# General
			'verb'=>1,'cpus'=>1,
		# Locations
			'dbprefix'=>'','genomes'=>[],
		# Groups construction
			'orthcriterion'=>'bsr',
	};
	bless $self, $class;
	$self->{'start_time'} = time;
	$self->msg(1,'Unus started');
	# Read the args
	$self->configure(@args);
   	$self->{'pm'} = new Parallel::ForkManager(1);
   	$self->{'pm'}->set_max_procs($self->{'cpus'});
	return $self;
}
sub configure {
	my ($self, @args) = @_;
	@args=@ARGV unless @args;
	my($help,$man)=(0,0);
	$self->{'configpath'} = ".::~/.unus:/etc/unus.d:/opt/unus/conf";
	Getopt::Long::GetOptionsFromArray(\@args,
		# General opts
			'help|?'	=> \$help,
			'man'		=> \$man,
			'verbosity=i'	=> \$self->{'verb'},
			'cpus=i'	=> \$self->{'cpus'},	
			'conf=s'	=> \$self->{'conf'},
			'dummy'		=> sub {  },
		# Locations
			'configpath=s'	=> \$self->{'configpath'},
			'basename=s'	=> \$self->{'basename'},
			'genomes=s{1,}'	=> sub { push @{$self->{'genomes'}},$_[1] if $_[1] },
			'genomesload'	=> \$self->{'genomesload'},
			'dbprefix=s'	=> \$self->{'dbprefix'},
			'blastdir=s'	=> \$self->{'blastdir'},
			'orthdir=s'	=> \$self->{'orthdir'},
			'alndir=s'	=> \$self->{'alndir'},
		# Taxa
			'nameslen=i'	=> \$self->{'nameslen'},
		# Groups construction
			'orthcriterion=s'	=> \$self->{'orthcriterion'},
			'orthloadpairs'		=> \$self->{'orthloadpairs'},
			'orthloadgroups'	=> \$self->{'orthloadgroups'},
			'orthloadseqs'		=> \$self->{'orthloadseqs'},
			'bsrtolerance=f'	=> \$self->{'bsrtolerance'},
			'bsrloadhistogram'	=> \$self->{'bsrloadhistogram'},
			'bsrstatic'		=> \$self->{'bsrstatic'},
			'bsrfixed=f'		=> \$self->{'bsrfixed'},
			'bsrx=i'		=> \$self->{'bsrx'},
			'bsrwins=i'		=> \$self->{'bsrwins'},
			'bsrpolygonratio=i'	=> \$self->{'bsrpolygonratio'},
			'rbhidentity=f'		=> \$self->{'rbhidentity'},
			'rbhlength=f'		=> \$self->{'rbhlength'},
			'rbhscorethreshold=f'	=> \$self->{'rbhscorethreshold'},
			'rbhxdropoff=i'		=> \$self->{'rbhxdropoff'},
			'rbhnucpenalty=i'	=> \$self->{'rbhnucpenalty'},
			'rbhfilter=s'		=> \$self->{'rbhfilter'},
			'rbhnoextratables'	=> \$self->{'rbhnoextratables'},
		# BLAST
			'tblastx'		=> \$self->{'tblastx'},
			'evalue=f'		=> \$self->{'evalue'},
			'identity=f'		=> \$self->{'identity'},
			'similarity=f'		=> \$self->{'similarity'},
			'blastbins=s'		=> \$self->{'blastbins'},
			'blastresults=i'	=> \$self->{'blastresults'},
			'blastcpus=i'		=> \$self->{'blastcpus'},
			'blastoutformat=s'	=> \$self->{'blastoutformat'},
		# Alignments
			'alnmethod=s'	=> \$self->{'alnmethod'},
			'alnload'	=> \$self->{'alnload'},
		# Tests
			'recombinationtest=s'	=> \$self->{'recombinationtest'},
			'recombinationtestload'	=> \$self->{'recombinationtestload'},
			'phitestbin=s'		=> \$self->{'phitestbin'},
			'phitestsignificance=f'	=> \$self->{'phitestsignificance'},
			'modeltest=s'		=> \$self->{'modeltest'},
		# Output files
			'nonexus'		=> \$self->{'nonexus'},
			'outloadfastaaln'	=> \$self->{'outloadfastaaln'},
			'outloadnexus'		=> \$self->{'outloadnexus'},
			'nexnomb'		=> \$self->{'nexnomb'},
			'nexnopaup'		=> \$self->{'nexnopaup'},
			'nophylip'		=> \$self->{'nophylip'},
			'outloadphylip'		=> \$self->{'outloadphylip'},
			'noraxcoords'		=> \$self->{'noraxcoords'},
			'outloadraxcoords'	=> \$self->{'outloadraxcoords'},
	) or Pod::Usage::pod2usage(2);
	Pod::Usage::pod2usage(1) if $help;
	Pod::Usage::pod2usage(-exitstatus => 0, -verbose => 2) if $man;
	if($self->{'conf'}){
		my @arr = $self->load_conf($self->{'conf'});
		$self->configure(@arr);
	}
	LOGDIE "You must provide a value for -basename" unless $self->{'basename'};
	LOGDIE "You must provide a value for -genomes" unless $self->{'genomes'};
	$self->{'genes'} = [$self->{'genomes'}->[0]];
}
sub load_conf {
	my($self,$conf,@opts) = @_;
	$conf = $self->{'conf'} unless defined $conf;
	PATH:for my $dir ( split ':', $self->{'configpath'} ){
		$dir =~ s{/+$}{};
		my $c = $dir . "/" . $conf;
		if ( -s $c || -s $c . ".conf") {
			$conf = $dir . "/" . $conf . (-s $c?"":".conf");
			last PATH;
		}
	}
	$self->msg(3,"Reading configuration file '".$conf."'");
	LOGDIE "I can't read the configuration file '".$conf."': $!" unless -r $conf;
	my @arr = ('-dummy');
	open CNF, "<", $conf or LOGDIE "I can't open the file '".$conf."': $!";
	my @cnf = <CNF>;
	close CNF;
	for(@cnf){
		chomp;
		s/\s*#.*$//;
		next if m/^\s*$/;
		if(m/^load\s+(.*)/){
			#$self->{'conf'} = $1;
			push @arr, $self->load_conf($1);
		}elsif(m/^([^=]+)=(.*)/){
			next if $1 eq 'conf';
			push @arr, ("-$1",$2);
			$self->msg(5,"Setting '$1' from file.");
		}else{
			LOGDIE "I can't parse the line $. of '".$conf."': $_";
		}
	}
	$self->{'conf'}="";
	return wantarray ? @arr : \@arr;
}
sub run {
	my ($self,@opts) = @_;
	## [0] Clean data
	$self->clean_data(@opts);
	# [1] Calculate the orthology groups
	$self->calculate_orthologs(@opts);
	# [2] Align sequences
	$self->align(@opts);
	# [3] Tests
	$self->tests(@opts);
	# [4] Phylo files
	$self->phylo_files(@opts);
	# [5] Inference
}
sub clean_data {
	my($self,@opts) = @_;
	require Unus::Fasta;
	my $fasta = Unus::Fasta->new(\$self);
	$self->{'ori_genomes'} = $self->{'genomes'};
	$self->{'genomes'} = $fasta->fasta_clean(@{$self->{'genomes'}});
	$self->{'genes'} = [$self->{'genomes'}->[0]];
	$self->{'number_of_genes'} = $fasta->fasta_count(@{$self->{'genes'}});
}
sub calculate_orthologs {
	my ($self,@opts) = @_;
	my $groups;
	switch ( lc $self->{'orthcriterion'} ){
		case 'bsr' {
			$self->msg(2,"Selected criterion: Bit Score Ratio");
			require Unus::Orth::Bsr;
			my $bsr = Unus::Orth::Bsr->new($self);
			$groups = $bsr->run;
		}
		case 'rbh' {
			$self->msg(2,"Selected criterion: Reciprocal Best Hit");
			require Unus::Orth::Rbh;
			my $rbh = Unus::Orth::Rbh->new($self);
			$groups = $rbh->run;
		}
		case 'orthomcl' {
			$self->msg(2,"Selected criterion: OrthoMCL [Li, Stoeckert & Roos 2003 Genome Res 13(9):2178-89]");
			Pod::Usage::pod2usage({-exitval=>1, -msg=>"Value for -orthcriterion still unimplemented"});
		}
		case /none|noorth/ {
			$self->msg(3,"No orthology search, reading further input");
		}
		else {
			Pod::Usage::pod2usage({-exitval=>1, -msg=>"Bad value for -orthcriterion: ".$self->{'orthcriterion'}});
		}
	}
	if($groups){
		require Unus::Orth;
		my $orth = Unus::Orth->new($self);
		$self->{'orthdir'} = $orth->groups2fasta($groups);
	}
}
sub align {
	my ($self,@opts) = @_;
	-d $self->{'orthdir'} or
		LOGDIE "I can not find the '".$self->{'orthdir'}."' directory, use the -orthdir parameter or build the groups with -orthcriterion.";
	-s $self->{'orthdir'}.".manif" or
		LOGDIE "I can not find the '".$self->{'orthdir'}.".manif' file, use the -orthdir parameter or build the groups with -orthcriterion.";
	switch ( lc $self->{'alnmethod'} ){
		case 'muscle' {
			$self->msg(3,"Selected alignment method: Muscle [Edgar 2004 NAR 32(5):1792-7]");
			require Unus::Align::Muscle;
			my $muscle = Unus::Align::Muscle->new($self);
			$self->{'alndir'} = $muscle->run(@opts);
		}
		case 'clustalw' {
			$self->msg(3,"Selected alignment method: ClustalW [Larkin et al 2007 Bioinformatics 23(21):2947-8]");
			require Unus::Align::ClustalW;
			my $clustalw = Unus::Align::ClustalW->new($self);
			$self->{'alndir'} = $clustalw->run(@opts);
		}
		case /none|noaln/ {
			$self->msg(3,"No alignment execution, reading further input");
		}
		else {
			Pod::Usage::pod2usage({-exitval=>1, -msg=>"Bad value for -alnmethod: ".$self->{'alnmethod'}});
		}
	}
}
sub tests {
	my ($self,@opts) = @_;
	-d $self->{'alndir'} or
		LOGDIE "I can not find the '".$self->{'alndir'}."' directory, use the -alndir parameter or build the alignments with -alnmethod.";
	-s $self->{'alndir'}.".manif" or
		LOGDIE "I can not find the '".$self->{'alndir'}.".manif' file, use the -alndir parameter or build the alignments with -alnmethod.";
	switch ( $self->{'recombinationtest'} ) {
		case 'phi' {
			$self->msg(3,"Selected test: PHI-Test for recombination detection [Bruen, Philippe & Bryant 2006 Genetics 172(4):2665-81]");
			require Unus::Test::Phi;
			my $phitest = Unus::Test::Phi->new($self);
			#$self->{'alndir'} =
			$phitest->run(@opts);
		}
		case /|none/ {
			$self->msg(3,"Any recombination test selected, ignoring recombination.");
		}
		else {
			Pod::Usage::pod2usage({-exitval=>1, -msg=>"Bad value for -recombinationtest: ".$self->{'recombinationtest'}});
		}
	}
	switch ( $self->{'modeltest'} ) {
		case 'modeltest' {
			$self->msg(2,"This -modeltest value is still unimplemented, ignoring...");
		}
		case 'mrmodeltest' {
			$self->msg(2,"This -modeltest value is still unimplemented, ignoring...");
		}
		case 'jmodeltest' {
			$self->msg(2,"This -modeltest value is still unimplemented, ignoring...");
		}
		case /|none/ {
			$self->msg(2,"Any substitution model test selected, ignoring substitution models.");
		}
		else {
			Pod::Usage::pod2usage({-exitval=>1, -msg=>"Bad value for -modeltest: ".$self->{'modeltest'}});
		}
	}
}
sub phylo_files {
	my ($self,@opts) = @_;
	-d $self->{'alndir'} or
		LOGDIE "I can not find the '".$self->{'alndir'}."' directory, use the -alndir parameter or build the alignments with -alnmethod.";
	-s $self->{'alndir'}.".manif" or
		LOGDIE "I can not find the '".$self->{'alndir'}.".manif' file, use the -alndir parameter or build the alignments with -alnmethod.";
	unless ( $self->{'nonexus'} ){
		# Nexus
		require Unus::Out::Nexus;
		my $nexus = Unus::Out::Nexus->new(\$self);
		$self->{'finalnexus'} = $nexus->create(@opts);
	}
	unless ( $self->{'nophylip'} ){
		# Phylip
		require Unus::Out::Phylip;
		my $phylip = Unus::Out::Phylip->new(\$self);
		$self->{'finalphylip'} = $phylip->create(@opts);
	}
	unless ( $self->{'noraxcoords'} ){
		# RAxML coords file
		require Unus::Out::RaxCoords;
		my $raxcoords = Unus::Out::RaxCoords->new(\$self);
		$self->{'finalraxcoords'} = $raxcoords->create(@opts);
	}
}
## S Y S T E M - W I D E   F U N C T I O N S ##
sub genomes {
	my ($self,@opts) = @_;
	my %args = @opts;
	my @out = @{$self->{'genomes'}};
	if(!$args{'-suffix'} && !$args{'-noprefix'}){
		for(0 .. $#out){
			$out[$_] = $self->{'dbprefix'}.$out[$_];
		}
	}
	if($args{'-onlyname'}){
		for(0 .. $#out){
			$out[$_] = File::Basename::basename($out[$_]);
		}
	}
	return wantarray ? @out : \@out;
}
sub per_genome {
	my ($self,$val) = @_;
	my @out = ();
	push @out, $val for (0 .. $#{$self->{'genomes'}});
	return wantarray ? @out : \@out;
}
sub msg {
	my ($self,$verb,$msg) = @_;
	print "".($self->{'verb'}>1?"[".sec2hr(time-$self->{'start_time'})."]":"").(" "x$verb)." $msg\n" if $verb <= $self->{'verb'};
	#print #($self->{'verb'}>1?"[".sec2hr(time - $self->{'start_time'})."] ":"").
	#	("·"x$verb)." ".$msg."\n" if $verb <= $self->{'verb'};
}
sub open_progress {
	my($self,$task,$N,$parallel,@opts) = @_;
	$self->{'progress'} = 0;
	$self->{'progress_task'} = $task;
	$self->{'progress_size'} = $N;
	$self->{'progress_start'} = time;
	$self->{'progress_task'}.= " (Impossible to calculate progress in parallel mode - ".$self->{'cpus'}." cpus)" if $parallel && $self->{'cpus'}>1;
	$self->{'progress_parallel'} = ($parallel && $self->{'cpus'}>1);
	return 0 unless $self->{'verb'} == 1;
	$| = 1;
	print " ".$self->{'progress_task'}."...";
}
sub add_progress {
	my($self,@opts) = @_;
	$self->{'progress'}++;
	return 0 unless $self->{'verb'} == 1;
	print "\r ".$self->{'progress_task'};
	return 0 if $self->{'progress_parallel'};
	if($self->{'progress_size'}){
		my $box = 100; # -> it would be nice to guess the terminal size here!, try Term::Screen
		$box-= length($self->{'progress_task'}) + 10 + 15;
		my $m = int(1000*($self->{'progress'})/$self->{'progress_size'});
		my $p = int($box*$m/1000);
		print " [".("="x$p).">".(" "x($box-$p))."] ".sprintf("%.1f",$m/10)."%  Left: ".
			sec2hr((time-$self->{'progress_start'})*($self->{'progress_size'}-$self->{'progress'})/($self->{'progress'}))."     ";
	}else{
		print ": ".$self->{'progress'}.".  Elapsed: ".(time - $self->{'progress_start'})."s     ";
	}
}
sub close_progress {
	my($self,@opts) = @_;
	if($self->{'verb'} == 1){
		my $o = " ".$self->{'progress_task'}."... done.  Elapsed time: ".sec2hr(time-$self->{'progress_start'});
		my $box = 100;
		print "\r$o".(" "x($box-length($o)))."      \n";
		$| = 0;
	}
	$self->{'progress'} = 0;
	$self->{'progress_task'} = '';
	$self->{'progress_size'} = 0;
	$self->{'progress_start'} = 0;
	$self->{'progress_parallel'} = 0;
}
sub sec2min {
	my $sec = shift;
	my $min_dec = $sec/60;
	my $min = int($min_dec);
	my $sec_dec = $min_dec-$min;
	$sec = int($sec_dec*60);
	return wantarray ? ($min,$sec) : "$min:".("0"x(2-length("$sec")))."$sec";
}
sub sec2hr {
	my $sec = shift;
	my ($min, $sec) = sec2min($sec);
	my ($hr, $min) = sec2min($min);
	return wantarray ? ($hr,$min,$sec) : "$hr:".("0"x(2-length "$min"))."$min:".("0"x(2-length "$sec"))."$sec";
}

1;

