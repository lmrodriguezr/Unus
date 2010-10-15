package Unus::Out::Nexus;
use strict;
use Log::Log4perl qw(:easy);
use Bio::AlignIO;
use File::Basename;

sub new {
	my ($class,$unus) = @_;
	my $self = {
		# Defaults
		'unus'=>$$unus,
		'nexnomb'=>0,
		'nexnopaup'=>0,
		'outloadnexus'=>0,
		'outgroup'=>'',
	};
	bless $self, $class;
	# Overwrite defaults where defined
	for ( keys %{$self} ) { $self->{$_}=$self->{'unus'}->{$_} if defined $self->{'unus'}->{$_} }
	$self->{'genomes'} = $self->{'unus'}->genomes;
	return $self;
}
sub create {
	my($self,@opts) = @_;
	# Fasta first
	unless( -s $self->{'unus'}->{'finalfastaaln'} ){
		require Unus::Out::FastaAln;
		my $fasta_aln = Unus::Out::FastaAln->new($self->{'unus'});
		$self->{'unus'}->{'finalfastaaln'} = $fasta_aln->create(@opts);
	}
	my $alndir = $self->{'unus'}->{'alndir'};
	$alndir ||= $self->{'unus'}->{'basename'}.".aln";
	my $outfile = $alndir.".nex";
	my $coords = $alndir.".coords";
	return $outfile if $self->{'outloadnexus'} && -s $outfile;
	
	$self->{'unus'}->msg(2,"Building the Nexus output");
	
	$self->{'unus'}->msg(3,"Writing sequences");
	my $fasta_i = Bio::AlignIO->new(-file=>$self->{'unus'}->{'finalfastaaln'}, '-format'=>"fasta");
	my $nexus_o = Bio::AlignIO->new(-file=>">$outfile", '-format'=>"nexus", -show_symbols=>0);
	while(my $fasta_aln=$fasta_i->next_aln){$nexus_o->write_aln($fasta_aln);}
	unless($self->{'nexnomb'}){
		$self->{'unus'}->msg(3,"Writing MrBayes block");
		open MB, ">>", $outfile or LOGDIE "I can not re-open the '$outfile' file: $!";
		print MB "BEGIN MrBayes;\n  Set AutoClose=yes NoWarn=yes;\n";
		# Partitions
		open COORDS, "<", $coords or LOGDIE "I can not open the '$coords' file: $!";
		my($l,$len) = (0,0);
		while(<COORDS>){
			chomp;
			LOGDIE "I can not parse the line $. of '$coords': '$_'" unless m/^(\d+)\.\.(\d+)$/;
			print MB "  CharSet Gene$.=$1-$2;\n";
			$l = $.; $len = $2+0;
		}
		print MB "  Partition by_gene = $l:Gene1";
		print MB ",Gene$_" for ( 2 .. $l);
		print MB ";\n  CharSet all_aln=1-$len;\n  Partition all = 1:all_aln;\n";
		close COORDS;
		print MB "  Set Partition=by_gene;\n";
		# Taxa
		print MB "  TaxSet ".basename($self->{'genomes'}->[$_])."=".($_+1).";\n" for (0 .. $#{$self->{'genomes'}});
		print MB "  OutGroup ".$self->{'unus'}->{'outgroup'}.";\n" if $self->{'unus'}->{'outgroup'};
		# Run and close
		print MB "  PrSet RatePr=Variable AAModel=Mixed;\n  Unlink Shape=(all) AAModel=(all);\n  MCMCp NGen=1000000 NChains=1 SampleFreq=1000 SaveBrLens=Yes;\n";
		print MB "  Quit;\nEND;\n";
		close MB;
	}
	return $outfile;
}
1;
