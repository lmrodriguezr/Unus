package Unus::Out::Phylip;
use strict;
use Log::Log4perl qw(:easy);
use Bio::AlignIO;

sub new {
	my ($class,$unus) = @_;
	my $self = {
		# Defaults
		'unus'=>$$unus,
		'outloadphylip'=>0,
	};
	bless $self, $class;
	# Overwrite defaults where defined
	for ( keys %{$self} ) { $self->{$_}=$self->{'unus'}->{$_} if defined $self->{'unus'}->{$_} }
	#$self->{'genomes'} = $self->{'unus'}->genomes;
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
	my $outfile = $alndir.".phylip";
	my $coords = $alndir.".coords";
	return $outfile if $self->{'outloadphylip'} && -s $outfile;
	
	$self->{'unus'}->msg(2,"Building the Phylip output");
	my $fasta_i = Bio::AlignIO->new(-file=>$self->{'unus'}->{'finalfastaaln'}, '-format'=>"fasta");
	my $nexus_o = Bio::AlignIO->new(-file=>">$outfile", '-format'=>"phylip");
	while(my $fasta_aln=$fasta_i->next_aln){$nexus_o->write_aln($fasta_aln);}
	return $outfile;
}
1;
