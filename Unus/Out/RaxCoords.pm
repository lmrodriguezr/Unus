package Unus::Out::RaxCoords;
use strict;
use Log::Log4perl qw(:easy);
use Bio::AlignIO;

sub new {
	my ($class,$unus) = @_;
	my $self = {
		# Defaults
		'unus'=>$$unus,
		'outloadraxcoords'=>0,
		'outgroup'=>'',
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
	my $outfile = $alndir.".raxcoords";
	my $coords = $alndir.".coords";
	return $outfile if $self->{'outloadraxcoords'} && -s $outfile;
	
	$self->{'unus'}->msg(2,"Building the RAxML coordinates");
	
	open COORDS, "<", $coords or LOGDIE "I can not open the '$coords' file: $!";
	open RAXCOORDS, ">", $outfile or LOGDIE "I can not open the '$outfile' file: $!";
	while(<COORDS>){
		chomp;
		LOGDIE "I can not parse the line $. of '$coords': '$_'" unless m/^(\d+)\.\.(\d+)$/;
		print RAXCOORDS "DNA, gene$.=$1-$2\n";
	}
	close COORDS;
	close RAXCOORDS;
	return $outfile;
}
1;
