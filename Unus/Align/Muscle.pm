package Unus::Align::Muscle;
use strict;
use Bio::Tools::Run::Alignment::Muscle;
use File::Basename;
use Log::Log4perl qw(:easy);
use Unus::Align;

sub new {
	my ($class,$unus) = @_;
	my $self = {
		# Defaults
		'unus'=>$unus,
		'alnload'=>0,
	};
	bless $self, $class;
	# Overwrite defaults where defined
	for ( keys %{$self} ) { $self->{$_}=$unus->{$_} if defined $unus->{$_} }
	$self->{'genomes'} = $unus->genomes;
	return $self;
}
sub run {
	my($self,@opts) = @_;
	my $align = Unus::Align->new($self->{'unus'});
	$align->start_run(@opts);
	ORTHFILE:while(my $orthfile = <$orthdir/*.fasta>){
		unless($self->{'unus'}->{'cpus'}==1){$self->{'unus'}->{'pm'}->start and next}
		my $alnfile = $alndir."/".File::Basename::basename($orthfile);
		$alnfile =~ s/\.fasta$/.aln/;
		my $muscle_factory = Bio::Tools::Run::Alignment::Muscle->new(-outfile_name=>$alnfile);
		$muscle_factory->quiet(1);
		$muscle_factory->align($orthfile);
		$align->step_run;
		$self->{'unus'}->{'pm'}->finish unless $self->{'unus'}->{'cpus'}==1;
	} # ORTHFILE
	$self->{'unus'}->{'pm'}->wait_all_children unless $self->{'unus'}->{'cpus'}==1;
	$align->finnish_run(@opts);
	return $alndir;
}
1;
