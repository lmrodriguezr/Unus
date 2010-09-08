package Unus::ClustalW;
use strict;
use Bio::Tools::Run::Alignment::Clustalw;
use File::Basename;
use Log::Log4perl qw(:easy);

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
	my $alndir = $self->{'unus'}->{'basename'}.".aln";
	unless(-d $alndir){ mkdir $alndir or LOGDIE "I can not create the '$alndir' directory: $!" }
	my $d=0; $d++ while(<$alndir/*.aln>);
	return $alndir if $d && $self->{'alnload'} && ((!$self->{'unus'}->{'orthgroups'}) || $self->{'unus'}->{'orthgroups'}==$d);
	my $orthdir = $self->{'unus'}->{'orthdir'};
	$orthdir ||= $self->{'unus'}->{'basename'}.".orth";
	unless($self->{'unus'}->{'orthgroups'}){
		$self->{'unus'}->{'orthgroups'} = 0;
		$self->{'unus'}->{'orthgroups'}++ while(my $orthfile = <$orthdir/*.fasta>);
	}
	LOGDIE "Unable to find the aligned sequences and the orthologs groups to align." unless $self->{'unus'}->{'orthgroups'};
	$self->{'unus'}->msg(2,"Aligning orthologs");
	$self->{'unus'}->open_progress('Aligning orthologs', $self->{'unus'}->{'orthgroups'}, 1);
	ORTHFILE:while(my $orthfile = <$orthdir/*.fasta>){
		unless($self->{'unus'}->{'cpus'}==1){$self->{'unus'}->{'pm'}->start and next}
		my $alnfile = $alndir."/".File::Basename::basename($orthfile);
		$alnfile =~ s/\.fasta$/.aln/;
		my $clustalw_factory = Bio::Tools::Run::Alignment::Clustalw->new(-outfile_name=>$alnfile);
		$clustalw_factory->quiet(1);
		$clustalw_factory->align($orthfile);
		$self->{'unus'}->add_progress;
		$self->{'unus'}->{'pm'}->finish unless $self->{'unus'}->{'cpus'}==1;
	} # ORTHFILE
	$self->{'unus'}->{'pm'}->wait_all_children unless $self->{'unus'}->{'cpus'}==1;
	$self->{'unus'}->close_progress;
	return $alndir;
}
1;
