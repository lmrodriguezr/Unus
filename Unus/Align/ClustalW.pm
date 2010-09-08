package Unus::Align::ClustalW;
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
	my $manif = $alndir.".manif";
	unless(-d $alndir){ mkdir $alndir or LOGDIE "I can not create the '$alndir' directory: $!" }
	my $d=0; $d++ while(<$alndir/*.aln>);
	return $alndir if -s $manif && $d && $self->{'alnload'} && ((!$self->{'unus'}->{'orthgroups'}) || $self->{'unus'}->{'orthgroups'}==$d);
	my $orthdir = $self->{'unus'}->{'orthdir'};
	$orthdir ||= $self->{'unus'}->{'basename'}.".orth";
	my $orthmanif = $orthdir.".manif";
	open ORTHMANIF, "<", $orthmanif or LOGDIE "I can not read the '$orthmanif' file.";
	$self->{'unus'}->{'orthgroups'} = 0;
	while(<ORTHMANIF>){
		chomp;
		LOGDIE "I can not find the '$_' file listed in '$orthmanif'." unless -e $_;
		$self->{'unus'}->{'orthgroups'}++
	}
	close ORTHMANIF;
	LOGDIE "Unable to find the aligned sequences and the ortholog groups to align." unless $self->{'unus'}->{'orthgroups'};
	$self->{'unus'}->msg(2,"Aligning orthologs");
	$self->{'unus'}->open_progress('Aligning orthologs', $self->{'unus'}->{'orthgroups'}, 1);
	if(-s $manif) { unlink $manif or LOGDIE "I can not delete the '$manif' file.\n"; }
	open ORTHMANIF, "<", $orthmanif or LOGDIE "I can not read the '$orthmanif' file.";
	ORTHFILE:while(my $orthfile = <ORTHMANIF>){
		unless($self->{'unus'}->{'cpus'}==1){$self->{'unus'}->{'pm'}->start and next}
		chomp $orthfile;
		my $alnfile = $alndir."/".File::Basename::basename($orthfile);
		$alnfile =~ s/\.fasta$/.aln/;
		$self->align($orthfile,$alnfile);
		open MANIF, ">>", $manif or LOGDIE "I can not write in the '$manif' file.";
		print MANIF $alnfile,"\n";
		close MANIF;
		$self->{'unus'}->add_progress;
		$self->{'unus'}->{'pm'}->finish unless $self->{'unus'}->{'cpus'}==1;
	} # ORTHFILE
	$self->{'unus'}->{'pm'}->wait_all_children unless $self->{'unus'}->{'cpus'}==1;
	$self->{'unus'}->close_progress;
	close ORTHMANIF;
	return $alndir;
}
sub align {
	my ($self,$in,$out,@opts) = @_;
	my $factory = Bio::Tools::Run::Alignment::Clustalw->new(-outfile_name=>$out);
	$factory->quiet(1);
	$factory->align($in,@opts);
	return;
}
1;
