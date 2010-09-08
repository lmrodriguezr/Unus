package Unus::Orth::Rbh;
use strict;
use Unus::Blast;
use Unus::Fasta;
use Log::Log4perl qw(:easy);

sub new {
	my ($class,$unus) = @_;
	my $self = {
		# Defaults
		'unus'=>$unus,
		'rbhidentity'=>0.7,
		'rbhlength'=>0.7,
		'rbhscorethreshold'=>0.95,
		'rbhxdropoff'=>150,
		'rbhnucpenalty'=>-1,
		'rbhfilter'=>'F',
		'orthloadpairs'=>1,
	};
	bless $self, $class;
	# Overwrite defaults where defined
	for ( keys %{$self} ) { $self->{$_}=$unus->{$_} if defined $unus->{$_} }
	$self->{'genomes'} = $unus->genomes;
	return $self;
}
sub run {
	my($self,@opts) = @_;
	my $orthref_file = $self->build_orthref_file(@opts);
	require Unus::Orth;
	my $orth = Unus::Orth->new($self->{'unus'});
	my $orthtable_file = $orth->filter_orthref_file($orthref_file,@opts);
	return $orthtable_file;
}
sub build_orthref_file {
	my ($self,@opts) = @_;
	my $orthref_file = $self->{'unus'}->{'basename'}.".orthref";
	return $orthref_file if -s $orthref_file && $self->{'orthloadpairs'};
	$self->{'unus'}->msg(3,"Building the orthpairs file");
	# Run
	if(-s $orthref_file){ unlink $orthref_file or LOGDIE "I can't delete the '$orthref_file' file" }
	$self->{'unus'}->open_progress('Building orthology groups', $self->{'unus'}->{'number_of_genes'}, 1);
	my $tmp_fasta = Unus::Fasta->new($self->{'unus'});
	GENESFILE:for my $file (@{$self->{'unus'}->{'genes'}}){
		my $seqIO = Bio::SeqIO->new(-file=>$file,-format=>'Fasta');
		my $blast = Unus::Blast->new($self->{'unus'});
		if($blast->{'blastdir'}){
			for my $genome (@{$self->{'genomes'}}){
				unless($self->{'unus'}->{'cpus'}==1){$self->{'unus'}->{'pm'}->start and next}
				$blast->run($file,$genome);
				$self->{'unus'}->{'pm'}->finish unless $self->{'unus'}->{'cpus'}==1;
			}
			$self->{'unus'}->{'pm'}->wait_all_children unless $self->{'unus'}->{'cpus'}==1;
		}
		GENE:while(my $seq = $seqIO->next_seq){
			$self->{'unus'}->msg(4,"Querying genomes with ".$seq->display_id);
			unless($self->{'unus'}->{'cpus'}==1){$self->{'unus'}->{'pm'}->start and next}
			my $refScore = 0;
			GENOME:for my $genome (0 .. $#{$self->{'genomes'}}){
				$blast->run($seq,$self->{'genomes'}->[$genome]);
				if(my $hit = $blast->best_hit(-identity=>$self->{'rbhidentity'}, -length=>$seq->length*$self->{'rbhlength'})){
					my $back_blast = Unus::Blast->new($self->{'unus'});
					$back_blast->run(
							$tmp_fasta->get_sequence($self->{'genomes'}->[$genome],$hit->accession), $self->{'genomes'}->[0],
							-blastresults=>10,-tag=>'backwards_');
					my @back_hits=$back_blast->hits(
							-identity=>$self->{'rbhidentity'}, -length=>$self->{'rbhlength'},
							-bits=>$self->{'rbhscorethreshold'}*$hit->hsp->bits);
					if(@back_hits && $#back_hits==0){
						if($back_hits[0]->accession eq $seq->display_id){
							open ORTH, ">>", $orthref_file or LOGDIE "I can't write on the '$orthref_file' file: $!";
							print ORTH $seq->display_id."\t$genome\t".$hit->accession."\n";
							close ORTH;
						}
					}
				}
			} # GENOME
			$self->{'unus'}->{'pm'}->finish unless $self->{'unus'}->{'cpus'}==1;
			$self->{'unus'}->add_progress;
		} # GENE
	} # GENESFILE
	$self->{'unus'}->{'pm'}->wait_all_children unless $self->{'unus'}->{'cpus'}==1;
	$self->{'unus'}->close_progress;
	return $orthref_file;
}
1;
