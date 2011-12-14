package Unus::Orth::Bsr;
use strict;
use Unus::Blast;
use Log::Log4perl qw(:easy);

sub new {
	my ($class,$unus) = @_;
	my $self = {
		# Defaults
		'unus'=>$unus,
		'bsrstatic'=>0,
		'bsrtolerance'=>0,
		'bsrfixed'=>'no',
		'orthloadpairs'=>1,
	};
	bless $self, $class;
	# Overwrite defaults where defined
	for ( keys %{$self} ) { $self->{$_}=$unus->{$_} if defined $unus->{$_} }
	$self->{'genomes'} = $unus->genomes;
	$self->{'thresholds'} = [];
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
	# Select the algorithm
	if($self->{'bsrfixed'} ne 'no'){
		$self->{'thresholds'} = $self->{'unus'}->per_genome($self->{'bsrfixed'}+0);
		$self->{'unus'}->msg(3,"BSR static threshold fixed at ".($self->{'bsrfixed'}+0));
	}else{
		require Unus::Orth::BsrAuto;
		my $bsrAuto = Unus::Orth::BsrAuto->new($self->{'unus'});
		$self->{'thresholds'} = $bsrAuto->thresholds;
		if($self->{'bsrstatic'}){
			@{ $self->{'thresholds'} } = $self->{'unus'}->per_genome($bsrAuto->{'mean_thr'});
		}
	}
	LOGDIE "Unexpected error setting the thresholds" unless $self->{'thresholds'};
	$self->{'unus'}->msg(3,"Building the orthpairs file");
	$self->{'unus'}->msg(4,"BSR thresholds: ".join(", ",@{$self->{'thresholds'}}));
	# Run
	if(-s $orthref_file){ unlink $orthref_file or LOGDIE "I can't delete the '$orthref_file' file: $!" }
	$self->{'unus'}->open_progress('Building orthology groups', $self->{'unus'}->{'number_of_genes'}, 1);
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
				my @hits = $blast->hits();
				if(@hits){
					$refScore = $hits[0]->hsp->bits unless $refScore;
					HITS:for my $hit (@hits){
						my $val = $hit->bits/$refScore;
						if($val >= $self->{'thresholds'}->[$genome]){
							$hit->accession or $hit->accession($hit->name);
							$hit->accession or LOGDIE "I can not find the accession nor the name, try description (".$hit->description.")";
							open ORTH, ">>", $orthref_file or LOGDIE "I can't write on the '$orthref_file' file: $!";
							print ORTH $seq->display_id."\t$genome\t".$hit->accession."\n";
							close ORTH;
						}
					} # HSP
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
sub mean_threshold {
	my $self = shift;
	return $self->{'mean_threshold'} if defined $self->{'mean_threshold'};
	return 0 unless $#{$self->{'thresholds'}}>=0;
	my $mean = 0;
	$mean+= ${$self->{'thresholds'}}[$_] for ( 0 .. $#{$self->{'thresholds'}} );
	$self->{'mean_threshold'} = $mean/($#{$self->{'thresholds'}}+1);
	return $self->mean_threshold;
}
1;
