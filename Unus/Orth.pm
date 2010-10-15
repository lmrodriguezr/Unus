package Unus::Orth;
use strict;
use Unus::Fasta;
use Bio::SeqIO;
use Log::Log4perl qw(:easy);

sub new {
	my ($class,$unus) = @_;
	my $self = {
		# Defaults
		'unus'=>$unus,
		'orthloadgroups'=>1,
		'orthloadseqs'=>0,
	};
	bless $self, $class;
	# Overwrite defaults where defined
	for ( keys %{$self} ) { $self->{$_}=$unus->{$_} if defined $unus->{$_} }
	$self->{'genomes'} = $unus->genomes;
	return $self;
}
sub filter_orthref_file {
	my($self,$orthref_file,@opts) = @_;
	my $orthtable_file = $self->{'unus'}->{'basename'}.".orthtable";
	return $orthtable_file if -s $orthtable_file && $self->{'orthloadgroups'};
	$self->{'unus'}->msg(3,"Filtering the orthologous pairs and building groups");
	$self->{'unus'}->open_progress('Building filtered orthology groups', $self->{'unus'}->{'number_of_genes'}, 1);
	if(-s $orthtable_file){ unlink $orthtable_file or LOGDIE "I can't delete the '$orthtable_file' file: $!" }
	GENESFILE:for my $file (@{$self->{'unus'}->{'genes'}}){
		my $seqIO = Bio::SeqIO->new(-file=>$file,-format=>'Fasta');
		GENE:while(my $seq = $seqIO->next_seq){
			unless($self->{'unus'}->{'cpus'}==1){$self->{'unus'}->{'pm'}->start and next}
			my $id = $seq->display_id;
			my @group = ();
			my $ok = 1;
			open ORTH, "<", $orthref_file or LOGDIE "I can't read the '$orthref_file' file: $!";
			# Fetch pairs of the reference gene
			PAIR:while(<ORTH>){
				next PAIR unless m/^$id\t/;
				chomp;
				my ($s,$genome,$o) = split /\t/;
				if($group[$genome]){
					$self->{'unus'}->msg(7,"Discarded '$id', duplicated (at least '".$group[$genome]."' and '$o') in the genome $genome");
					$ok = 0;
					last PAIR;
				}
				$group[$genome] = $o;
			} # PAIR
			close ORTH;
			# Check possible holes
			GENOME:for my $genome (0 .. $#{$self->{'genomes'}}){
				unless($group[$genome]){
					$self->{'unus'}->msg(7,"Discarded '$id', absent in the genome $genome");
					$ok = 0;
					last GENOME;
				}
			} # GENOME
			# Save results if OK
			if($ok){
				open TABLE, ">>", $orthtable_file or LOGDIE "I can't write on the the '$orthtable_file' file: $!";
				print TABLE join("\t",@group)."\n";
				close TABLE;
			}
			$self->{'unus'}->{'pm'}->finish unless $self->{'unus'}->{'cpus'}==1;
			$self->{'unus'}->add_progress;
		} #GENE
	} # GENESFILE
	$self->{'unus'}->{'pm'}->wait_all_children unless $self->{'unus'}->{'cpus'}==1;
	$self->{'unus'}->close_progress;
	$self->{'unus'}->{'orthgroups'} = 0;
	open TABLE, "<", $orthtable_file or LOGDIE "I can't read the the '$orthtable_file' file: $!";
	$self->{'unus'}->{'orthgroups'}++ while(<TABLE>);
	close TABLE;
	return $orthtable_file;
}
sub groups2fasta {
	my($self,$orthtable_file,@opts) = @_;
	my $orthdir = $self->{'unus'}->{'basename'}.".orth";
	my $manif = $orthdir.".manif";
	unless(-d $orthdir){ mkdir $orthdir or LOGDIE "I can not create the '$orthdir' directory: $!" }
	my $d=0; $d++ while(<$orthdir/*.fasta>);
	return $orthdir if -s $manif && $d && $self->{'orthloadseqs'} && ((!$self->{'unus'}->{'orthgroups'}) || $self->{'unus'}->{'orthgroups'}==$d);
	$self->{'unus'}->msg(3,"Extracting orthologous genes and compiling fasta files");
	open TABLE, "<", $orthtable_file or LOGDIE "I can't read the '$orthtable_file' file: $!";
	my $lines=0;
	$lines++ while(<TABLE>);
	close TABLE;
	$self->{'unus'}->msg(4,"$lines groups found");
	$self->{'unus'}->open_progress('Extracting genes', $lines, 1);
	my $fasta = Unus::Fasta->new($self->{'unus'});
	open TABLE, "<", $orthtable_file or LOGDIE "I can't read the '$orthtable_file' file: $!";
	if(-s $manif) { unlink $manif or LOGDIE "I can not delete the '$manif' file."; }
	TABLE:while(my $group=<TABLE>){
		unless($self->{'unus'}->{'cpus'}==1){$self->{'unus'}->{'pm'}->start and next}
		chomp $group;
		my @row = split /\t/, $group;
		my $id = $row[0];
		$self->{'unus'}->msg(6,"Building group for '$id'");
		$id=~s/[^A-Za-z0-9:]/_/g;
		my $out = Bio::SeqIO->new(-file=>">".$orthdir."/$id.fasta");
		GENOME:for my $genome (0 .. $#row){
			$out->write_seq($fasta->get_sequence($self->{'genomes'}->[$genome], $row[$genome]));
		}
		open MANIF, ">>", $manif or LOGDIE "I can not write in the '$manif' file.";
		print MANIF $orthdir."/$id.fasta\n";
		close MANIF;
		$self->{'unus'}->add_progress;
		$self->{'unus'}->{'pm'}->finish unless $self->{'unus'}->{'cpus'}==1;
	}
	$self->{'unus'}->{'pm'}->wait_all_children unless $self->{'unus'}->{'cpus'}==1;
	close TABLE;
	$self->{'unus'}->close_progress;
	return $orthdir;
}
1;
