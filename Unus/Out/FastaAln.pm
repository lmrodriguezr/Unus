package Unus::Out::FastaAln;
use strict;
use Log::Log4perl qw(:easy);
use Bio::SeqIO;
use Bio::AlignIO;
use File::Basename;

sub new {
	my ($class,$unus) = @_;
	my $self = {
		# Defaults
		'unus'=>$unus,
		'outloadfastaaln'=>0,
	};
	bless $self, $class;
	# Overwrite defaults where defined
	for ( keys %{$self} ) { $self->{$_}=$unus->{$_} if defined $unus->{$_} }
	$self->{'genomes'} = $unus->genomes;
	return $self;
}
sub create {
	my($self,@opts) = @_;
	my $alndir = $self->{'unus'}->{'alndir'};
	$alndir ||= $self->{'unus'}->{'basename'}.".aln";
	my $manif = $alndir.".manif";
	my $coords = $alndir.".coords";
	my $outfile = $alndir.".fasta";
	-d $alndir and -s $manif or LOGDIE "I can not find the alignments";
	return $outfile if $self->{'outloadfastaaln'} && -s $outfile && -s $coords;
	
	open MANIF, "<", $manif or LOGDIE "I can not read the '$manif' file: $!";
	$self->{'unus'}->{'alngroups'} = 0;
	while(<MANIF>){
		chomp;
		LOGDIE "I can not find the '$_' file listed in '$manif'." unless -e $_;
		$self->{'unus'}->{'alngroups'}++
	}
	close MANIF;
	LOGDIE "Unable to find the aligned sequences." unless $self->{'unus'}->{'alngroups'};
	$self->{'unus'}->msg(2,"Building the FastA output");
	my $entireAlnIO = Bio::SeqIO->new(-file => ">$outfile", -format=>"fasta");
	$self->{'unus'}->open_progress('Building the FastA output',$self->{'unus'}->{'alngroups'}*(1+$#{$self->{'genomes'}}), 0);
	GENOME:for my $genome (0 .. $#{$self->{'genomes'}}){
		$self->{'unus'}->msg(3,"Concatening alignment for ".$self->{'genomes'}->[$genome]);
		my $id = basename($self->{'genomes'}->[$genome]);
		my $genomeSeq = Bio::Seq->new();
		$genomeSeq->id($id);
		my $coord = 0;
		open MANIF, "<", $manif or LOGDIE "I can not read the '$manif' file: $!";
		open COORDS, ">", $coords or LOGDIE "I can not write in the '$coords' file: $!";
		open COORDSID, ">", "$coords.$id" or LOGDIE "I can not write in the '$coords.$id' file: $!";
		ALN:while(my $aln = <MANIF>){
			my $alnio = Bio::AlignIO->new(-file=>$aln, -format=>'fasta');
			my $alni = $alnio->next_aln;
			SEQ:foreach my $seq ( $alni->each_seq() ){
				if($seq->display_id() =~ m/^$id:\d+/){
					$genomeSeq->seq($genomeSeq->seq().$seq->seq());
					$coord++;
					my $c = "$coord..".($coord+=$seq->length-1);
					print COORDS "$c\n";
					print COORDSID $seq->display_id().":$c\n";
					last SEQ;
				}
			}
			$self->{'unus'}->add_progress();
		}
		close COORDSID;
		close COORDS;
		close MANIF;
		$entireAlnIO->write_seq($genomeSeq);
	}
	return $outfile;
}
1;
