package Unus::Fasta;
use strict;
use Log::Log4perl qw(:easy);
use Bio::SeqIO;
use File::Basename;

sub new {
	my ($class,$unus) = @_;
	my $self = {
		# Defaults
		'unus'=>$unus,
	};
	bless $self, $class;
	# Overwrite defaults where defined
	for ( keys %{$self} ) { $self->{$_}=$unus->{$_} if defined $unus->{$_} }
	$self->{'genomes'} = $unus->genomes;
	return $self;
}
sub fasta_clean {
	my ($self,@in) = @_;
	mkdir $self->{'unus'}->{'basename'}.".in" unless -d $self->{'unus'}->{'basename'}.".in";
	my @out=();
	# First code the names
	my (%coded_taxa,%ct);
	for my $ori (@in){
		my $cod = $self->taxon($ori);
		if($ct{$cod}){
			my $old_cod = $cod;
			my $ori1 = $ct{$cod};
			my $ori2 = $ori;
			my ($cod1,$cod2) = $self->solve_taxa($old_cod, $ori1, $ori2, \%ct);
			$coded_taxa{$ori1} = $cod1;
			$coded_taxa{$ori2} = $cod2;
			$ct{$cod1} = $ori1;
			$ct{$cod2} = $ori2;
			$ct{$cod} = 0;
		}else{
			$coded_taxa{$ori}=$cod;
			$ct{$cod}=$ori;
		}
	}
	$self->{'coded_taxa'} = \%coded_taxa;
	# And then clean the sequences
	for my $genome (@in){
		$self->{'unus'}->msg(3,"Cleaning the genome '$genome' (".$coded_taxa{$genome}.")");
		my $i = 0;
		my $in = Bio::SeqIO->new(-file=>$genome, -format=>'Fasta');
		my $taxon = $coded_taxa{$genome};
		my $out = Bio::SeqIO->new(-file=>">".$self->{'unus'}->{'basename'}.".in/$taxon", -format=>'Fasta');
		while(my $seq = $in->next_seq){
			$seq->display_id("$taxon:".(++$i));
			$seq->desc('');
			$out->write_seq($seq);
		}
		push @out, $self->{'unus'}->{'basename'}.".in/$taxon";
	}
	return wantarray ? @out : \@out;
}
sub taxon {
	my ($self,$dirty,$len,@opts) = @_;
	$len = 4 unless defined $len;
	my $out = File::Basename::basename($dirty);
	return $out if length($out)<=$len;
	$out=~s/\.[^\.]+$//;
	return $out if length($out)<=$len;
	$out=~s/[^A-Za-z0-9]//g;
	return $out if length($out)<=$len;
	$out =~ s/[0-9]//g;
	return $out if length($out)<=$len;
	return substr $out, 0, 4;
}
sub solve_taxa {
	my ($self, $old_cod, $ori1, $ori2, $ctr) = @_;
	$self->{'unus'}->msg(5,"Solving conflict name for '$old_cod'");
	my $pcod1 = $self->taxon($ori1,20);
	my $pcod2 = $self->taxon($ori2,20);
	$self->{'unus'}->msg(6,"Attempting to differentiate '$pcod1' and '$pcod2'");
	if($pcod1 ne $pcod2){
		my $pref_l = 0;
		for my $l (1 .. 8){ $pref_l = $l if substr($pcod1,0,$l) eq substr($pcod2,0,$l) }
		$pcod1 = substr($pcod1,0,3).substr($pcod1,$pref_l,1);
		$pcod2 = substr($pcod2,0,3).substr($pcod2,$pref_l,1);
		$self->{'unus'}->msg(6,"Found different codes: '$pcod1' and '$pcod2'");
	}else{
		$pcod1 = substr($old_cod,0,3)."0";
		$pcod2 = substr($old_cod,0,3)."1";
	}
	$pcod1 = $self->check_taxon($pcod1,$ctr);
	${$ctr}{$pcod1} = 1;
	my @out = ($pcod1, $self->check_taxon($pcod2,$ctr));
	return wantarray ? @out : \@out;
}
sub check_taxon {
	my ($self, $cod, $ctr) = @_;
	my %ct = %{$ctr};
	return $cod unless $ct{$cod};
	for ( 0 .. 9, 'a' .. 'z' ){
		my $recod = substr($cod,0,-1).$_;
		return $recod unless $ct{$recod};
	}
}
sub fasta_count {
	my ($self,@files) = @_;
	my $out = 0;
	for(@files){
		my $seqIO = Bio::SeqIO->new(-file=>$_,-format=>'Fasta');
		$out++ while($seqIO->next_seq);
	}
	return $out;
}
sub get_sequence {
	my ($self,$fasta,$gene_id,@opts) = @_;
	my $seqio = Bio::SeqIO->new(-file=>$fasta,-format=>'Fasta');
	SEQ:while(my $seq = $seqio->next_seq){
		return $seq if($seq->display_id eq $gene_id);
	} # SEQ
	LOGDIE "I can't find the sequence '$gene_id' within '$fasta'.";
	return;
}
1;
