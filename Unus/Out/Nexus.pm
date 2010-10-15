package Unus::Out::Nexus;
use strict;
use Log::Log4perl qw(:easy);
use Bio::SeqIO;
use File::Basename;

sub new {
	my ($class,$unus) = @_;
	my $self = {
		# Defaults
		'unus'=>$unus,
		'nexnomb'=>$nexnomb,
		'nexnopaup'=>$nexnopaup,
	};
	bless $self, $class;
	# Overwrite defaults where defined
	for ( keys %{$self} ) { $self->{$_}=$unus->{$_} if defined $unus->{$_} }
	#$self->{'genomes'} = $unus->genomes;
	return $self;
}
sub create {
	my($self,@opts) = @_;

}
1;
