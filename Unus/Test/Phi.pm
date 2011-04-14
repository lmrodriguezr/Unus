package Unus::Test::Phi;
use strict;
use Log::Log4perl qw(:easy);
use File::Copy;
use File::chdir;
use File::Basename;
use Cwd 'abs_path';

sub new {
	my ($class,$unus) = @_;
	my $self = {
		# Defaults
		'unus'=>$unus,
		'phitestbin'=>'Phi',
		'phitestsignificance'=>0.05,
		'recombinationtestload'=>0,
	};
	bless $self, $class;
	# Overwrite defaults where defined
	for ( keys %{$self} ) { $self->{$_}=$unus->{$_} if defined $unus->{$_} }
	$self->{'genomes'} = $unus->genomes;
	return $self;
}
sub run {
	my($self,@opts) = @_;
	my $alndir = $self->{'unus'}->{'alndir'};
	$alndir ||= $self->{'unus'}->{'basename'}.".aln";
	-d $alndir and -s "$alndir.manif" or LOGDIE "I can not find the alignments";
	if(-s "$alndir.manif.prephi") {
		if($self->{'recombinationtestload'} && -s "$alndir.manif.postphi"){
			copy "$alndir.manif.postphi", "$alndir.manif" or LOGDIE "I can not copy '$alndir.manif.postphi' into '$alndir.manif': $!";
			return "$alndir.manif" 
		}else{
			copy "$alndir.manif.prephi", "$alndir.manif" or LOGDIE "I can not copy '$alndir.manif.prephi' into '$alndir.manif': $!";
		}
	}else{
		copy "$alndir.manif", "$alndir.manif.prephi" or LOGDIE "I can not copy '$alndir.manif' into '$alndir.manif.prephi': $!";
	}
	return $self->test_recombination("$alndir.manif.prephi",$alndir,@opts);
}
sub test_recombination {
	my ($self,$oldmanif,$alndir,@opts) = @_;
	my $newmanif = "$alndir.manif";
	if(-s $newmanif){ unlink $newmanif or LOGDIE "I can not delete the '$newmanif' file: $!" }
	my $prealns = 0;
	open OLDMANIF, "<", $oldmanif or LOGDIE "I can not read the '$oldmanif' file: $!";
	$prealns++ while(<OLDMANIF>);
	close OLDMANIF;
	open OLDMANIF, "<", $oldmanif or LOGDIE "I can not read the '$oldmanif' file: $!";
	$self->{'unus'}->msg(3,"Running PHI tests");
	$self->{'unus'}->open_progress('Running PHI tests', $prealns, 1);
	OLD:while(my $alnrel=<OLDMANIF>){
		chomp $alnrel;
		my $aln = abs_path $alnrel;
		my $phi = $self->{'phitestbin'}=~m/\// ? abs_path $self->{'phitestbin'} : $self->{'phitestbin'};
		mkdir "$aln.phi" or LOGDIE "I can not create the '$aln.phi' directory: $!";
		my $olddir = Cwd::abs_path;
		$self->{'unus'}->msg(4,"Testing recombination in $aln");
		unless($self->{'unus'}->{'cpus'}==1){$self->{'unus'}->{'pm'}->start and next}
		
		chdir "$aln.phi";
		system "$phi -f '$aln' &>$aln.phi/Phi.out";
		LOGDIE "I can not execute the PHI test: $phi:\n$!" unless $?==0;
		-s "Phi.log" or LOGDIE "PHI test ($phi -f '$aln' &>$aln.phi/Phi.out) returned an empty output, please check '$aln.phi/Phi.out' to inspect the reasons.";
		copy "Phi.inf.sites", "$aln.phi.inf.fasta";
		copy "Phi.log","$aln.phi.log";
		unlink "Phi.inf.sites", "Phi.inf.list", "Phi.log","Phi.out";
		chdir $olddir;
		
		rmdir "$aln.phi" or LOGDIE "I can not delete the '$aln.phi' directory: $!";
		open PHI, "<", "$aln.phi.log" or LOGDIE "I can not read the PHI report '$aln.phi.log': $!";
		while(<PHI>){
			if(m/^PHI \(Normal\): *([\d\.eE-]+)/){
				if($1 ne "--" && $1+0<=$self->{'phitestsignificance'}){
					$self->{'unus'}->msg(5,"Discarded alignment for recombination suspect (p=".($1+0)."): ".basename($aln));
				}else{
					open NEWMANIF, ">>", $newmanif or LOGDIE "I can not write in the '$newmanif' file: $!";
					print NEWMANIF $alnrel,"\n";
					close NEWMANIF;
				}
			}
		}
		close PHI;
		$self->{'unus'}->{'pm'}->finish unless $self->{'unus'}->{'cpus'}==1;
		$self->{'unus'}->add_progress;
	}
	$self->{'unus'}->{'pm'}->wait_all_children unless $self->{'unus'}->{'cpus'}==1;
	$self->{'unus'}->close_progress;
	close OLDMANIF;
	copy $newmanif, "$alndir.manif.postphi" or LOGDIE "I can not copy '$newmanif' into '$alndir.manif.postphi': $!";
	
	return $newmanif;
}
1;
