package Unus::Blast;
use strict;
use Log::Log4perl qw(:easy);
use Bio::SearchIO;
use File::Basename;
use Digest::MD5 qw(md5_hex);

sub new {
	my ($class,$unus,@opts) = @_;
	my $self = {
		# Defaults
		'unus'=>$unus,
		'tblastx'=>0,
		'evalue'=>10,
		'identity'=>0,
		'similarity'=>0,
		'blastcpus'=>1,
		'blastdir'=>'',
		'blastbins'=>'',
		'blastresults'=>300,
		'tag'=>'',
	};
	my %args = @opts;
	bless $self, $class;
	# Overwrite defaults where defined
	for ( keys %{$self} ) { $self->{$_} = $unus->{$_} if defined $unus->{$_} }
	for ( keys %args ) { $self->{$_} = $args{$_} }
	$self->{'blastdir'} = "/tmp/unus-$$" if -d '/tmp' && (!$self->{'blastdir'});
	mkdir $self->{'blastdir'} if $self->{'blastdir'} && ! -d $self->{'blastdir'};
	LOGDIE "I can't create the BLAST directory '".$self->{'blastdir'}."': $!" if $self->{'blastdir'} && ! -d $self->{'blastdir'};
	$self->{'program'} = $self->{'tblastx'} ? "tblastx" : "blastn";
	$self->{'formatdb'} = ( $self->{'blastbins'} ? $self->{'blastbins'}."/" : "" ) . 'formatdb';
	return $self;
}
sub run {
	my ($self,$query,$db,@opts) = @_;
	$self->{'query'} = $query if $query;
	$self->{'db'} = $db if $db;
	LOGDIE "You must provide a Database to run BLAST" unless defined $self->{'db'};
	LOGDIE "You must provide a Query to run BLAST" unless defined $self->{'query'};
	# Reset parameters
	my %new = @opts;
	my %opt = %{ $self };
	for my $k(keys %new){ my $kn = $k; $kn=~s/^-(.*)/$1/; $opt{$kn}=$new{$k}; }
	# Compile databse if unavailable
	if(($self->{'program'} eq 'blastn' || $self->{'program'} eq 'tblastx' || $self->{'program'} eq 'tblastn') &&
							!(-e $self->{'db'}.".nin" || -e $self->{'db'}.".nal" )){
		system($self->{'formatdb'}, '-i', $self->{'db'}, '-o', 'T', '-p', 'F') == 0
			or LOGDIE "Error running formatdb: $!";
	}elsif(($self->{'program'} eq 'blastp' || $self->{'program'} eq 'blastx') &&
							!(-e $self->{'db'}.".pin" || -e $self->{'db'}.".pal" )){
		system($self->{'formatdb'}, '-i', $self->{'db'}, '-o', 'T', '-p', 'T') == 0
			or LOGDIE "Error running formatdb: $!";
	}
	my $file = "";
	$self->{'report'} = undef;
	$self->{'report_file'} = "";
	if($self->{'blastdir'}){
		$file = $self->{'blastdir'}."/".$opt{'tag'}.
			# ToDo: enhance class detection here
			(-s $self->{'query'} ? md5_hex($self->{'query'}) : $self->{'query'}->display_id).
			"__".basename($self->{'db'}).".blast";
		$self->{'blast_dir'} = $file;
		if(-s $file){
			$self->{'unus'}->msg(6,"Loading BLAST at $file");
			$self->{'report'} = Bio::SearchIO->new(
					-file=>$file,
					-format=>"blastxml");
		}
	}
	unless($self->{'report'}){
		# Load it only if mandatory
		require Bio::Tools::Run::StandAloneBlast;
		my $factory = Bio::Tools::Run::StandAloneBlast->new(
				-database=>$self->{'db'},
				-e=>$opt{'evalue'}, -p=>$self->{'program'},
				-a=>$self->{'blastcpus'},
				-v=>$opt{'blastresults'},
				-b=>$opt{'blastresults'}, -m=>7);
		$factory->o($file) if $file;
		$self->{'unus'}->msg(6,"Running BLAST".($file?" and saving output at $file":""));
		$self->{'report'} = $factory->blastall($self->{'query'});
	}
	if($self->{'report'} && -s $self->{'query'} && $file && $self->{'blastdir'}){
		$self->{'unus'}->msg(6,"Splitting BLAST report saved in '$file'");
		#while(my $result = $self->{'report'}->next_result){
		#	my $f = $self->{'blastdir'}."/".$result->query_accession."__".basename($self->{'db'}).".blast";
		#	my $sO = Bio::SearchIO->new(
		#			-file=>">$f", -output_format=>'blastxml')
		#			or LOGDIE "I can't create the report '$f': $!";
		#	$sO->write_result($result)
		#			or LOGDIE "I can't write the result from query '".$result->query_accession.
		#			"', file '$file': $!";
		#}
		my $xml_header = "";
		my $xml_footer = "  </BlastOutput_iterations>\n</BlastOutput>\n";
		my $xml_oniter = 0;
		my $xml_query  = "";
		my $xml_iter   = "";
		open XML, "<", $file or LOGDIE "I can't read '$file': $!";
		while(<XML>){
			if($xml_oniter){
				if(m{</Iteration>}){
				   print XMLOUT $_;
				   print XMLOUT $xml_footer;
				   close XMLOUT;
				   $xml_query = "";
				   $xml_oniter = 0;
				}elsif($xml_query){
				   print XMLOUT $_;
				}elsif(m{<Iteration_query-def>(.*?)</Iteration_query-def>}){
				   $xml_query = $1;
				   my $Ofile = $self->{'blastdir'}."/".$xml_query."__".basename($self->{'db'}).".blast";
				   open XMLOUT, ">", $Ofile or LOGDIE "I can't create the report '$Ofile': $!";
				   print XMLOUT $xml_header;
				   print XMLOUT $xml_iter;
				   print XMLOUT $_;
				   $xml_iter  = "";
				}else{
				   $xml_iter.= $_;
				}
			}elsif(m{<Iteration>}){
			   $xml_oniter = 1;
			   $xml_iter   = $_;
			}else{
			   $xml_header.= $_;
			}
		}
		close XML;
	}
}
sub hsps {
	my ($self,@opts) = @_;
	my @out=();
	return (wantarray?@out:\@out) unless $self->{'report'};
	
	# Set search parameters
	my %new = @opts;
	my %a = (	-evalue=>$self->{'evalue'},
			-score=>0,
			-bits=>0,
			-length=>0,
			-identity=>$self->{'similarity'},
			-positives=>$self->{'similarity'},
			-onlybesthsp=>0,
			-returnhits=>0);
	for my $k(keys %new){ my $kl = $k; $kl=~s/^([^-].*)/-$1/; $a{$kl}=$new{$k}+0; }
	
	# Create safety file
	if($self->{'report_file'}){copy $self->{'report_file'}, $self->{'report_file'}.".ub" or
			LOGDIE "I can not copy '".$self->{'report_file'}
			."' into '".$self->{'report_file'}.".ub' file: $!"}
	
	# Run
	RESULT:while(my $result = $self->{'report'}->next_result){
		HIT:while(my $hit = $result->next_hit){
			next HIT unless $hit->num_hsps;
			HSP:while(my $hsp = $hit->next_hsp){
				my $ret = $a{'-returnhits'} ? $hit : $hsp;
				push @out, $ret if
					$hsp->evalue<=$a{'-evalue'} &&
					$hsp->score>=$a{'-score'} &&
					$hsp->bits>=$a{'-bits'} &&
					$hsp->length('-query')>=$a{'-length'} &&
					$hsp->frac_identical>=$a{'-identity'} &&
					$hsp->frac_conserved>=$a{'-positives'};
				next HIT if $a{'-onlybesthsp'};
			}
		}
	}
	
	# Remove safety file
	if($self->{'report_file'}){unlink $self->{'report_file'}.".ub" or
			LOGDIE "I can not delete the '".$self->{'report_file'}.".ub' file: $!"}
	return wantarray ? @out : \@out;
}
sub hits {
	my ($self,@opts) = @_;
	return $self->hsps("-onlybesthsp",1,"-returnhits",1,@opts);
	my @out=();
	return (wantarray?@out:\@out) unless $self->{'report'};
	my %new = @opts;
	my %a = (	-evalue=>$self->{'evalue'},
			-score=>0,
			-bits=>0,
			-length=>0,
			-identity=>$self->{'similarity'},
			-positives=>$self->{'similarity'});
	for my $k(keys %new){ my $kl = $k; $kl=~s/^([^-].*)/-$1/; $a{$kl}=$new{$k}+0; }
	if($self->{'report_file'}){copy $self->{'report_file'}, $self->{'report_file'}.".ub" or
			LOGDIE "I can not copy '".$self->{'report_file'}
			."' into '".$self->{'report_file'}.".ub' file: $!"}
	RESULT:while(my $result = $self->{'report'}->next_result){
		HIT:while(my $hit = $result->next_hit){
			my $hsp = $hit->hsp;
			push @out, $hit if
					$hsp->evalue<=$a{'-evalue'} &&
					$hsp->score>=$a{'-score'} &&
					$hsp->bits>=$a{'-bits'} &&
					$hsp->length('-query')>=$a{'-length'} &&
					$hsp->frac_identical>=$a{'-identity'} &&
					$hsp->frac_conserved>=$a{'-positives'};
		}
	}
	if($self->{'report_file'}){unlink $self->{'report_file'}.".ub" or
			LOGDIE "I can not delete the '".$self->{'report_file'}.".ub' file: $!"}
	return wantarray ? @out : \@out;
}
sub best_hit {
	my ($self,@opts) = @_;
	my @hits = $self->hits(@opts);
	return shift @hits;
}
sub best_hsp {
	my ($self,@opts) = @_;
	return $self->best_hit(@opts)->hsp;
}
1;
