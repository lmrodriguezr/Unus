package Unus::Orth::BsrAuto;
use strict;
use Bio::SeqIO;
use Unus::Blast;
use Log::Log4perl qw(:easy);

sub new {
	my ($class,$unus,@opts) = @_;
	my $self = {
		# Defaults
		'unus'=>$unus,
		'bsrtolerance'=>0,
		'bsrx'=>100,
		'bsrwins'=>20,
		'bsrloadhistogram'=>0,
		'bsrpolygonratio'=>10,
	};
	bless $self, $class;
	# Overwrite defaults where defined
	for ( keys %{$self} ) { $self->{$_} = $unus->{$_} if defined $unus->{$_} }
	$self->{'genomes'} = $unus->genomes;
	$self->{'thresholds'} = [];
	return $self;
}
sub thresholds {
	my ($self,@opts) = @_;
	return $self->{'thresholds'} if $#{$self->{'thresholds'}}==$#{$self->{'genomes'}};
	# Build the histogram:
	$self->{'unus'}->msg(2,"Building the BSR histogram");
	my $data_file = $self->extract_values(@opts);
	$self->build_histogram($data_file,@opts);
	return $self->{'thresholds'};
}
sub extract_values {
	my ($self,@opts) = @_;
	my $hist_file = $self->{'unus'}->{'basename'}.".histogram";
	return $hist_file if -s $hist_file && $self->{'bsrloadhistogram'};
	$self->{'unus'}->msg(3,"Extracting values to build the BSR histogram");
	$self->{'unus'}->open_progress('Building the BSR histogram', $self->{'unus'}->{'number_of_genes'}, 1);
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
				my @hsps = $blast->hsps();
				$self->{'unus'}->msg(5,"Found ".($#hsps+1)." HSPs between '".$seq->display_id."' and genome '$genome'");
				if(@hsps){
					$refScore = $hsps[0]->bits unless $refScore;
					HSP:for(@hsps){
						my $val = int($_->bits*$self->{'bsrx'}/($refScore))/$self->{'bsrx'};
						open HISTOGRAM, ">>", $hist_file or LOGDIE "I can't write in the histogram file '$hist_file': $!";
						print HISTOGRAM "$genome\t$val\n";
						close HISTOGRAM;
					}
				}else{
					open HISTOGRAM, ">>", $hist_file or LOGDIE "I can't write in the histogram file '$hist_file': $!";
					print HISTOGRAM "$genome\t0\n";
					close HISTOGRAM;
				}
			} # GENOME
			$self->{'unus'}->{'pm'}->finish unless $self->{'unus'}->{'cpus'}==1;
			$self->{'unus'}->add_progress;
		} # GENE
	} # GENESFILE
	$self->{'unus'}->{'pm'}->wait_all_children unless $self->{'unus'}->{'cpus'}==1;
	$self->{'unus'}->close_progress;
	return $hist_file;
}
sub build_histogram {
	my ($self,$hist_file,@opts) = @_;
	$self->{'bsrdata'} = ();
	$self->{'bsrhist'} = ();
	$self->{'unus'}->msg(3,"Reading histogram and setting thresholds");
	open HISTOGRAM, "<", $hist_file or LOGDIE "I can't read the histogram file '$hist_file': $!";
	while(<HISTOGRAM>){
		chomp;
		m/^(\d+)\t([\d\.]+)$/ or LOGDIE "I can't parse the $. line from '$hist_file': $_";
		$self->{'bsrdata'}->[$1] = [] unless defined $self->{'bsrdata'}->[$1];
		my $bsratio = $2>1 && $2<($self->{'unus'}->{'tblastx'}?3:1.5) ? 1 : $2+0; # <- This is to correct a small BLAST bug
		push @{$self->{'bsrdata'}->[$1]}, $bsratio;
		unless(defined $self->{'bsrhist'}->[$1]){ $self->{'bsrhist'}->[$1]->[$_] = 0 for (0 .. $self->{'bsrx'}) }
		LOGDIE "BSR '$bsratio' out of range '0-".$self->{'bsrx'}."'" unless defined $self->{'bsrhist'}->[$1]->[$bsratio*$self->{'bsrx'}];
		$self->{'bsrhist'}->[$1]->[$bsratio*$self->{'bsrx'}]++;
	}
	close HISTOGRAM;
	$self->{'unus'}->msg(3,"Analyzing the BSR distribution in the search for thresholds");
	$self->{'windows'} = [];
	$self->{'polygons'} = [];
	$self->{'known_topology'} = [];
	$self->{'win_thrs'} = [];
	$self->{'mean_thr'} = 0;
	open BSRGRAPH, ">", $self->{'unus'}->{'basename'}.".bsrgraph.R" or LOGDIE "I can't open the bsrgraph R script file '".
			$self->{'unus'}->{'basename'}.".bsrgraph.R': $!";
	print BSRGRAPH "idx <- seq(from=0,to=1,by=".(1/$self->{'bsrx'}).");\n";
	print BSRGRAPH "idx_win <- seq(from=0,to=1,by=".(1/$self->{'bsrwins'}).");\n";
	GENOME: for my $genome (0 .. $#{ $self->{'unus'}->{'genomes'} }){
		$self->{'windows'}->[$genome] = [];
		for my $i(0 .. $self->{'bsrwins'}){
			$self->{'windows'}->[$genome]->[$i] = 1;
			my $a = ($i-1)/$self->{'bsrwins'};
			my $b = ( $i )/$self->{'bsrwins'};
			for my $j (@{ $self->{'bsrdata'}->[$genome] }){
				$self->{'windows'}->[$genome]->[$i]++ if $j>$a and $j<=$b;
			}
		}
		print BSRGRAPH "genome$genome <- c(".join(",",@{$self->{'bsrhist'}->[$genome]}).");\n";
		print BSRGRAPH "win_genome$genome <- c(".join(",",@{$self->{'windows'}->[$genome]}).");\n";
		my $polygonsides = 4;
		my @polygon = (1,$self->{'bsrwins'});
		my @win_mean_log = ();
		POLYSIDE: while($#polygon+1<$polygonsides){
			my @best_dot = (0,0); # (index,path-len)
			POLYDOT: for my $dot (1 .. $self->{'bsrwins'}){
				$win_mean_log[$dot] = -1*$self->{'bsrpolygonration'} unless $self->{'windows'}->[$genome]->[$dot]>0;
				$win_mean_log[$dot] = log($self->{'windows'}->[$genome]->[$dot])*$self->{'bsrpolygonratio'}
						unless defined $win_mean_log[$dot];
				my $path_len = 1;
				POLYEDGE: for my $edge (@polygon){
					$path_len*= sqrt( ($dot-$edge)**2 + ($win_mean_log[$dot]-$win_mean_log[$edge])**2 ); # Euclidean distance
				}
				@best_dot = ($dot, $path_len) if $path_len>$best_dot[1];
			}
			push(@polygon, $best_dot[0]);
		}
		if ($polygon[2]>$polygon[3]){
			my $tmp = $polygon[2];
			$polygon[2]=$polygon[3];
			$polygon[3]=$tmp;
		}
		# Detect the lowest window within the range
		my $minwin = $polygon[2];
		for my $win($polygon[2] .. $polygon[3]){
			$minwin=$win if $self->{'windows'}->[$genome]->[$win] < $self->{'windows'}->[$genome]->[$minwin];
		}
		$self->{'win_thrs'}->[$genome] = $minwin/$self->{'bsrwins'};
		# Eval topology
		if($win_mean_log[$polygon[1]]<=$win_mean_log[$polygon[3]]
						&& $win_mean_log[$polygon[3]]<=$win_mean_log[$polygon[2]]
						&& $win_mean_log[$minwin]>=$win_mean_log[$polygon[3]]){
			$self->{'known_topology'}->[$genome] = 0;
			$self->{'unus'}->msg(3,"Unknown topology for genome $genome, setting threshold to the average");
		}else{
			$self->{'known_topology'}->[$genome] = 1;
			# Detect the lowest point within the window
			my $b = ($minwin/$self->{'bsrwins'}) * $self->{'bsrx'};
			my $a = $b - ($self->{'bsrx'}/$self->{'bsrwins'});
			my $minx = $a;
			for my $x ($a .. $b){ $minx = $x if $self->{'bsrhist'}->[$genome]->[$x] < $self->{'bsrhist'}->[$genome]->[$minx] }
			$self->{'thresholds'}->[$genome] = $minx/$self->{'bsrx'};
			$self->{'mean_thr'} += $self->{'thresholds'}->[$genome];
			$self->{'unus'}->msg(3,"Setting threshold for genome $genome to ".$self->{'thresholds'}->[$genome]);
			$self->{'unus'}->msg(4,"Internal range: ".$polygon[2]."-".$polygon[3]);
		}
		$self->{'polygons'}->[$genome] = [];
		push @{ $self->{'polygons'}->[$genome] }, $_/$self->{'bsrwins'} for @polygon;
	} # END GENOME
	my $N = 0;
	$N+=$_ for(@{ $self->{'known_topology'} });
	$self->{'mean_thr'} = $self->{'mean_thr'}/$N;
	for my $genome (0 .. $#{ $self->{'unus'}->{'genomes'} }){
		$self->{'thresholds'}->[$genome] = $self->{'mean_thr'} unless $self->{'known_topology'}->[$genome];
	}
	my $win_size = 1/$self->{'bsrwins'};
	my $dot_size = 1/$self->{'bsrx'};
	print BSRGRAPH "pdf(file='".$self->{'unus'}->{'basename'}.".bsrgraph.pdf');\n";
	for my $genome (0 .. $#{ $self->{'unus'}->{'genomes'} }){
		print BSRGRAPH 	"plot(idx,log(genome$genome),xlab='BitScore Ratio',ylab='".$self->{'genomes'}->[$genome].
						"',type='l',col='cornflowerblue');\n";
		print BSRGRAPH	"polygon(c(".($self->{'win_thrs'}->[$genome]-$win_size).",".$self->{'win_thrs'}->[$genome].",".
						$self->{'win_thrs'}->[$genome].",".($self->{'win_thrs'}->[$genome]-$win_size)."),".
						"c(-10,-10,max(log(genome$genome)*2),max(log(genome$genome)*2)),".
						",col='azure2',border=NA);\n" if $self->{'known_topology'}->[$genome];
		print BSRGRAPH	"lines(idx_win,log(win_genome$genome),col='darkorchid4');\n".
				"abline(v=".$self->{'thresholds'}->[$genome].",col='chartreuse3');\n".
				"abline(v=".$self->{'mean_thr'}.",col='darkorchid1',lty=2);\n".
				"abline(v=c(".join(",",@{ $self->{'polygons'}->[$genome] })."),col='grey');\n";
	}
	print BSRGRAPH "dev.off();\n";
	close BSRGRAPH;
	return;
}
1;
