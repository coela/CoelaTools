package TPI;
use strict;
use Data::Dumper;
use G;
use tRNA;

sub tpi {
	my $self = shift;
	my %opts = ( 
			debug => 0,
			darwin_exec => "~/projects/darwinmac64/darwin",
			darwin_lib => "~/projects/darwinmac64/lib",
			tmp_sgml => "/tmp/genomecds$$.sgml",
			tmp_darwin => "/tmp/darwinrap$$.R",
			tmp_result => "/tmp/darwintpi$$.tab",
			darwin_dir => "/tmp/genomecds$$.sgml",
			@_,
			);
	
	open my $SGML, ">$opts{tmp_sgml}" or die;
	system('unlimit stacksize');

	my $i;
	if ($sequence){
		print $SGML "<E>";
		print $SGML "<ID>Hoge</ID>";
		print $SGML "<DNA>". uc $sequence."</DNA>";
		print $SGML "<SEQ>". uc $sequence."</SEQ>";
		print $SGML "</E>\n";
	}
	else{
		for my $cds ($self->cds()){
			next unless $self->{$cds}->{transl_table} == 11;
			if (substr ($self->{$cds}->{translation},-1,1) eq "\/"){
				chop $self->{$cds}->{translation};
			}

#omit cds with ! seqs <-> codon
			next if $self->{$cds}->{translation} =~ /[^ARNDCEQGHILKMFPSTWYV]/;

			my $nuc_seq = $self->get_geneseq($cds);
			unless ($nuc_seq){
				warn "$cds has no geneseq";
				next;
			}

			my $gene_length = length $nuc_seq;
			next if $gene_length % 3;
			next if $nuc_seq =~ /[^atgc]/;
			$i++;
			print $SGML "<E>";
			print $SGML "<ID>$cds</ID>";
			print $SGML "<DNA>". uc $nuc_seq."</DNA>";
			print $SGML "<SEQ>". uc $nuc_seq."</SEQ>";
			print $SGML "</E>\n";
		}
	}
	close($SGML);
	my $opt{tmp_darwin} = "/tmp/darwinrap$$.R";
	my $opt{tmp_result} = "/tmp/darwintpi$$.tab";
	open(my $SCRIPT, ">$opt{tmp_darwin}");

	print $SCRIPT <<EOF;
		DB = ReadDb('/tmp/genomecds$$.sgml');
		SetuptRNA(Bacteria);
		output := CreateArray(1..$i,[]);
		for i from 1 by 1 to $i do 
			tpi := ComputeTPI(Entry(i));
			output[i] := [ID(Entry(i)),tpi[1],tpi[2]];
		od;
		WriteData(output,'$opt{tmp_result}','\t');
		quit();
EOF

	close($SCRIPT);
	if ($opts{debug} == 1){
		print 'script start\n';
		system ("$opts{darwin_exec} -l $opts{darwin_lib} -i $opt{tmp_darwin}");
	}
	else{
		system ("$opts{darwin_exec} -l $opts{darwin_lib} -i $opt{tmp_darwin} -o /dev/null");
	}
	my $Dstat;
	($Dstat = $?/256) && die "Aborted in Darwin with status $Dstat.?n";
	unlink "$opts{tmp_sgml}.map";
	unlink "$opts{tmp_sgml}.tree";
	unlink $opt{tmp_darwin};
	unlink ($opts{tmp_sgml});

	open my $TPI, $opt{tmp_result} or die;
	while(<$TPI>){
		chomp;
		my ($cds,$tpi,$tpi2) = split /\t/;
		$cds =~ /ID\((.+)\)/;
		$self->{$1}->{'tpi'} = $tpi;
		$self->{$1}->{'tpi_l'} = $tpi2;
	}
	unlink $opt{tmp_result};
	return $self;
}

1;
