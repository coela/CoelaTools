use 5.12.4;
use Fatal qw( open close );
use Data::Dumper;

package main;
use PDL;
use PDL::IO::Dumper;
use PDL::LiteF;
use PDL::NiceSlice;
use PDL::Stats::Basic;
use PDL::MatrixOps;
use Data::Dumper;
use G;

use File::Find::Rule;
use File::Copy;
use File::Basename;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

use lib './lib';
use tRNA;
use CoOccurence;
use TPI;
use FreeEnergy;

local $| = 1;

my $hogehoge;
my %oc;
my %t_table;

### create codon list ####
my @codons;
my @nc = ('A','T','G','C');
for my $a (@nc){
	for my $b(@nc){
		for my $c(@nc){
			push @codons, $a.$b.$c;
			$t_table{$a.$b.$c} = translate($a.$b.$c);
		}
	}
}
###########################


### read tRNA table ####
my %check;
open my $TRNA, $ARGV[0] or die;
while(<$TRNA>){
	chomp;
	my @line = split /\t/;
	push @{$check{$line[3]}}, $_; 
}
#########################

my @genome_dir = $ARGV[1] || ('/Users/coela/projects/kaiseki2_genomes');
my $rule =  File::Find::Rule->new;
$rule->file;
$rule->name( '*.gbk' );

my @files = $rule->in( @genome_dir );

my $heatmap;
for my $gf (@files){
###
	my $genome = basename $gf;
	my $gb;
	if (-f "/Users/coela/projects/codon_co_occurrence/store/$genome"){
		$gb = Storable::retrieve("/Users/coela/projects/codon_co_occurrence/store/$genome");
	}
	else{
		$gb = load($gf,"no msg");
		$gb = TPI::add_tpi($gb);
#	go_table() if $genome eq "NC_000913.gbk";
		cai($gb, -output=>"/dev/null");
		Storable::nstore $gb, "/Users/coela/projects/codon_co_occurrence/store/$genome";
	}
###

	my @tpi_array;
#	next unless defined $check{$gb->{FEATURE0}->{organism}};
	my $gcsi = gcsi($gb);
#say $gb->{FEATURE0}->{organism};
#say $gcsi;
	$gb->{FEATURE0}->{organism} = "Escherichia coli K12" if basename($gf) eq "NC_000913.gbk";

	next if $oc{$gb->{FEATURE0}->{organism}};
	$oc{$gb->{FEATURE0}->{organism}} = 1;
	next if $gb->{FEATURE0}->{chromosome};

	my $codon_co_occurrence;
	my $codon_tables = ();
	my $codon_check;

# count codon co occurences
	for my $cds ( $gb->feature('CDS') ) {
		next unless $gb->{$cds}->{transl_table} == 11;
		push @tpi_array, $gb->{$cds}->{tpi} if defined $gb->{$cds}->{tpi};

		next unless defined $gb->{$cds}->{translation};
		my $nuc_seq = $gb->get_geneseq($cds);
		my $translate = $gb->{$cds}->{translation};

		if ($shuffle) {
			my @words = $nuc_seq =~ m/.{1,3}/g;
			$nuc_seq = join "", shuffle @words;
			$translate = translate($nuc_seq);
		}

		my $i = 0;
		my $previous_codon_by_aa;

		for my $aa ( split //, $translate) {
			$aa = uc $aa;
		my $codon = uc substr $nuc_seq, $i * 3, 3;
		next if $codon =~ /[^ATGC]/;
		next unless $t_table{$codon} eq $aa;

		unless ($codon) {
			warn
				"error in codon position:\n >$cds\n $nuc_seq\n\n$gb->{$cds}->{translation}\n\n";
			next;
		}

		if ( defined $previous_codon_by_aa->{$aa} and $i != 0 ) {
			$codon_co_occurrence->{$aa}->{ $previous_codon_by_aa->{$aa} }->{$codon}++;
			$codon_check->{$aa}->{$codon} = 1;
		}
		$previous_codon_by_aa->{$aa} = $codon;
		$i++;
	}
	}

# codon count;
	for my $aminoacid ( sort keys %{$codon_co_occurrence} ) {
		my @co;
		my @codons = sort ( keys %{ $codon_check->{$aminoacid} } );
		for my $codon1 (@codons) {
			my @codon_line;
			for my $codon2 (@codons) {
				push @codon_line,
						 (
							$codon_co_occurrence->{$aminoacid}->{$codon1}->{$codon2}
							or 0
						 );
			}
			push @co, \@codon_line;
		}
		my $codon_co = CoOccurence->new(
				colnames  => \@codons,
				aminoacid => $aminoacid,
				table     => pdl @co,
				);

		push @$codon_tables, $codon_co;
# say $aminoacid;
#	say join " ", @{$codon_co->colnames};
#	say $codon_co->sd_from_expected;
#	say $codon_co->table;
#	say $codon_co->oe;
	}


# trna_count;
# init tRNA information from tRNACE-DB
	my @trnas;
	my $trnacount;
	for ( @{$check{$gb->{FEATURE0}->{organism}}} ) {
		chomp;
		my @line = split /\t/;
		my $trna = tRNA->new(
				aminoacid => G::Seq::AminoAcid::three2one( $line[6] ),
				anticodon => $line[7],
				);
		$trnacount->{$line[7]}++;
		push @trnas, $trna;
	}

	my %hash;
	my @aalist = ("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V");
#	my @aalist = ("S");

	my %skip;
	for my $aa (@aalist){
		my @trnas_match = grep {$_->aminoacid eq $aa} @trnas;
		my @codons_match = grep {$_->aminoacid eq $aa } @$codon_tables;
		my $codon_co = shift @codons_match;	
		my @trnas_match_nr;
		my %trna_hash;

		for my $trna (@trnas_match){
			push @trnas_match_nr, $trna unless $trna_hash{$trna->anticodon};
			$trna_hash{$trna->anticodon} = 1;
		}

		my %match;
		for my $codon (@{$codon_co->colnames}) {
			my @direct_trna_partners = grep {$codon eq $_->translating_codon } @trnas_match_nr;
			my @w_trnas = grep { grep{$_ eq $codon} @{$_->wobble} }@trnas_match_nr;
#			say "$aa $codon" if scalar @w_trnas >= 2;
			if (@direct_trna_partners) {
				my $direct_trna_partner = shift @direct_trna_partners;
				push @{$match{$direct_trna_partner->anticodon}}, $codon;
			}
			else{
				die if (scalar @w_trnas >= 2);
				my $wobbling_trna_partner = shift @w_trnas;
				unless ($wobbling_trna_partner){
					$skip{$aa} = 1;
					next;
				}
				push @{$match{$wobbling_trna_partner->anticodon}}, $codon;
			}
		}
		if ($skip{$aa}){
#say "$aa skipped";
#next;
		}

		my $i = -1;
		my %position;
		for my $codon (@{$codon_co->colnames}) {
			$i++;
			$position{$codon} = $i;
		}
		my @trna_table;
		for my $trna (keys %match) {
			my $hoge;
			for (@{$match{$trna}}){
				$hoge += $codon_co->table->slice(":,($position{$_})");
			}
			push @trna_table, $hoge;
		}

		my $trna_pdl = transpose pdl @trna_table;
		my @trna_table2;
		for my $trna (keys %match) {
			my $hoge;
			for (@{$match{$trna}}){
				$hoge += $trna_pdl->slice(":,($position{$_})");
			}
			push @trna_table2, $hoge;
		}

		next if (scalar keys %match) == 1;
		my $trna_co = CoOccurence->new(
				colnames  => [keys %match],
				aminoacid => $aa,
				table     => (transpose pdl @trna_table2),
				);

#		say '-----------------';
		my @anticodon_to_codon;
		for my $anti (@{$trna_co->colnames}){
			my $old_anti = $anti;
			$anti =~ tr/ATGC/TACG/;
			$anti = reverse $anti;
			push @anticodon_to_codon, $anti;
#			say $anti, $trnacount->{$old_anti};
		}
		for (@{$codon_co->colnames}){
			my $v = $codon_co->sd_from_expected->diagonal(0,1)->index("$position{$_}")->sum;
			die join "\t" ,($gb->{FEATURE0}->{organism},$codon_co->aminoacid,$_,$heatmap->{$gb->{FEATURE0}->{organism}}->{$_},$v) if defined $heatmap->{$gb->{FEATURE0}->{organism}}->{$_};
			$heatmap->{$gb->{FEATURE0}->{organism}}->{$_} = $v;
		}
		say $codon_co->aminoacid;
		say Dumper $codon_co->colnames;
#	
##		say join "\t", ($gb->{FEATURE0}->{organism},$gcsi,$trna_co->sd_from_expected->diagonal(0,1)->average,$trna_co->sd_from_expected->diagonal(0,1)->stdv);
#		say $trna_co->table;
#	say $trna_co->sd_from_expected;
#		say "e";
		say $codon_co->aminoacid;
		say Dumper $codon_co->colnames;
#		say "OE";
#		say $codon_co->oe;
#		say "table";
#		say $trna_co->table;
#		say "e";
#		say $trna_co->e;
#		say "OE";
#		say $trna_co->oe;
	}
=c
		my $result = "/Users/coela/projects/codon_co_occurrence/tmp/$genome.tpi";
	open my $TPILIST, '>' , $result;
	say $TPILIST $_ for @tpi_array;
	close $TPILIST;

	my $script = "/tmp/rscript$$.R";
	open(SCRIPT, ">$script");

	print SCRIPT <<EOF;
	tpi = read.table("$result")
		png("/Users/coela/projects/codon_co_occurrence/fig/tpi/$genome.png")
		hist(tpi[,1],main="TPI $genome")
		EOF

		close(SCRIPT);
	system("R --vanilla --slave < $script");
	my $Rstat;
	($Rstat = $?/256) && die "Aborted in Darwin with status $Rstat.?n";
	unlink $script;
	=cut
		$hogehoge++;
	say STDERR $hogehoge;
}
print "\t";
say join "\t",@codons;
for my $s (keys %{$heatmap}){
	my @vals = map {$heatmap->{$s}->{$_}||'NA'} @codons;
	print "$s\t";
	say join "\t", @vals;
}
