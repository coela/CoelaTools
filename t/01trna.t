#t/01trna.t
use strict;
use Test::More;
plan (tests => 5);
use_ok("CoelaTools::tRNA");

my $trna = CoelaTools::tRNA->new( anticodon=>'TCA', aminoacid =>'S');
isa_ok($trna,'CoelaTools::tRNA');
is ($trna->complementary_codon, 'TGA',"complementary_codon for TGA OK");
my @methods = (
		'translation_possible'
		);
can_ok('CoelaTools::tRNA', @methods);
is_deeply ($trna->translatable_codons, ['TGA','TGG'], "wobbling check OK");

TODO: {
#tRNAとかの全部チェック
			}

