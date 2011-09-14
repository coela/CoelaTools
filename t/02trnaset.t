#t/01trna.t
use Data::Dumper;
use strict;
use Test::More;
plan (tests => 5);

use_ok("CoelaTools::tRNASet");
my $hoge;
while(<DATA>){
	$hoge .= $_;
}

my $trnaset = CoelaTools::tRNASet->new( source => $hoge, string_data => 1);
my @methods = (
		'copynumbers_by_aa',
		'translatable_tRNAs_by_codon',
		'translatable_tRNAs_by_codon_nr',
		'tRNAs_by_aa_nr',
		'tRNAs_by_aa');
isa_ok($trnaset,'CoelaTools::tRNASet');
can_ok('CoelaTools::tRNASet', @methods);
is($trnaset->copynumbers_by_aa('S'), 5, 'function copynumbers_by_aa');
is(scalar @{$trnaset->tRNAs}, 86, 'number of tRNAs in $obj->tRNAs');

__DATA__
>C008587	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	225381	225457	Ile	GAT
>C008588	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	225500	225575	Ala	TGC
>C008589	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	228928	229004	Asp	GTC
>C008590	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	236931	237007	Asp	GTC
>C008591	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	262095	262170	Thr	CGT
>C008592	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	563946	564022	Arg	TCT
>C008593	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	779777	779852	Lys	TTT
>C008594	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	779988	780063	Val	TAC
>C008595	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	780066	780141	Lys	TTT
>C008596	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	780291	780366	Val	TAC
>C008597	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	780370	780445	Lys	TTT
>C008598	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	780592	780667	Lys	TTT
>C008599	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	780800	780875	Lys	TTT
>C008600	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	1744459	1744535	Val	GAC
>C008601	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	1744540	1744616	Val	GAC
>C008602	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2042573	2042648	Asn	GTT
>C008603	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2057875	2057950	Asn	GTT
>C008604	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2060284	2060359	Asn	GTT
>C008605	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2284233	2284309	Pro	GGG
>C008606	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2464331	2464405	Arg	CCT
>C008607	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2518953	2519028	Val	TAC
>C008608	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2519073	2519148	Val	TAC
>C008609	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2519195	2519270	Val	TAC
>C008610	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2519275	2519350	Lys	TTT
>C008611	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2945409	2945485	Met	CAT
>C008612	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2945519	2945595	Met	CAT
>C008613	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2945629	2945705	Met	CAT
>C008614	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	3108388	3108463	Phe	GAA
>C008615	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	3213620	3213695	Met	CAT
>C008616	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	3834245	3834335	SeC	TCA
>C008617	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	3941458	3941533	Glu	TTC
>C008618	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	3944895	3944971	Asp	GTC
>C008619	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	3944980	3945055	Trp	CCA
>C008620	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	3980398	3980474	Arg	CCG
>C008621	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	3980533	3980608	His	GTG
>C008622	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	3980629	3980715	Leu	CAG
>C008623	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	3980758	3980834	Pro	TGG
>C008624	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	4035164	4035240	Ile	GAT
>C008625	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	4035283	4035358	Ala	TGC
>C008626	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	4166395	4166470	Glu	TTC
>C008627	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	4173411	4173486	Thr	TGT
>C008628	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	4173495	4173579	Tyr	GTA
>C008629	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	4173696	4173770	Gly	TCC
>C008630	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	4173777	4173852	Thr	GGT
>C008631	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	4207797	4207872	Glu	TTC
>C008632	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	4390383	4390458	Gly	GCC
>C008633	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	4390495	4390570	Gly	GCC
>C008634	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	4390606	4390681	Gly	GCC
>C008635	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	4494428	4494512	Leu	CAA
>C008636	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	4604424	4604338	Leu	CAG
>C008637	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	4604309	4604223	Leu	CAG
>C008638	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	4604188	4604102	Leu	CAG
>C008639	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	4360649	4360574	Phe	GAA
>C008640	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	3706715	3706639	Pro	CGG
>C008641	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	3425174	3425098	Ile	GAT
>C008642	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	3425055	3424980	Ala	TGC
>C008643	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	3421677	3421602	Thr	GGT
>C008644	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	3320180	3320094	Leu	GAG
>C008645	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	3316311	3316235	Met	CAT
>C008646	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2997079	2997006	Gly	CCC
>C008647	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2816667	2816575	Ser	GCT
>C008648	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2816571	2816495	Arg	ACG
>C008649	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2816296	2816220	Arg	ACG
>C008650	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2816157	2816081	Arg	ACG
>C008651	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2815882	2815806	Arg	ACG
>C008652	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2783859	2783784	Met	CAT
>C008653	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2727466	2727391	Glu	TTC
>C008654	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2516253	2516178	Ala	GGC
>C008655	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2516138	2516063	Ala	GGC
>C008656	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2056126	2056051	Asn	GTT
>C008657	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	2041581	2041492	Ser	CGA
>C008658	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	1990141	1990066	Gly	GCC
>C008659	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	1990011	1989938	Cys	GCA
>C008660	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	1989925	1989839	Leu	TAA
>C008661	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	1286845	1286761	Tyr	GTA
>C008662	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	1286551	1286467	Tyr	GTA
>C008663	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	1096875	1096788	Ser	GGA
>C008664	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	1030935	1030848	Ser	TGA
>C008665	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	925194	925107	Ser	GGA
>C008666	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	696356	696280	Met	CAT
>C008667	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	696270	696186	Leu	TAG
>C008668	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	696162	696088	Gln	TTG
>C008669	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	696053	695979	Gln	TTG
>C008670	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	695963	695887	Met	CAT
>C008671	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	695839	695765	Gln	CTG
>C008672	Ecol_K12_MG1655	Gammaproteobacteria	Escherichia coli K12	695727	695653	Gln	CTG
