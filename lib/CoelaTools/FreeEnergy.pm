package FreeEnergy;
use strict;
use G;

sub add_free_energy{
	my $self = shift;

	my %opts = (
			exec_command => 'hybrid-ss-min',
			before => 4,
			after => 37,
			@_,
			);

	for my $cds ($self->feature('CDS')){
		my $before = $self->before_startcodon($cds,$opts{before});
		my $after = substr ($self->get_geneseq($cds), 0 , $opts{after});
		my $sequence = $before.$after;
		my $freeenergy = `echo "$sequence" | $opts{exec_command} --stream`;
		chomp $freeenergy;
		$self->{$cds}->{'free_energy'} = $freeenergy;
	}
	return $self;
}

1;

