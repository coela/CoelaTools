package CoelaTools::CoOccurence;
use strict;
use PDL;
use Moose;

has 'colnames'         => ( isa => 'ArrayRef[Str]', is => 'ro' );
has 'table'            => ( isa => 'PDL',           is => 'ro',		required => 1 );
has 'sd_from_expected' => ( isa => 'PDL',           is => 'ro',		builder => '_build_sd_from_expected',	lazy_build => 1);
has 'pd_from_expected' => ( isa => 'PDL',           is => 'ro',		builder => '_build_pd_from_expected',	lazy_build => 1);
has 'oe'							 => ( isa => 'PDL',           is => 'ro',		builder => '_build_oe', lazy_build => 1 ); 
has '_index_position_by_colname' => ( isa => 'HashRef' , is => 'ro', builder => '_build_index_position_by_colname');

sub _build_oe {
	my $self = shift;
	warn "there is only one column in this piddle" unless ( scalar @{ $self->colnames } >= 2 );
	my @piddle;
	my $sum = sum $self->table;
	my $col = sumover( mv( $self->table, 1, 0 ) );
	my $row = transpose( sumover( $self->table ) );
	my $e = $sum * ($col / $sum) * $row / $sum ;
	my $oe = $self->table / ( $sum * $col / $sum * $row / $sum );
	return $oe;
}

sub _build_sd_from_expected {
	my $self = shift;
	warn "there is only one column in this piddle" unless ( scalar @{ $self->colnames } >= 2 );
	my @piddle;

	my $sum = sum $self->table;
	my $col = sumover( mv( $self->table, 1, 0 ) );
	my $row = transpose( sumover( $self->table ) );

	my $sd = ( $self->table - $sum * $col / $sum * $row / $sum ) / (
			sqrt(
				$sum * $col / $sum * $row / $sum * ( 1 - $col / $sum * $row / $sum )
				)
			);
	return $sd;
}

sub _build_pd_from_expected {
	my $self = shift;
	warn "there is only one column in this piddle" if ( scalar @{ $self->colnames } >= 2 );
	my @piddle;

	my $sum = sum $self->table;
	my $col = sumover( mv( $self->table, 1, 0 ) );
	my $row = transpose( sumover( $self->table ) );

	my $sd =
		100 *
		( $self->table - $sum * $col / $sum * $row / $sum ) /
		( $sum * $col / $sum * $row / $sum );
	return $sd;
}


sub _build_index_position_by_colname {
	my $self = shift;
	my $i = 0;
	my %position;

	for (@{$self->colnames}){
		$position{$_} = $i;
		$i++;
	}
	return \%position;
}

sub index_position_by_colname {
	my $self = shift;
	my $colname = shift;
	return ${$self->_index_position_by_colname}{$colname};
}

sub log {
	my $self = shift;
	say G::Seq::AminoAcid::one2three( $self->aminoacid ) . ":";
	say sum $self->table;
	say $self->table;
	if ( scalar @{ $self->colnames } >= 2 ) {
		say "SD from expected";
		say $self->sd_from_expected;
		say "PD from expected";
		say $self->pd_from_expected;
	}
}

__PACKAGE__->meta->make_immutable;
no Moose;

