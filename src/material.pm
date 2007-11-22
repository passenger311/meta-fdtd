
package material;
require Exporter;

our @ISA       = qw(Exporter);
our @EXPORT    = qw(epsilon,background,sphere,cylinder,background,layer,write_epsilon);
our $VERSION   = 1.00;


BEGIN {
    $background = 0;
    $num_layer = 0;
    $num_sphere = 0;
    $num_cylinder = 0; 
    $property{epsilon} = 1.0;
    $layer[0] = 0;
};

sub epsilon {
    $property{epsilon} = shift @_;
}

sub write_epsilon_1 {

    my $out = shift;

    $eps = $background->{prop}{epsilon};
    print $out "#INF\n";                    
    print $out "$eps\t\t! background epsilon\n";
    print $out "$num_layer\t\t! number of layers\n";
    print $out "$num_cylinder\t\t! number of cylinders\n";
    print $out "$num_sphere\t\t! number of spheres\n";
    print $out "\n";

    print $out "#LAY\n";                    
    for ($i=1; $i<=$num_layer; $i++){    ! print layers
	$lay = $layer[$i];
	$i0 = $lay->{i0};
	$i1 = $lay->{i1};
	$j0 = $lay->{j0};
	$j1 = $lay->{j1};
	$k0 = $lay->{k0};
	$k1 = $lay->{k1};
	$eps = $lay->{prop}{epsilon};
	print $out "$i0 $i1 $j0 $j1 $k0 $k1 $eps\n";
    }
};

sub background {
    my( @p ) = @_;
    
    $background = { prop => { %property } };
};


sub layer {

    my ($i0,$i1,$j0,$j1,$k0,$k1) = @_;

    $num_layer++;
    %box= ( i0 => $i0, i1 => $i1,
	    j0 => $j0, j1 => $j1,
	    k0 => $k0, k1 => $k1);
    $lay = {
	prop => { %property },
	box  => { %box },
	i0 => $i0, i1 => $i1,
	j0 => $j0, j1 => $j1,
	k0 => $k0, k1 => $k1,
    };
    $layer[$num_layer] = $lay;

};


sub sphere {



};


sub cylinder {



};




1;
