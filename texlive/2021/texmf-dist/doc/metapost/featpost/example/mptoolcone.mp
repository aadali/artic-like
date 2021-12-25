% mptoolcone.mp
% L. Nobre G. 
% Y2K

input featpost3Dplus2D;

%input mp-tool;

	defaultscale := 1.1;
%	defaultfont := "putr8r";
%	defaultfont := "bchr8r";

beginfig(0);
    color orig, axxc, ayyc, azzc;
    color nalfa, gammacnt;
    path gamma, cone, mybox;
    numeric theta, radius, xxsize, yysize, zzsize;
    picture mypict;

    f 			:= (3,4,2);
    Spread 		:= 100;
    xxsize 		= 0.8;
    yysize 		= 0.8;
    zzsize 		= 1.0;
    theta 		= 25;

    orig 		= (0,0,0);
    mypict = thelabel("Texto com acentuação.",rp(orig));
    mybox = bbox mypict;
    draw mypict;
    draw mybox;
    axxc 		= (xxsize,0,0);
    ayyc 		= (0,yysize,0);
    azzc 		= (0,0,zzsize);
    nalfa 		= (0,sind(theta),cosd(theta));
    gammacnt 	 	= cosd(theta)*cosd(theta)*nalfa;
    radius 		:= cosd(theta)*sind(theta);
    gamma 		= rigorouscircle( gammacnt, nalfa, radius );
    cone      = rigorouscone(true, gammacnt, nalfa, radius, (xxsize,yysize,0) );
    draw gamma dashed evenly;
    draw cone;
    drawarrow rp(orig)..rp(axxc);
    drawarrow rp(orig)..rp(ayyc);
    drawarrow rp(orig)..rp(azzc);
endfig;

end;


