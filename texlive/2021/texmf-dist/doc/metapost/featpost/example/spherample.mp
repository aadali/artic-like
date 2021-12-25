% spherample.mp
% L. Nobre G. 
% 2002

input featpost3Dplus2D;

SphericalDistortion := true;

beginfig(1);
    color orig, axxc, ayyc, azzc;
    color nalfa, gammacnt, vertex;
    numeric altura, grids, radius, xxsize, aux,
            yysize, zzsize, i, numc, slope;
    path sfear;
    
    f 			:= (5,4,2);
    Spread 		:= 120;
    ShiftV		:= (200,200);
    altura 		= 1;
    radius		= 0.3;
    grids		= 1;
    numc		= 9;
    slope		= 3;

    for i=1 upto numc:
        aux		:= (i-1-0.5*(numc-1))*grids;
        gammacnt 	:= (0,aux,0);
        vertex		:= (aux/slope,aux,altura);
	nalfa           := vertex - gammacnt;
	sfear		:= rigorousfearpath(gammacnt,radius);
        unfill sfear;
	draw sfear; 
	draw spatialhalfcircle(gammacnt,nalfa,radius,true);
    endfor;
    xxsize 		= radius;
    yysize 		= 0.5*(numc-1)*grids;
    zzsize 		= altura;
    orig 		= (0,0,0);
    axxc 		= (xxsize,0,0);
    ayyc 		= (0,yysize,0);
    azzc 		= (0,0,zzsize);
%    drawarrow rp(orig)..rp(axxc);
%    drawarrow rp(orig)..rp(ayyc);
%    drawarrow rp(orig)..rp(azzc);
endfig;

end;


