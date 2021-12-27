% rigorouscone.mp
% L. Nobre G. 
% 2003

input featpost3Dplus2D;

beginfig(1);
    color nalfa, gammacnt, vertex;
    path gamma, cone;
    numeric altura, grids, radius, aux, i, numc, slope;
    
    f 			:= (5,4,2);
    Spread 		:= 60;
    ShiftV		:= (200,200);
    altura 		= 1;
    radius		= 0.3;
    grids		= 1;
    numc		= 9;
    slope		= 3;
    nalfa 		= (0,0,1);

    for i=1 upto numc:
        aux		:= (i-0.5*(numc+1))*grids;
        gammacnt 	:= (0,aux,0);
        vertex		:= (aux/slope,aux,altura);
        gamma 		:= rigorouscircle( gammacnt, nalfa, radius );
        cone            :=   rigorouscone(true, gammacnt, nalfa, radius, vertex );
	unfill cone;
        draw gamma dashed evenly;
        draw cone;
    endfor;
endfig;

end;


