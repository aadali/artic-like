% rigorousdiscSD.mp
% L. Nobre G. 
% 2003

input featpost3Dplus2D;

f := 0.8*(6,3,4);

SphericalDistortion := true;
%%%%%Spread	     := 6;      

%SphericalDistortion := false;
%%%%%%Spread	     := 20; 

beginfig(1);
    color nalfa, gammacnt, vertex;
    numeric altura, grids, radius, aux, auy, i, j, numc;
    
    altura 		= 1;
    radius		= 0.5;
    grids		= 1.4;
    numc		= 7;

    for i=1 upto numc:
        aux := (i-0.5*(numc+1))*grids;
	for j=1 upto numc:
	    auy	:= (j-0.5*(numc+1))*grids;
	    if normaldeviate > 0:
	        gammacnt := (aux,auy,0);
		vertex	 := (aux,auy,altura);
	    else:
		if normaldeviate > 0:
	            gammacnt := (aux-altura/2,auy,radius);
		    vertex   := (aux+altura/2,auy,radius);
		else:
	            gammacnt := (aux,auy-altura/2,radius);
		    vertex   := (aux,auy+altura/2,radius);
		fi;
	    fi;
	    nalfa := vertex - gammacnt;
            rigorousdisc( 0, true, gammacnt, radius, nalfa );
	endfor;
    endfor;
    produce_auto_scale;
endfig;

end;


