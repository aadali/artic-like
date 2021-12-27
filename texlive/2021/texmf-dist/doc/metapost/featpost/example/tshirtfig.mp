% tshirtfig.mp
% L. Nobre G.
% 2004

input featpost3Dplus2D;

beginfig(1);
	f := (5,4,3);
	Spread := 70;
	pickup pencircle scaled 2.4pt;
	numeric gridstep, cylength, cyradius, numcylis;
	numeric i, j, curcoord, curangle, curnorma, tmp;
	numeric othcoord, othangle, zcoord, coef, othnorma;
	color base, lenvec;
	cylength := 0.2;
	cyradius := 0.3;
	gridstep := cyradius*2.2;
	numcylis := 7;
	coef := 0.22;
	for i=1 upto numcylis:			% saddle
	    curnorma := i-0.5*(numcylis+1);
	    curcoord := curnorma*gridstep;
	    curangle := -curcoord*2*coef;
	    for j=1 upto numcylis:
		othnorma := j-0.5*(numcylis+1);
		othcoord := othnorma*gridstep;
		othangle := othcoord*2*coef;
		zcoord:=coef*(othcoord**2-curcoord**2);
		base := (curcoord,othcoord,zcoord);
		lenvec := cylength*N((-curangle,-othangle,1));
		tmp := zcoord/2.4+0.8;
		if tmp > 1: show tmp; fi;
		background := tmp*TableC2;
		rigorousdisc( 0, true,base,cyradius,lenvec);
	    endfor;
	endfor;
endfig;

end.