% nembiaxi.mp
% L. Nobre G.
% 2002

input featpost3Dplus2D;

 

beginfig(1);
	f := 1.3*(4,5,3);
	Spread := 70;

	numeric hstep, hmax, hsdev, basestep, basemax, basesdev;
	numeric i, j, k, angsdev, cylength, cyradius, basen, hn;
	numeric vcx, vcy, vcz, hcurr, xcurr, ycurr, hbase, aone;
	numeric atwo, zcurr, counter;
	color lenvec, currpos, xnoise, tmpcol;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Nb
	cylength   := 0.45;
	cyradius   := 0.1;
	basen	   := 5;
	hn	   := 3;
	basestep   := cylength*1.65;
	hstep	   := cylength;
	basesdev   := basestep*0.1;
	hsdev	   := hstep*0.1;
	hbase	   := -0.8;
	angsdev	   := 12;

	basemax	   := basen*basestep;
	hmax	   := hn*hstep;
	counter    := 0;
	for i=1 upto basen:
	    for j=1 upto basen:
		hcurr := hbase;
		forever:
		    exitunless hcurr < hmax+hbase;
		    hcurr   := hcurr + hstep;
		    zcurr   :=         hcurr + hsdev*normaldeviate;
		    xcurr   := i*basestep + basesdev*normaldeviate;
		    ycurr   := j*basestep + basesdev*normaldeviate;
		    aone    :=              uniformdeviate( 10 )-5;
		    atwo    :=               angsdev*normaldeviate;
		    vcz     := cosd( atwo );
		    vcy     := sind( atwo )*sind( aone );
		    vcx     := sind( atwo )*cosd( aone );
		    currpos := ( xcurr, ycurr, zcurr );
		    xnoise  := ( sind( angsdev*normaldeviate ), 0, 0 );
		    xnoise  := green + xnoise;
		    tmpcol  := ncrossprod( ( vcx, vcy, vcz ), xnoise );
		    lenvec  := cylength*tmpcol;
%		    rigorousdisc( 0, true, currpos-0.5*lenvec, cyradius, lenvec);
		    counter := incr( counter );
		    generateonebiax( counter, aone, atwo, cylength,
                                              90, 0.35*cylength, currpos );
		endfor;
	    endfor;
	endfor;
	NL := counter;
	director_invisible( false, 7, true );
endfig;

end;
