% nemuniax.mp
% L. Nobre G.
% 2003

input featpost3Dplus2D;

beginfig(1);
	f := 2*(5,3,4);
	Spread := 70;

	numeric hstep, hmax, hsdev, basestep, basemax, basesdev;
	numeric i, j, k, angsdev, cylength, cyradius, basen, hn;
	numeric vcx, vcy, vcz, hcurr, xcurr, ycurr, hbase, aone;
	numeric atwo, zcurr, counter;
	color lenvec, currpos;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Nu
	cylength   := 0.6;
	cyradius   := 0.1;
	basen	   := 12;
	hn	   := 5;
	basestep   := cyradius*3.1;
	hstep	   := cylength*1.5;
	basesdev   := cyradius*0.1;
	hsdev	   := hstep*0.2;
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
		    aone    :=               uniformdeviate( 360 );
		    atwo    :=               angsdev*normaldeviate;
		    vcz     := cosd( atwo );
		    vcy     := sind( atwo )*sind( aone );
		    vcx     := sind( atwo )*cosd( aone );
		    currpos := ( xcurr, ycurr, zcurr );
		    lenvec  := cylength*( vcx, vcy, vcz ); 
%		    rigorousdisc( 0, true, currpos-0.5*lenvec, cyradius, lenvec);
		    counter := incr( counter );
		    generatedirline( counter, aone, 90-atwo, cylength, currpos );
%		    generateonebiax( counter, aone, 90-atwo, cylength,
%                                              -60, 0.25*cylength, currpos );
		endfor;
	    endfor;
	endfor;
	NL := counter;
	director_invisible( true, 5, false );
%	director_invisible( false, 5, true );
endfig;

end;
