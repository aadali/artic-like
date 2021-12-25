% nemquira.mp
% L. Nobre G.
% 2002

input featpost3Dplus2D;

 

beginfig(1);
	f := 7*(1,2,1);
	Spread := 70;

	numeric hstep, hmax, hsdev, basestep, basemax, basesdev;
	numeric i, j, k, angsdev, cylength, cyradius, basen, hn;
	numeric vcx, vcy, vcz, hcurr, xcurr, ycurr, hbase, aone;
	numeric atwo, zcurr, angstep, counter;
	color lenvec, currpos;

	angstep    := 20; 
	cylength   := 0.45;
	cyradius   := 0.1;
	basen	   := 4.5;
	hn	   := 1;
	basestep   := cylength*1.37;
	hstep	   := basestep;
	basesdev   := basestep*0.05;
	hsdev	   := basesdev;
	hbase	   := -0.8;
	angsdev	   := 4;

	basemax	   := basen*basestep;
	hmax	   := hn*hstep;
	counter    := 0;
	for i=-basen upto basen:
	    for j=-basen upto basen:
		hcurr := hbase;
		forever:
		    exitunless hcurr < hmax+hbase;
		    hcurr   := hcurr + hstep;
		    zcurr   :=                hcurr + hsdev*normaldeviate;
		    xcurr   :=        i*basestep + basesdev*normaldeviate;
		    ycurr   :=        j*basestep + basesdev*normaldeviate;
		    aone    :=                      angsdev*normaldeviate;
		    atwo    :=        j*angstep  +  angsdev*normaldeviate;
		    vcz     := sind( atwo+90 );
		    vcy     := cosd( atwo+90 )*sind( aone );  
		    vcx     := cosd( atwo+90 )*cosd( aone );   
		    currpos := ( xcurr, ycurr, zcurr );
		    lenvec  := cylength*( vcx, vcy, vcz ); 
%		    rigorousdisc( 0, true, currpos-0.5*lenvec, cyradius, lenvec);
		    counter := incr(counter);
		    generatedirline( counter, aone, 90+atwo, cylength, currpos );
		endfor;
	    endfor;
	endfor;
	NL := counter;
	director_invisible( true, 5, false );
endfig;

end;
