% nsmetica.mp
% L. Nobre G.
% 2002

input featpost3Dplus2D;

 

beginfig(1);
	f := 2*(5,3,2);
	Spread := 70;
	NF := 0;

	numeric hstep, hmax, hsdev, basestep, basemax, basesdev;
	numeric i, j, k, angsdev, cylength, cyradius, basen, hn;
	numeric vcx, vcy, vcz, hcurr, xcurr, ycurr, hbase, aone;
	numeric atwo, zcurr, counter;
	color lenvec, currpos;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sma
	cylength   := 0.45;
	cyradius   := 0.1;
	basen	   := 11;
	hn	   := 3;
	basestep   := cyradius*2.4;
	hstep	   := cylength*2.1;
	basesdev   := cyradius*0.3;
	hsdev	   := hstep*0.04;
	hbase	   := -0.8;
	angsdev	   := 7;

	basemax	   := basen*basestep;
	hmax	   := hn*hstep;
	hcurr      := hbase;
	counter    := 0;
	for k=1 upto hn:
	  hcurr    := hcurr + hstep;
          for i=1 upto basen:
	    for j=1 upto basen:
	      zcurr   :=         hcurr + hsdev*normaldeviate;
	      xcurr   := (i-1)*basestep + uniformdeviate( basestep );
	      ycurr   := (j-1)*basestep + uniformdeviate( basestep );
	      aone    :=               uniformdeviate( 360 );
	      atwo    :=               angsdev*normaldeviate;
	      vcz     := cosd( atwo );
	      vcy     := sind( atwo )*sind( aone );
	      vcx     := sind( atwo )*cosd( aone );
	      currpos := ( xcurr, ycurr, zcurr );
	      lenvec  := cylength*( vcx, vcy, vcz ); 
	      counter := incr( counter );
	      generatedirline( counter, aone, 90-atwo, cylength, currpos );
%	      rigorousdisc( 0, true, currpos-0.5*lenvec, cyradius, lenvec);
	    endfor;
	  endfor;
	endfor;
	NL := counter;
	director_invisible( true, 5, false );
endfig;

end;
