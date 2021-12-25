% nembends.mp
% L. Nobre G.
% 2002

input featpost3Dplus2D;

%input mp-tool;

beginfig(1);
	f := 2*(5,3,4);
	Spread := 70;

	numeric hsdev, basestep, basxstep, basesdev;
	numeric i, j, k, coesdev, cylength, cyradius, basen;
	numeric vcx, vcy, vcz, hcurr, xcurr, ycurr, aone;
	numeric atwo, zcurr, coestep, angsdev, basxsdev;
	color lenvec, currpos;


	coestep    := -0.15; 
	cylength   := 0.6;
	cyradius   := 0.1;
	basen	   := 3;
	basestep   := cylength*1.6;
	basxstep   := cyradius*4;
	basesdev   := basestep*0.001;
	basxsdev   := basxstep*0.001;
	hsdev	   := 0.005;
	coesdev	   := 0.001;
	angsdev    := 0.02;
	for i=-basen upto basen:
	    for j=-basen upto basen:
		    xcurr   :=        i*basxstep + basxsdev*normaldeviate;
		    ycurr   :=        j*basestep + basesdev*normaldeviate;
		    zcurr   :=   ycurr**2*coestep + coesdev*normaldeviate;
		    aone    :=                      angsdev*normaldeviate;
		    atwo :=angle(1,ycurr*2*coestep)+angsdev*normaldeviate;
		    vcz     := cosd( aone )*sind( atwo );
		    vcy     := cosd( aone )*cosd( atwo ); 
		    vcx     := sind( aone ); 
		    currpos := ( xcurr, ycurr, zcurr );
		    lenvec  := cylength*( vcx, vcy, vcz ); 
		    rigorousdisc(0,true,currpos-0.5*lenvec,cyradius,lenvec);
	    endfor;
	endfor;
endfig;

end;
