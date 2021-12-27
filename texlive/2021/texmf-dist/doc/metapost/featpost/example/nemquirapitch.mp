% nemquirapitch.mp
% L. Nobre G.
% 2002

input featpost3Dplus2D;

verbatimtex 
\documentclass{article} 
\usepackage{beton} 
\usepackage{concmath} 
\usepackage{ccfonts} 
\begin{document} 
etex

beginfig(1);
	f := 7*(1,2,1);
	Spread := 70;
	ahlength := 3mm;
	ahangle := 48;

	numeric hsdev, basestep, basemax, basesdev, xpos;
	numeric i, j, angsdev, cylength, cyradius, basen;
	numeric vcx, vcy, vcz, hcurr, xcurr, ycurr, aone;
	numeric atwo, zcurr, angstep, counter, linelen;
	color lenvec, currpos;

	angstep    := 20; 
	cylength   := 0.45;
	cyradius   := 0.1;
	basen	   := 4.5;
	basestep   := cylength*1.37;
	basesdev   := basestep*0.05;
	hsdev	   := basesdev;
	angsdev	   := 4;
	linelen    := cylength;

	basemax	   := basen*basestep;
	xpos       := -(basen+1)*basestep;
	draw rp((-xpos-0.5*cylength,-basemax,0))--rp((-xpos+0.5*cylength,-basemax,0));
	draw rp((-xpos-0.5*cylength,basemax,0))--rp((-xpos+0.5*cylength,basemax,0));
	drawdblarrow rp((-xpos,-basemax,0))--rp((-xpos,basemax,0));
	label.urt( btex $\frac{p}{2}$ etex scaled 2.5, rp((-xpos,0,0)) );
	counter    := 0;
	for i=-basen upto basen:
	    for j=-basen upto basen:
		    zcurr   :=                        hsdev*normaldeviate;
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
	NL := counter;
	director_invisible( true, 5, false );
endfig;

verbatimtex 
\end{document} 
etex

end;
