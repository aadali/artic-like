% nemdirector.mp
% L. Nobre G.
% 2013

input featpost3Dplus2D;

verbatimtex 
\documentclass{article} 
\usepackage{lucidabr} 
\begin{document} 
etex

beginfig(1);
	f := 2*(5,3,4);
	Spread := 70;

	numeric cylength, cyradius, zoneradi, zoneleng;
	numeric i, angsdev, halfnum, directlen, arrthick;
	numeric vcx, vcy, vcz, xcurr, ycurr, thickfact;
	numeric atwo, zcurr, aone, apos, anul;
	color lenvec, currpos, tippos;
	path arrpath;
	picture tmppict;
        anul       := 0;
	angsdev	   := 7;
	cylength   := 0.45;
	cyradius   := 0.1;
	zoneradi   := 2;
	zoneleng   := 4;
	directlen  := 5.6;
	halfnum    := 130;
	thickfact  := 5;
	arrthick   := Spread*ps( (0,0,0), thickfact );
	tippos     := (0,0,0.5*directlen);
	arrpath    := rp(-tippos)--rp(tippos);
	ahlength   := 3*arrthick;
	ahangle    := 40;
	for i=1 upto halfnum:
	  zcurr   := -0.5*zoneleng+uniformdeviate( zoneleng );
	  apos    := uniformdeviate( 180 );
	  xcurr   := -zoneradi*cosd( apos )*uniformdeviate( 1 );
	  ycurr   := zoneradi*sind( apos )*uniformdeviate( 1 );
	  aone    := uniformdeviate( 360 );
	  atwo    := anul+angsdev*normaldeviate;
	  vcz     := cosd( atwo );
	  vcy     := sind( atwo )*sind( aone );
	  vcx     := sind( atwo )*cosd( aone );
	  currpos := ( xcurr, ycurr, zcurr );
	  lenvec  := cylength*( vcx, vcy, vcz ); 
	  generatedirline( i, aone, 90-atwo, cylength, currpos );
	endfor;
	NL:= halfnum;
	director_invisible( true, thickfact, false );
	drawarrow arrpath withpen pencircle scaled arrthick;
	for i=1 upto halfnum:
	  zcurr   := -0.5*zoneleng+uniformdeviate( zoneleng );
	  apos    := uniformdeviate( 180 );
	  xcurr   := zoneradi*cosd( apos )*uniformdeviate( 1 );
	  ycurr   := zoneradi*sind( apos )*uniformdeviate( 1 );
	  aone    := uniformdeviate( 360 );
	  atwo    := anul+angsdev*normaldeviate;
	  vcz     := cosd( atwo );
	  vcy     := sind( atwo )*sind( aone );
	  vcx     := sind( atwo )*cosd( aone );
	  currpos := ( xcurr, ycurr, zcurr );
	  lenvec  := cylength*( vcx, vcy, vcz ); 
	  generatedirline( i, aone, 90-atwo, cylength, currpos );
	endfor;
	NL:= halfnum;
	director_invisible( true, thickfact, false );
	label.urt(btex $\vec{n}$ etex scaled 8.5, rp(tippos) );
%	label.ulft(btex $\vec{n}$ etex scaled 8.5 rotated 90, rp(tippos) );
%	tmppict := currentpicture rotated -90;
%	currentpicture := nullpicture;
%	draw tmppict;
endfig;

verbatimtex
\end{document}
etex

end;
