% optest.mp
% L. Nobre G. 
% 2000

beginfig(1);
%	linecap := squared;
	pickup pencircle scaled 0.1mm;
	numeric salto, i, shif;
	pair cen;
	shif = 2cm;
	cen = (10.5cm,14.845cm);
	draw ((-shif,0)--(shif,0)) shifted cen;
	draw ((0,-shif)--(0,shif)) shifted cen;
	pickup pencircle scaled 0.05mm;
	salto = 1mm;
	for i = 20 upto 22:
		draw fullcircle scaled (i*salto) shifted cen;
	endfor;
	for i = 27 upto 28:
		draw unitsquare shifted (-0.5,-0.5) 
				scaled (i*salto) 
				rotated 45 shifted cen;
	endfor;
	for i = 29 upto 31:
		draw fullcircle scaled (i*salto) shifted cen;
	endfor;
	for i = -10 step 2 until 10:
		draw (origin--shif*(1-0.0125*abs(i))*dir(i)) shifted cen;
	endfor;
endfig;

beginfig(2);
	linecap := squared;
	pickup pensquare scaled 0.15mm;
	numeric salto, ix, iy;
	path face;
	picture descrit;
	salto = 0.3in;
	face = unitsquare xscaled 1.33 scaled 2;
	fill face withcolor red;
	draw face;
	descrit = currentpicture;
	currentpicture := nullpicture;
	for iy = 1 upto 5:
	    for ix = 1 upto 5:
		draw descrit shifted ( ix*salto, iy*salto );
	    endfor;
	endfor;
endfig;

end;
