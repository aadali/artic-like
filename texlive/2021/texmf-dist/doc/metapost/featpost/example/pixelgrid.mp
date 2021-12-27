% pixelgrid.mp
% L. Nobre G. 
% 2001

%input mp-tool;
path fullsquare;
fullsquare := unitsquare shifted - center unitsquare ;

beginfig(1); 						
	numeric u, ix, iy, s, h, v, i, j, n, f, r, l;
	numeric last, angl, hori;
	path frame;
	pair cpos;
	u = 2mm;
	s = 1.2u;
	n = 1;		% 6
	f = 30;		% 5
	h = f*n;
	v = ceiling(h/sqrt(2));
	frame = fullsquare scaled 1u;
%	for ix=0 upto h:
%	    for iy=-v upto v:
%		draw frame shifted (ix*s,iy*s);
%	    endfor;
%	endfor;
	filldraw frame;
	for i=1 upto n:
	    r := i*f;
	    l := round(r/sqrt(2));
	    last := 0;
	    for j=0 upto l:
		hori := round(sqrt(r*r-j*j));
		cpos := (hori*s,j*s);
		angl := angle(cpos);
		filldraw frame shifted cpos;
		show angl-last;
		last := angl;
%		draw origin--cpos;
	    endfor;
	endfor;
%	draw fullcircle scaled (2*r*s) withcolor red;
%	draw boundingbox currentpicture;
endfig;

end;
