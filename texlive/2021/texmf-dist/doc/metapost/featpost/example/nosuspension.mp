% nosuspension.mp
% the hard way to write

beginfig(0);
	picture l[], auxpic;
	numeric gridsp, u, pensz, yprop, 
	        altura, largura, small, big, thisjump, thatjump;
	pen thispen, thatpen;
	path auxpat, thispat, thatpat;
	pair currpos, reforigin;
	u = 1mm;
	pensz = 0.5u;
	yprop = 1.5;
	gridsp = 6u;
	altura = 3gridsp;
	largura = 2gridsp;
	small = 0.5gridsp;
	big = 0.6667gridsp;
	thisjump = 3gridsp;
	thatjump = 4gridsp;
	reforigin = (10cm,15cm) + left*9.5cm;
	auxpat = (0,yprop*u)--(0.4u,yprop*u){right}..{down}(1u,0.9u)--(1u,0);
	auxpat := auxpat--(auxpat reflectedabout(origin,right));
	thispat = auxpat--(auxpat reflectedabout(origin,up));
	thatpat = thispat rotated 90;
	thispen = makepen(thispat) scaled pensz;
	thatpen = makepen(thatpat) scaled pensz;
	draw (big,0)--origin
		      --(0,altura)
		      --(largura,0)
		      --(largura,altura)
		      --(largura-big,altura) withpen thatpen;
	l1 = currentpicture;
	currentpicture := nullpicture;
	draw (largura,largura)--(0,largura)
		      --origin
		      --(altura,largura)
		      --(altura,0)
		      --(altura-largura,0) withpen thatpen;
	l2 = currentpicture;
	currentpicture := nullpicture;
	draw (altura,largura-small)--(altura,largura)
		      --(0,largura)
		      --(altura,0)
		      --origin
		      --(0,small) withpen thispen;
	l3 = currentpicture;
	currentpicture := nullpicture;
	draw (small,altura)--(0,altura)
		      --origin
		      --(largura,0)
		      --(largura,altura)
		      --(largura-small,altura) withpen thispen;
	l4 = currentpicture;
	currentpicture := nullpicture;
	draw (small,0)--origin
		      --(0,altura)
		      --(largura,altura)
		      --(largura,altura-gridsp)
		      --(0,altura-gridsp) withpen thispen;
	l5 = currentpicture;
	currentpicture := nullpicture;
	draw (largura,small)--(largura,0)
		      --origin
		      --(0,altura)
		      --(largura,altura)
		      --(largura,altura-small) withpen thispen;
	draw (0,largura)--(largura-small,largura) withpen thispen;
	l6 = currentpicture;
	currentpicture := nullpicture;
	draw (-small,0)--origin
		      --(0,altura)
		      --(-small,altura) withpen thispen;
	l7 = currentpicture;
	currentpicture := nullpicture;
	draw (largura-small,largura)--(largura-small,altura)
		      --(0,altura)
		      --origin
		      --(largura,0)
		      --(largura,largura)
		      --(0,largura) withpen thispen;
	l8 = currentpicture;
	currentpicture := nullpicture;
	draw origin--(small,0)
		      --(small,largura)
		      --(0,largura) withpen thatpen;
	draw (small,largura-big)--(largura,largura)
		      --(largura,largura-small) withpen thatpen;
	draw (altura,small)--(altura,0)
		      --(largura-big,largura-small) withpen thatpen;
	l9 = currentpicture;
	currentpicture := nullpicture;
	currpos = (thisjump,thatjump+gridsp)+reforigin;
	draw l1 scaled 2 shifted currpos;
	currpos := currpos+(2thisjump,0);
	draw l2 scaled 2 shifted currpos withcolor 0.875red;
	currpos := reforigin;
	draw l3 shifted currpos;
	currpos := currpos+(thatjump,0);
	draw l4 shifted currpos;
	currpos := currpos+(thisjump,gridsp);
	draw l3 shifted currpos;
	currpos := currpos+(thatjump,-gridsp);
	draw l5 shifted currpos;
	currpos := currpos+(thisjump,0);
	draw l6 shifted currpos;
	currpos := currpos+(thisjump,gridsp);
	draw l1 shifted currpos;
	currpos := currpos+(thisjump,-gridsp);
	draw l3 shifted currpos;
	currpos := currpos+(thatjump+small,0);
	draw l7 shifted currpos;
	currpos := currpos+(gridsp,gridsp);
	draw l2 shifted currpos;
	currpos := currpos+(thatjump,-gridsp);
	draw l1 shifted currpos;
	currpos := reforigin+(2thatjump,-thatjump);
	draw l8 shifted currpos;
	currpos := currpos+(thisjump+small,-gridsp);
	draw l7 shifted currpos;
	currpos := currpos+(gridsp,gridsp);
	draw l9 shifted currpos;
	currpos := currpos+(thatjump,-gridsp);
	draw l6 shifted currpos;
	currpos := currpos+(thisjump,gridsp);
	draw l3 shifted currpos;
endfig;

end;


