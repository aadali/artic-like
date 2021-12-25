% some2Dvecs.mp
% L. Nobre G.
% 2002

verbatimtex
\documentclass{article}
\begin{document}
etex

beginfig(1);
	numeric hsize, qsize, nsize, ldist;
	numeric alray, beray, alang, beang;
	numeric stampnum, initang, i, stepang, stampoffset;
	pair hvec, qvec, nvec;
	path Hvec, Qvec, Nvec, alarc, bearc, vertic;
	picture labelstr, labelpic;
	color mycolor;
	picture stamp;
	bboxmargin := 4mm;
	stampnum = 36;
	stampoffset = 0.42mm;
	alang = 6;
	beang = 17;
	hsize = 6mm;
	qsize = 3*hsize;
	nsize = 4*hsize;
	ldist = 3.17*hsize;
	alray = 2*0.735*nsize;
	beray = 2*0.6*qsize;
	mycolor = (0.8,0.7,1);
%	drawoptions( withcolor mycolor );
	hvec = hsize*left;
	qvec = qsize*(-sind(beang),cosd(beang));
	nvec = nsize*(-sind(alang),cosd(alang));
	Hvec = origin--hvec;
	Qvec = origin--qvec;
	Nvec = origin--nvec;
	vertic = origin--ldist*up;
	alarc = fullcircle scaled alray 
		cutbefore vertic cutafter Nvec;
	bearc = fullcircle scaled beray 
		cutbefore vertic cutafter Qvec;
	drawarrow Hvec;
	drawarrow Qvec;
	drawarrow Nvec;
	draw alarc;
	draw bearc;
	draw vertic dashed evenly;
	draw unitsquare scaled (0.2*hsize);
	labelstr = btex 
		   $2,\!4\times 1,\!8\,\mbox{mm}^2$
		   etex;
	label.rt( labelstr, ldist*right );
%	label.lft( btex $\vec{n}_0$ etex, nvec );
	label.lft( btex $x$ etex, nvec );
	label.lft( btex $\vec{q}$ etex, qvec );
	label.top( btex $\vec{B}$ etex, hvec );
	label.urt( btex $-0.48^\circ$ etex, point 1 of alarc );
	label.lft( btex $6^\circ$ etex, point 1 of bearc );
	stamp = currentpicture;
	fill bbox stamp withcolor (0.85,0.85,0);
	initang = 180/stampnum;
	stepang = 2*initang;
	drawoptions();
	for i=0 upto (stampnum-1):
	    draw stamp 
		 shifted (stampoffset*dir(initang+i*stepang))
		 withcolor white;
	endfor;
	draw stamp; 
endfig;

verbatimtex
\end{document}
etex

end.
