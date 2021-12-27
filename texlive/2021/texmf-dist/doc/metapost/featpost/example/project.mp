% project.mp
% L. Nobre G.
% 2013

verbatimtex \documentclass{article}\begin{document} etex

beginfig(1);
	ahangle := 30;
	numeric angthe, angalf, angpsi, xtopo, ytopo, xpoin;
	numeric radthe, radalf, radpsi, extenfac;
	numeric ysize, xsize, margs, units, parac, i, dim;
	pair vzer, incl, hori, topo;
	path parab, plani, veloc, xaxis, yaxis, auxil;
	pen thick, usual, vthin;
angpsi = -13;
extenfac = 1.2;
	dim = 50;
	units = 1mm;
	margs = 8units;
	xsize = 70units;
	ysize = 40units;
	radthe = 34units;
	radalf = 46units;
	radpsi = 22units;
	thick = pencircle scaled 0.5units;
	usual = pencircle scaled 0.2units;
	vthin = pencircle scaled 0.1units;
	xtopo = xsize/2-margs;
	ytopo = ysize-2*margs;
	parac = ytopo/(xtopo*xtopo);
	yaxis = ( 0, -margs/2 )--( 0, ysize-margs );
	xaxis = ( -margs, 0 )--( xsize-margs, 0 );
	plani = (margs*( -cosd( angpsi ), -sind( angpsi ) ) )
	--(extenfac*( 2*xtopo, 2*xtopo*sind( angpsi )/cosd( angpsi ) ) );
	topo = ( xtopo, ytopo );
	parab = origin;
	for i = 1 upto extenfac*dim:
	    xpoin := 2*i*xtopo/dim;
	    parab := parab..(xpoin,ytopo-parac*(xpoin-xtopo)*(xpoin-xtopo));
	endfor;
	ypart( vzer ) = ytopo;
	vzer = whatever*( direction 0 of parab );
	veloc = origin--vzer;
	hori = ( 2*xtopo, 0 );
	incl = plani intersectionpoint reverse parab;
	auxil = origin--topo;
	pickup thick;
	draw parab cutafter reverse plani;
	drawarrow veloc;
	pickup usual;
	drawarrow xaxis;
	drawarrow yaxis;
	draw plani;
	draw parab cutbefore reverse plani dashed evenly;
%	draw auxil dashed evenly;
	pickup vthin;
	draw fullcircle scaled radthe cutafter veloc;
%	draw fullcircle scaled radalf cutafter auxil;
	draw fullcircle scaled radpsi cutafter plani;

	dotlabel.urt( btex $H$ etex, hori );
	dotlabel.urt( btex $I$ etex, incl );
%	dotlabel.top( btex $T$ etex, topo );
	label.bot( btex $x$ etex, point 1 of xaxis );
	label.lft( btex $y$ etex, point 1 of yaxis );
	label.lft( btex $\vec{v}_0$ etex, vzer );
	label.urt( btex $\Theta$ etex, ( radthe/2, 0 ) );
%	label.urt( btex $\alpha$ etex, ( radalf/2, 0 ) );
	label.urt( btex $\Psi$ etex, ( radpsi/2, 0 ) );
endfig;

verbatimtex \end{document} etex

end;
