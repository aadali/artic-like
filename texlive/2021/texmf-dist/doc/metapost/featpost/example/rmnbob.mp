% rmnbob.mp
% L. Nobre G.
% probably before 2004

input featpost3Dplus2D;

Spread := 70;
f := (5,7,10);

beginfig(1);
	numeric i, r[], stepr, len, Numr, Minn, Maxn;
	numeric fws, NumLines, linex;
	path aro[], oneline;
	color centre, axe, baseone, basetwo;
	
	stepr := 0.1;
	len := 3;
	Numr := 5;
	Minn := 5;
	centre := (0,0,0);
	axe := (1,0,0);
	fws := 2;
	NumLines := 30;
	for i=0 upto NumLines:
	    linex := fws*(i/NumLines-0.5);
	    oneline := rp((linex,0.5*fws,-0.5*fws))
                     --rp((linex,-0.5*fws,-0.5*fws))
                     --rp((linex,-0.5*fws,0.5*fws));
	    draw oneline withcolor (0.7*background);
	endfor;
	baseone := centre+0.5*len*axe;
	basetwo := centre-0.5*len*axe;
	Maxn := Minn+Numr;
	for i=Minn upto Maxn:
	    r[i] := i*stepr;
	    aro[i] := rigorouscircle(baseone,axe,r[i]);
	endfor;
	rigorousdisc( 0, true, baseone, r[Maxn], -len*axe );
	drawoptions( dashed evenly );
	rigorousdisc( 0, true, baseone, r[Minn], -len*axe );
	for i=Maxn downto Minn:
	    filldraw aro[i] withcolor ((0.5+0.3*normaldeviate)*background);
	endfor;
	drawoptions( );
	draw aro[Maxn];
endfig;

end.
