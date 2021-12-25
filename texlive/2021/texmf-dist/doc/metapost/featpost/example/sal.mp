% Cesium Cloride Crystal Structure and EBG fields together
% L. Nobre G. 
% 2002

input featpost3Dplus2D;

verbatimtex 
\documentclass{article} 
\begin{document} 
etex

beginfig(1);
    numeric sfr, radfac, ind;
    color xv, yv, zv;
    path p[];
    f := 3*(5,2,3);
    Spread := 50;
    sfr := 0.2;
    radfac := 2.0;
    xv := (sfr,0,0);
    yv := (0,sfr,0);
    zv := (0,0,sfr);
    V1 := (-1,-1,-1);
    V2 := (-1,1,-1);
    V3 := (1,-1,-1);
    V4 := (1,1,-1);
    V5 := (-1,-1,1);
    V6 := (-1,1,1);
    V7 := (1,-1,1);
    V8 := (1,1,1);
    V9 := (0,0,0);
    for ind=1 upto 8:
	p[ind] := rigorousfearpath( V[ind], sfr );
    endfor;
    p9 := rigorousfearpath( V9, sfr*radfac );
    unfill p1;
    draw p1;
    draw rp(V1+yv)--rp(V2-yv);
    draw rp(V1+xv)--rp(V3-xv);
    unfill p2;
    draw p2;
    unfill p3;
    draw p3;
    draw rp(V4-yv)--rp(V3+yv);
    draw rp(V4-xv)--rp(V2+xv);
    unfill p4;
    draw p4;
    draw rp(V1+zv)--rp(V5-zv);
    draw rp(V2+zv)--rp(V6-zv);
    unfill p9;
    draw p9;
    draw rp(V3+zv)--rp(V7-zv);
    draw rp(V4+zv)--rp(V8-zv);
    unfill p5;
    draw p5;
    draw rp(V5+yv)--rp(V6-yv);
    draw rp(V5+xv)--rp(V7-xv);
    unfill p6;
    draw p6;
    unfill p7;
    draw p7;
    draw rp(V8-yv)--rp(V7+yv);
    draw rp(V8-xv)--rp(V6+xv);
    unfill p8;
    draw p8;
    label( btex Cs${}^{+}$ etex, rp( V8 ) );
    label( btex Cl${}^{-}$ etex, rp( V9 ) );
endfig;

beginfig(2);
	numeric inc, ind, tep, larg, altr, lene;
	path efld, velv, arch, auxi;
	picture bfld;
	pair vpos, cent;
	tep = 1cm;
	larg = 4;
	altr = 3;
	lene = larg*tep;
	bfld = btex $\times$ etex;
%	bfld := bfld shifted -(center bfld);
	efld = origin--(lene,0);
	for ind = 0 upto altr:
	    drawarrow efld shifted (0,ind*tep);
	endfor;
	auxi = efld shifted (0,tep);
	for ind = 1 upto altr:
	    for inc = 1 step 2 until 2*larg:
		label( bfld, (0.5*inc*tep,(ind-0.5)*tep) );
	    endfor;
	endfor;
	vpos = (tep,0.5tep);
	label.lft( btex $\vec{B}$ etex, (0.5tep,(altr-0.5)*tep) );
	label.lft( btex $\vec{E}$ etex, origin );
	label.top( btex $\vec{v}$ etex, vpos );
	velv = vpos--((larg-1)*tep,(altr-0.5)*tep);
	drawarrow velv;
	cent = auxi intersectionpoint velv;
	arch = ((fullcircle scaled (altr*tep*0.87)) shifted cent) cutafter velv;
        draw arch;
	label.rt( btex $\theta$ etex, point 0.5*length arch of arch );
endfig;

verbatimtex \end{document} etex

end;

