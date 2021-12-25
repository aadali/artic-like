% Crossed Electric and Magnetic fields
% L. Nobre G. 
% 2002

%        defaultfont := "putr8r";
%        defaultscale := 6.5;

verbatimtex 
\documentclass{article} 
\usepackage{palatino} 
\usepackage{concmath} 
\usepackage{ccfonts}
\begin{document} 
etex

beginfig(1);
	numeric inc, ind, tep, larg, altr, lene;
	path efld, velv, arch, auxi;
	picture bfld;
	pair vpos, cent;
	tep = 1cm;
	larg = 4;
	altr = 3;
	lene = larg*tep;
	bfld = btex $\times$ etex;
%	bfld = thelabel("x",origin);
	bfld := bfld shifted -(center bfld);
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
%	label.lft( "B", (0.5tep,(altr-0.5)*tep) );
%	label.lft( "E", origin );
%	label.top( "v", vpos );
	velv = vpos--((larg-1)*tep,(altr-0.5)*tep);
	drawarrow velv;
	cent = auxi intersectionpoint velv;
	arch = ((fullcircle scaled (altr*tep*0.87)) shifted cent) cutafter velv;
        draw arch;
	label.rt( btex $\theta$ etex, point 0.5*length arch of arch );
%	label.rt( "T", point 0.5*length arch of arch );
endfig;

verbatimtex \end{document} etex

end;

