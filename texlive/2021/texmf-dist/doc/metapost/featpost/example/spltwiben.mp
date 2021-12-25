% spltwiben.mp
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

def zsurface( expr xc, yc ) = 
	 ((3*xc-5)*(3*xc+5)/(3*sqrt(3)*yc-5))
enddef;

f := (4,-5,3);
Spread := 120;

beginfig(1); 
	numeric i, j, nx, ny, lowx, higx, lowy, higy; 
	numeric posx, posy, posz, facz;
	color gridx[], pos, reforig, endbeta, endpsip;
	color endfiel, kdiv, krot, kben, akdiv, akrot, akben;
	path facex;
	facz = 7;
	nx = 15;
	ny = 50;
	lowx = 0;
	higx = 0.5*sqrt(3);
	reforig = (0,-0.5,0);
	endpsip = (0.3,-0.5,0);
	endbeta = (0,-0.2,0);
	endfiel = (0,-0.5,9/facz);
	   kdiv = (0,-0.5,8/facz);
	   krot = (0,-0.5,2/facz);
	   kben = (0,-0.5,5/facz);
	  akdiv = (higx,0.5,8/facz);
	  akrot = (higx,-0.5,2/facz);
	  akben = (0,0,5/facz);
	draw rp((0,0,0))--rp((0,0,5/facz));
	for i=1 upto nx:
	    posx := lowx + i*( higx - lowx )/nx;
       	    lowy := -sqrt(3)*posx/3;
	    higy := -lowy;
	    for j=0 upto ny:
		posy := lowy + j*( higy - lowy )/ny;
		posz := zsurface( posx, posy )/facz;
		gridx[j] := ( posx, posy, posz );
		% show gridx[j];
	    endfor;
	    pos := ( posx, lowy, 0 );
	    facex := rp( pos )--rp( gridx0 );
	    for j=1 upto ny:
		facex := facex--rp( gridx[j] );
	    endfor;
	    facex:=facex--rp((posx,higy,0))--cycle;
	    unfill facex;
	    draw facex;
	endfor;
	drawarrow rp(reforig)--rp(endpsip);
	drawarrow rp(reforig)--rp(endbeta);
	drawarrow rp(reforig)--rp(endfiel);
	label.ulft(btex $\mathcal{X}$ etex, rp(endbeta) );
	label.llft(btex $\mathcal{Z}$ etex, rp(endpsip) );
	label.llft(btex $\mathcal{H}$ etex, rp(endfiel) );
	dotlabel("",rp(reforig+(0,0,Z(gridx[ny]))));
	dotlabel("",rp(reforig+(0,0,Z(gridx0))));
	dotlabel("",rp(kben));
	dotlabel.rt(btex $\mathcal{K}_1$ etex,rp(gridx[ny]));
	dotlabel.lrt(btex $\mathcal{K}_2$ etex,rp(gridx0));
	dotlabel.lft(btex $\mathcal{K}_3$ etex,rp(akben));
endfig;

verbatimtex 
\end{document} 
etex

end;
