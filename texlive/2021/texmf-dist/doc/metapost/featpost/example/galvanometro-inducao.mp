% galvanometro-inducao.mp
% Pedro Sebastiao
% 2008

input featpost3Dplus2D;

% bash: export TEX=latex
% tcsh: setenv TEX latex

verbatimtex
\documentclass{article} 
\usepackage{beton,concmath,ccfonts} 
\begin{document} 
etex

beginfig(1);
    color orig, axxc, ayyc, azzc,hum,hdois,pnove,pdez;
    color field, nalfa, ngamma, alfacnt, gammacnt, nazy, nax,pum,pdois,ptres,pquatro,pcinco,pseis,psete,poito,zpos,zppos,gammacnt;
    path alfa, gamma, circ;
    numeric nalfad, theta, radius, gammang, xxsize, yysize, zzsize, phi,teta;
    numeric xxpos, yypos, zzpos, fieldpos, deltarad, deltapos,inc,resize,sum,sdois;
    numeric thetarad, thetapos, nalfapos, ngammapos,aux,alpha,ia,iz;
    numeric etalpharad, etalphapos, alphapos, gammapos;
    picture xxlabel, yylabel, zzlabel, fieldlabel, alphalabel;
    picture rlabel, ngammalabel, etalphalabel, philabel;
    picture gammalabel, thetalabel;
    
    f 			:= (2,2.5,2.0);
 %   f 			:= (2,2,20.0);
    Spread 		:= 100;
    xxsize 		= 1.1;
    xxpos 		= 6;
    yysize 		= 1.1;
    yypos 		= 6;
    zzsize 		= 1;
    zzpos 		= 4;
    field 		= (0.7,0,0.3);
    fieldpos 		= 2;
    deltarad 		= 0.45;
    deltapos 		= 3;
    nalfad 		= 1.0;
    theta 		= 25;
    thetarad 		= 0.3;
    thetapos 		= 2;
    nalfapos 		= 7;
    gammang 		= 65;
    ngammapos 		= 7;
    etalpharad 		= 0.4;
    etalphapos 		= 2;
    alphapos 		= 6;
    gammapos 		= 0;
    radius              =0.5;
    alpha := 60;
    
    xxlabel 	= btex $x$ etex;
    yylabel 	= btex $y$ etex;
    zzlabel 	= btex $z$ etex;
    alphalabel 	= btex $2\alpha$ etex;
    fieldlabel 	= btex $\mathbf{B}$ etex;
    thetalabel 	= btex $\theta$ etex;
    rlabel 	= btex $r$ etex;
    ngammalabel = btex $\vec{n}_\gamma$ etex;
    etalphalabel= btex $\eta_\alpha$ etex;
    philabel 	= btex $\phi$ etex;
    gammalabel 	= btex $\gamma$ etex;
    orig 		= (0,0,0);
    axxc 		= (xxsize,0,0);
    ayyc 		= (0,yysize,0);
    azzc 		= (0,0,zzsize);
    hum = 0.9*azzc;
    hdois = -0.9*azzc;
    resize:=0.85;
    

    pum := (resize*cosd(alpha),resize*sind(alpha),0.9);
    pdois := (-resize*cosd(alpha),resize*sind(alpha),0.9);
    ptres :=(resize*cosd(alpha),resize*sind(alpha),-0.9);
    pquatro :=(-resize*cosd(alpha),resize*sind(alpha),-0.9);
    draw rp(pum)..rp(ptres);
    draw rp(pdois)..rp(pquatro);
    pcinco := (resize*cosd(alpha),resize*sind(-alpha),0.9);
    pseis := (-resize*cosd(alpha),resize*sind(-alpha),0.9);
    psete :=(resize*cosd(alpha),resize*sind(-alpha),-0.9);
    poito :=(-resize*cosd(alpha),resize*sind(-alpha),-0.9);
    draw rp(pcinco)..rp(psete);
    draw rp(pseis)..rp(poito);

    anglinen( pum, pdois, hum, resize, alphalabel, 4 );	
    anglinen( ptres, pquatro, -hum, resize, "", deltapos );	
    anglinen( pcinco, pseis, hum, resize, "", deltapos );	
    anglinen( psete, poito, -hum, resize, "", deltapos );	

    draw rp(hum)..rp(pdois) dashed evenly;
   draw rp(hum)..rp(pcinco) dashed evenly;
    
    alpha:=60;
    resize:=0.80;

    pickup pencircle scaled 3;
    
    pum := (resize*cosd(alpha),resize*sind(alpha),0.9);
    pdois := (-resize*cosd(alpha),-resize*sind(alpha),0.9);
    ptres :=(resize*cosd(alpha+180),-resize*sind(alpha),-0.9);
    pquatro :=(resize*cosd(alpha),resize*sind(alpha),-0.9);
    draw rp(pum)..rp(pdois);
    draw rp(pdois)..rp(ptres);
    draw rp(ptres)..rp(pquatro);
    draw rp(pquatro)..rp(pum);

    pickup pencircle scaled 1.5;
    resize:=0.85;

    for ia=60 step 10 until (120):
      for iz=-0.0 step 1.4 until (0.71):
	pum := (resize*cosd(ia),resize*sind(ia),iz);
	pdois := (cosd(ia),sind(ia), iz);
	ptres := (resize*cosd(ia+180),resize*sind(ia+180),iz);
	pquatro := (cosd(ia+180),sind(ia+180), iz);
	drawarrow rp(pum)..rp(pdois) withcolor red;
	draw rp(pdois)..rp(pdois+(0,yysize,0)) withcolor red;
	drawarrow rp(pquatro)..rp(ptres) withcolor red;
	draw rp(pquatro-(0,yysize,0))..rp(pquatro) withcolor red;
	hdois:=pum;
	pcinco := (0.5*resize*cosd(ia),0.5*resize*sind(ia),iz);
	pseis := (0.5*resize*cosd(360-ia),0.5*resize*sind(360-ia),iz);
	pdez := (resize*cosd(360-ia),resize*sind(360-ia),iz);
	for i=0 step 0.01 until 1:
	  poito := ((1-i)**3)*pum+3*i*((1-i)**2)*pcinco+3*(i**2)*(1-i)*pseis+(i**3)*pdez;
	  draw rp(hdois)..rp(poito)  withcolor red;
	  hdois:=poito;
	endfor
      endfor
    endfor

    for ia=40 step 10 until (50):
      for iz=-0.0 step 1.4 until (0.71):
	pum := (resize*cosd(ia),resize*sind(ia),iz);
	pdois := (cosd(ia),sind(ia), iz);
	ptres := (resize*cosd(ia+180),resize*sind(ia+180),iz);
	pquatro := (cosd(ia+180),sind(ia+180), iz);
	drawarrow rp(pum)..rp(pdois) withcolor black;
	draw rp(pdois)..rp(pdois+(0,yysize,0)) withcolor black;
	drawarrow rp(pquatro)..rp(ptres) withcolor black;
	draw rp(pquatro-(0,yysize,0))..rp(pquatro) withcolor black;
	hdois:=pum;
	pcinco := (0.5*resize*cosd(ia),0.5*resize*sind(ia),iz);
	pseis := (0.5*resize*cosd(360-ia),0.5*resize*sind(360-ia),iz);
	pdez := (resize*cosd(360-ia),resize*sind(360-ia),iz);
	for i=0 step 0.01 until 1:
	  poito := ((1-i)**3)*pum+3*i*((1-i)**2)*pcinco+3*(i**2)*(1-i)*pseis+(i**3)*pdez;
	  draw rp(hdois)..rp(poito)  withcolor black;
	  hdois:=poito;
	endfor
      endfor
    endfor
    for ia=130 step 10 until (140):
      for iz=-0.0 step 1.4 until (0.71):
	pum := (resize*cosd(ia),resize*sind(ia),iz);
	pdois := (cosd(ia),sind(ia), iz);
	ptres := (resize*cosd(ia+180),resize*sind(ia+180),iz);
	pquatro := (cosd(ia+180),sind(ia+180), iz);
	drawarrow rp(pum)..rp(pdois) withcolor black;
	draw rp(pdois)..rp(pdois+(0,yysize,0)) withcolor black;
	drawarrow rp(pquatro)..rp(ptres) withcolor black;
	draw rp(pquatro-(0,yysize,0))..rp(pquatro) withcolor black;
	hdois:=pum;
	pcinco := (0.5*resize*cosd(ia),0.5*resize*sind(ia),iz);
	pseis := (0.5*resize*cosd(360-ia),0.5*resize*sind(360-ia),iz);
	pdez := (resize*cosd(360-ia),resize*sind(360-ia),iz);
	for i=0 step 0.01 until 1:
	  poito := ((1-i)**3)*pum+3*i*((1-i)**2)*pcinco+3*(i**2)*(1-i)*pseis+(i**3)*pdez;
	  draw rp(hdois)..rp(poito)  withcolor black;
	  hdois:=poito;
	endfor
      endfor
    endfor

    
    ia=60;
    iz = 0.7;
    pnove := (resize*cosd(ia),resize*sind(ia),iz);
    pdez := (resize*cosd(360-ia),resize*sind(360-ia),iz);
    pcinco := (0.5*resize*cosd(ia),0.5*resize*sind(ia),iz);
    pseis := (0.5*resize*cosd(360-ia),0.5*resize*sind(360-ia),iz);
    pickup pencircle scaled 5.5;
%    draw rp(pdois);
%    draw rp(pquatro);
%    draw rp(pnove);
%    draw rp(pdez);
    
    pickup pencircle scaled 0.5;
	
    draw rigorouscircle( (0,0,0), (0,0,1), 1) dashed evenly;
    resize:=1;
    alpha := 60;

    iz:=0.0;
    pum := (resize*cosd(alpha),resize*sind(alpha),iz);
    pdois := (-resize*cosd(alpha),resize*sind(alpha),iz);
    ptres :=(resize*cosd(alpha),resize*sind(alpha),-iz);
    pquatro :=(-resize*cosd(alpha),resize*sind(alpha),-iz);
    pcinco := (resize*cosd(alpha),resize*sind(-alpha),iz);
    pseis := (-resize*cosd(alpha),resize*sind(-alpha),iz);
    psete :=(resize*cosd(alpha),resize*sind(-alpha),-iz);
    poito :=(-resize*cosd(alpha),resize*sind(-alpha),-iz);
    
    anglinen( pum, pdois, hum-(0.9-iz)*azzc, resize, "", deltapos );	
    anglinen( ptres, pquatro, -hum+(0.9-iz)*azzc, resize, "", deltapos );	
    anglinen( pcinco, pseis, hum-(0.9-iz)*azzc, resize, "", deltapos );	
    anglinen( psete, poito, -hum+(0.9-iz)*azzc, resize, "", deltapos );	

    resize:=0.85;
    alpha := 60;
    pickup pencircle scaled 0.5;
    iz:=0.0;
    pum := (resize*cosd(alpha),resize*sind(alpha),iz);
    pdois := (-resize*cosd(alpha),resize*sind(alpha),iz);
    ptres :=(resize*cosd(alpha),resize*sind(alpha),-iz);
    pquatro :=(-resize*cosd(alpha),resize*sind(alpha),-iz);
    pcinco := (resize*cosd(alpha),resize*sind(-alpha),iz);
    pseis := (-resize*cosd(alpha),resize*sind(-alpha),iz);
    psete :=(resize*cosd(alpha),resize*sind(-alpha),-iz);
    poito :=(-resize*cosd(alpha),resize*sind(-alpha),-iz);
    
    anglinen( pum, pdois, hum-(0.9-iz)*azzc, resize, "", deltapos );	
    anglinen( ptres, pquatro, -hum+(0.9-iz)*azzc, resize, "", deltapos );	
    anglinen( pcinco, pseis, hum-(0.9-iz)*azzc, resize, "", deltapos );	
    anglinen( psete, poito, -hum+(0.9-iz)*azzc, resize, "", deltapos );	

    
    pickup pencircle scaled 1;
    drawarrow rp(-azzc)..rp(1.3*azzc);
%    drawarrow rp(-axxc)..rp(axxc);
%    drawarrow rp(-ayyc)..rp(ayyc);
    pickup pencircle scaled 1.5;
%    draw rp(alfacnt);
%    draw rp(gammacnt);
%   produce_auto_scale;
endfig;

verbatimtex 
\end{document} 
etex

end;
