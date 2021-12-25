% anglinerigorouscircle.mp
% L. Nobre G. 
% 2006

input featpost3Dplus2D;

% bash: export TEX=latex
% tcsh: setenv TEX latex

%prologues := 3;
prologues := 0;

verbatimtex
\documentclass{article} 
\usepackage{beton,concmath,ccfonts} 
\begin{document} 
etex

beginfig(1);
    color orig, axxc, ayyc, azzc;
    color field, nalfa, ngamma, alfacnt, gammacnt, nazy, nax;
    path alfa, gamma;
    numeric nalfad, theta, radius, gammang, xxsize, yysize, zzsize;
    numeric xxpos, yypos, zzpos, fieldpos, deltarad, deltapos;
    numeric thetarad, thetapos, nalfapos, ngammapos;
    numeric etalpharad, etalphapos, alphapos, gammapos;
    picture xxlabel, yylabel, zzlabel, fieldlabel, deltalabel;
    picture nalfalabel, ngammalabel, etalphalabel, alphalabel;
    picture gammalabel, thetalabel;
    
    f 			:= (3,4,2);
    Spread 		:= 100;
    xxsize 		= 0.8;
    xxpos 		= 6;
    yysize 		= 0.8;
    yypos 		= 6;
    zzsize 		= 1.3;
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

    xxlabel 	= btex $x_c$ etex;
    yylabel 	= btex $y_c$ etex;
    zzlabel 	= btex $z_c$ etex;
    deltalabel 	= btex $\Delta$ etex;
    fieldlabel 	= btex $\mathbf{B}$ etex;
    thetalabel 	= btex $\theta$ etex;
    nalfalabel 	= btex $\vec{n}_\alpha$ etex;
    ngammalabel = btex $\vec{n}_\gamma$ etex;
    etalphalabel= btex $\eta_\alpha$ etex;
    alphalabel 	= btex $\alpha$ etex;
    gammalabel 	= btex $\gamma$ etex;
    orig 		= (0,0,0);
    axxc 		= (xxsize,0,0);
    ayyc 		= (0,yysize,0);
    azzc 		= (0,0,zzsize);
    nalfa 		= (0,sind(theta),cosd(theta));
    alfacnt 		= (0,0,nalfad*cosd(theta));
    radius 		= nalfad*sind(theta);
    alfa 		= rigorouscircle( alfacnt, (0,0,1), radius );
    gammacnt 		= cosd(theta)*cosd(theta)*nalfa;
    radius 		:= cosd(theta)*sind(theta);
    gamma 		= rigorouscircle( gammacnt, nalfa, radius );
    nazy 		= radius*sind(gammang)*ncrossprod( nalfa, axxc );
    nax 		= radius*cosd(gammang)*ncrossprod( ayyc, nalfa );
    ngamma 		= gammacnt+nax+nazy;
    draw alfa dashed evenly;
    draw gamma dashed evenly;
    drawarrow rp(orig)..rp(axxc);
    drawarrow rp(orig)..rp(ayyc);
    drawarrow rp(orig)..rp(azzc);
    drawarrow rp(orig)..rp(field);
    drawarrow rp(orig)..rp(nalfa);
    drawarrow rp(orig)..rp(ngamma);
    drawarrow subpath (4,3) of alfa;
    drawarrow subpath (8,7) of gamma;
    anglinen( field, azzc, orig, deltarad, deltalabel, deltapos );	
    anglinen( azzc, nalfa, orig, thetarad, thetalabel, thetapos );	
    anglinen( ngamma, field, orig, etalpharad, etalphalabel, etalphapos );
    labeln( xxlabel, rp(axxc), xxpos );
    labeln( yylabel, rp(ayyc), yypos );
    labeln( zzlabel, rp(azzc), zzpos );
    labeln( fieldlabel, rp(field), fieldpos );
    labeln( nalfalabel, rp(nalfa), nalfapos );
    labeln( ngammalabel, rp(ngamma), ngammapos );
    labeln( alphalabel, point 3 of alfa, alphapos );
    labeln( gammalabel, point 7 of gamma, gammapos );
    pickup pencircle scaled 1.5;
    draw rp(alfacnt);
    draw rp(gammacnt);
%   produce_auto_scale;
endfig;

verbatimtex 
\end{document} 
etex

end;
