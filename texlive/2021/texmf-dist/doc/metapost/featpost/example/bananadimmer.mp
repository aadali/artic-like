% bananadimmer.mp
% Pedro J. Sebastiao.
% L. Nobre G. 
% 2006

input featpost3Dplus2D;

verbatimtex
\documentclass{article}
\usepackage{beton}
\usepackage{concmath}
\usepackage{ccfonts}
\begin{document}
etex

beginfig(1);
    color axxc, ayyc, azzc;
    color nalfa, ngamma;
    numeric radius, gammang, xxsize, yysize, zzsize;
    numeric xxpos, yypos, zzpos;
    numeric angulo;
    numeric nalfapos;
    picture xxlabel, yylabel, zzlabel, nalfalabel;
    
    f    	       := (3,4,2);
    Spread             := 60;

    gammang 		= 0;
    angulo              = -90;
    beta                = 11.8;
    nbeta               = 1.3;
    delta               = 35;
    nalfapos 		= 1;
    xxsize 		= 0.8;
    xxpos 		= 6;
    yysize 		= 1.3;
    yypos 		= 6;
    zzsize 		= 1.3;
    zzpos 		= 4;
    
    xxlabel 	= btex $\mbox{\protect \scriptsize light}$ etex;
    yylabel 	= btex $\mbox{\protect \scriptsize analyser}$ etex;
    zzlabel 	= btex $\mbox{\protect \scriptsize polariser}$ etex;
    nalfalabel 	= btex $\mathbf{n}_a$ etex;

    axxc        = (xxsize,0,0);
    ayyc 	= (0,yysize,0);
    azzc 	= (0,0,zzsize);
    nalfa       = (0,nbeta*cosd(beta)*sind(delta),nbeta*cosd(beta)*cosd(delta));
    radius 	= nbeta*sind(beta);
    
    ngamma     := 0.5*(nalfa+(radius*cosd(angulo),radius*sind(angulo),0));
    drawarrow rp(black)..rp(axxc);
    drawarrow rp(black)..rp(ayyc);
    drawarrow rp(black)..rp(azzc);
    drawarrow rp(black)..rp(nalfa);
    labeln( xxlabel, rp(axxc), xxpos );
    labeln( yylabel, rp(ayyc), yypos );
    labeln( zzlabel, rp(azzc), zzpos );
    labeln( nalfalabel, rp(nalfa), nalfapos );

    banana(
      ngamma,
      conorm(ngamma),
      ( 90, -26, 90 ),
      0.2,
      150);
    banana(
      ngamma,
      conorm(ngamma),
      ( 0, 0, 66 ),
      0.2,
      135);
endfig;

verbatimtex \end{document} etex

end.
