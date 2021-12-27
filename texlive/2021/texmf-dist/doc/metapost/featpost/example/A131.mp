% A131.mp
% P. J. Sebastiao
% 2012

input featpost3Dplus2D;

prologues:=3;

verbatimtex
  \documentclass{article}
  \usepackage{amsmath,amssymb}
  \begin{document}
etex

beginfig(1);
    color orig, axxc, ayyc, azzc,ngorig, ngorigb;
    color field, nalfa, ngamma[], alfacnt, gammacnt, nazy, nax,nalfaz, alfacntb, gammacntb,ngammab[],nazyb, naxb;
    path alfa, gamma, alfab, gammab;
    numeric nalfad, theta, radius, gradius, gammang, xxsize, yysize, zzsize, thetb, radiusb;
    numeric xxpos, yypos, zzpos, fieldpos, deltarad, deltapos;
    numeric thetarad, thetapos, nalfapos, ngammapos,veccpos;
    numeric etalpharad, etalphapos, alphapos, gammapos;
    picture xxlabel, yylabel, zzlabel, fieldlabel, deltalabel, xxcalabel,yycalabel, phimlabel, phiMlabel;
    picture nalfalabel, ngammalabel, etalphalabel, alphalabel;
    picture gammalabel, thetalabel,vecclabel, thetalabelb, corelabel;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ShiftV := origin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    f 			:= (3,4,2);
    Spread 		:= 100;
    xxsize 		:= 0.8;
    xxpos 		:= 6;
    yysize 		:= 0.8;
    yypos 		:= 6;
    zzsize 		:= 1.3;
    zzpos 		:= 4;
    field 		:= (0.56,0.3,0.5);
    fieldpos 		:= 2;
    deltarad 		:= 0.45;
    deltapos 		:= 3;
    nalfad 		:= 0.8;
    theta 		:= 25;
    thetb 		:= 18;
    thetarad 		:= 0.3;
    thetapos 		:= 2;
    nalfapos 		:= 7;
    gammang 		:= -20;
    ngammapos 		:= 4;
    etalpharad 		:= 0.4;
    etalphapos 		:= 4;
    alphapos 		:= 6;
    gammapos 		:= 6;
    veccpos             := 3;

    xxcalabel   = btex $x_{\mbox{\protect \scriptsize layer}}$ etex;
    yycalabel   = btex $y_{\mbox{\protect \scriptsize layer}}$ etex;
    xxlabel 	= btex $x_{\mbox{\protect \scriptsize} M}$ etex;
    yylabel 	= btex $y_{\mbox{\protect \scriptsize M}}$ etex;
    zzlabel 	= btex $z_{\mbox{\protect \scriptsize M}}$ etex;
    corelabel 	= btex $\mbox{core}$ etex;
    deltalabel 	= btex $\beta$ etex;
    fieldlabel 	= btex $\mathbf{n}$ etex;
    thetalabel 	= btex $\theta_1$ etex;
    thetalabelb	= btex $\theta_2$ etex;
    nalfalabel 	= btex $\mathbf{n}_\alpha$ etex;
    ngammalabel = btex $\mathbf{ch}$ etex;
    etalphalabel= btex $\eta_\alpha$ etex;
    alphalabel 	= btex $\alpha$ etex;
    phimlabel 	= btex $-\phi_{\mbox{\protect \scriptsize 0}}/2$ etex;
    phiMlabel 	= btex $\phi_{\mbox{\protect \scriptsize 0}}/2$ etex;
    gammalabel 	= btex $\gamma$ etex;
    vecclabel 	= btex $\mathbf{c}$ etex;
    orig 		= (0,0,0);
    axxc 		= (xxsize,0,0);
    ayyc 		= (0,yysize,0);
    azzc 		= (0,0,zzsize);
    nalfa 		= (-nalfad*sind(theta)*cosd(90),nalfad*sind(theta)*sind(90),nalfad*cosd(theta));
    nalfaz 		= (-nalfad*sind(thetb)*cosd(90),nalfad*sind(thetb)*sind(90),-nalfad*cosd(thetb));
    alfacnt 		= (0,0,nalfad*cosd(theta));
    alfacntb 		= (0,0,-nalfad*cosd(thetb));
    radius 		= nalfad*sind(theta);
    radiusb 		= nalfad*sind(thetb);
    gradius             = 0.5*radius;
    alfa 		= rigorouscircle( alfacnt, (0,0,1), radius );
    alfab 		= rigorouscircle( alfacntb, (0,0,1), radiusb );
    ngorig              = orig+nalfa;
    ngorigb             = orig+nalfaz;
    gammacnt 		= cosd(theta)*cosd(theta)*nalfa+0.75*nalfa;
    gammacntb 		= cosd(theta)*cosd(theta)*nalfaz+0.75*nalfaz;
    radius 		:= cosd(theta)*sind(theta);
    gamma 		= rigorouscircle( gammacnt, nalfa, gradius );
    gammab 		= rigorouscircle( gammacntb, nalfaz, gradius );
    
    pickup pencircle scaled 1;

    for i=0 step 40 until 0:
       nazy 		:= gradius*sind(gammang+i)*ncrossprod( nalfa, axxc );
       nax 		:= gradius*cosd(gammang+i)*ncrossprod( ayyc, nalfa );
       ngamma[i] 	:= gammacnt+nax+nazy;
       draw rp(ngorig)..rp(ngamma[i]);
       nazyb 		:= gradius*sind(gammang+i)*ncrossprod( nalfaz, axxc );
       naxb 		:= gradius*cosd(gammang+i)*ncrossprod( ayyc, nalfaz );
       ngammab[i] 	:= gammacntb+naxb+nazyb;
       draw rp(ngorigb)..rp(ngammab[i]);
    endfor;
    
    pickup pencircle scaled 0.5;

    draw alfa dashed evenly;
    draw alfab dashed evenly;

    drawarrow rp(orig)..rp(axxc);
    drawarrow rp(orig)..rp(ayyc);
    drawarrow rp(-azzc)..rp(azzc);
    drawarrow rp(orig)..rp(2*field);
    pickup pencircle scaled 2;
    draw rp(orig)..rp(nalfa);
    draw rp(orig)..rp(nalfa);
    draw rp(orig)..rp(nalfaz);
    pickup pencircle scaled 0.5;
    draw rp(nalfa)..rp(gammacnt) dashed evenly;
    draw rp(nalfaz)..rp(gammacntb) dashed evenly;
%    draw subpath (3.05,1.05) of alfa;
%    drawarrow subpath (4,3) of alfa;
%    draw subpath (0,1.3) of gamma;
%    draw subpath (7.5,8) of gamma;
%    drawarrow subpath (1,2) of gamma;
    anglinen( field, azzc, orig, deltarad, deltalabel, deltapos );	
    anglinen( azzc, nalfa, orig, thetarad*1.5, thetalabel, thetapos );	
    anglinen( nalfaz, -azzc, orig, thetarad*2, thetalabelb, thetapos );	
%    anglinen( nalfa, nalfaz, orig, thetarad*1.25, btex $\Delta$ etex, 5 );	
%    anglinen( ngamma0, field, orig, etalpharad, etalphalabel, etalphapos );
    labeln( xxlabel, rp(axxc), xxpos );
    labeln( yylabel, rp(ayyc), yypos );
    labeln( zzlabel, rp(azzc), zzpos );
    labeln( corelabel, rp(orig+(0.2,-0.1,-0.2)), 1 );
    labeln( fieldlabel, rp(2*field), fieldpos );
%    labeln( nalfalabel, rp(nalfa), nalfapos );
    labeln( ngammalabel, rp(ngamma0), ngammapos );
    labeln( ngammalabel, rp(ngammab0), 1 );
    labeln( phiMlabel, point 1 of alfa, 7 );
    labeln( phimlabel, point 4 of alfa, 5 );
    labeln( phiMlabel, point 0.4 of alfab, 7 );
    labeln( phimlabel, point 3.5 of alfab,4 );
%    labeln( gammalabel, point 8 of gamma, gammapos );
%    labeln( yycalabel, rp(gammacnt+1.3*radius*(-1,0,0)), 0);
%    labeln( xxcalabel, rp(gammacnt+1.3*radius*(0,1,-sind(theta))), 0);
    pickup pencircle scaled 1.5;
%    draw rp(alfacnt);
%    draw rp(gammacnt);
    pickup pencircle scaled 0.5;

    path arco, recta_um, recta_dois, recta_xlayer, path_ngamma;
    path arcob, recta_umb, recta_doisb ;
    pickup pencircle scaled 0.5;
    recta_um   = rp(alfacnt)--rp(alfacnt+nalfad*1.4*radius*(1,0,0));
    recta_dois = rp(alfacnt)--rp(alfacnt+nalfad*1.4*radius*(0,1,0));
    draw recta_um;  
    draw recta_dois;
    draw rp(nalfa)..rp((0,nalfad*sind(theta)*sind(90),-nalfad*cosd(thetb))) dashed evenly;
    
    recta_umb   = rp(alfacntb)--rp(alfacntb+nalfad*1.4*radius*(1,0,0));
    recta_doisb = rp(alfacntb)--rp(alfacntb+nalfad*1.4*radius*(0,1,0));
    draw recta_umb;  
    draw recta_doisb;

    pickup pencircle scaled 1;
    arco = (alfa cutbefore recta_dois) cutafter recta_um;
    draw arco;
    arcob = (alfab cutbefore recta_doisb) cutafter recta_umb;
    draw arcob;


    pickup pencircle scaled 0.5;
 %   drawarrow rp(gammacnt)--rp(gammacnt+1.3*radius*(-1,0,0));

    path_ngamma = rp(gammacnt)--rp(ngamma0)--rp((0,0,0));
%    draw path_ngamma;
%%% c vector
    draw rp(gammacnt)--rp(ngamma0) dashed evenly;
    draw rp(gammacntb)--rp(ngammab0) dashed evenly;

    recta_xlayer = rp(gammacnt)--rp(gammacnt+1.3*radius*(0,1,-sind(theta)));
%    drawarrow recta_xlayer;

    pickup pencircle scaled 0.5;

    arco := (gamma cutbefore recta_xlayer) cutafter path_ngamma; 
    draw gamma dashed evenly;
    draw gammab dashed evenly;
    undraw arco;
    
    pickup pencircle scaled 0.5;
    draw gamma dashed evenly;
    
%    draw subpath (0,0.9) of gamma;
%    draw subpath (7.6,8) of gamma;
%    draw subpath (1.5,8) of gamma;
    draw subpath (3.2,1.2) of alfa;
%    labeln( vecclabel, rp(ngamma0), veccpos);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    picture schemepic;
    schemepic = currentpicture;
  numeric initang, stampnum, stepang, stampoffset;
  stampnum = 15;
  stampoffset = 0.2mm;
  initang = 180/stampnum;
  stepang = 2*initang;
  drawoptions();
  for i=0 upto (stampnum-1):
    draw schemepic 
      shifted (stampoffset*dir(initang+i*stepang))
      withcolor (blue+green);
  endfor;
  draw schemepic; 
endfig;

verbatimtex \end{document} etex

end;
