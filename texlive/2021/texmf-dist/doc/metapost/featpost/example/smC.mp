% smC.mp
% P. J. Sebastiao

input featpost3Dplus2D;

% bash: export TEX=latex
% tcsh: setenv TEX latex

verbatimtex \documentclass{article} \usepackage{beton} \usepackage{concmath} \usepackage{ccfonts} \begin{document} etex

beginfig(1);
    color orig, axxc, ayyc, azzc;
    color field, nalfa, ngamma[], alfacnt, gammacnt, nazy, nax;
    path alfa, gamma;
    numeric nalfad, theta, radius, gammang, xxsize, yysize, zzsize;
    numeric xxpos, yypos, zzpos, fieldpos, deltarad, deltapos;
    numeric thetarad, thetapos, nalfapos, ngammapos,veccpos;
    numeric etalpharad, etalphapos, alphapos, gammapos;
    picture xxlabel, yylabel, zzlabel, fieldlabel, deltalabel, xxcalabel,yycalabel;
    picture nalfalabel, ngammalabel, etalphalabel, alphalabel;
    picture gammalabel, thetalabel,vecclabel;
    
    f 			:= (3,4,2);
    Spread 		:= 100;
    xxsize 		:= 0.8;
    xxpos 		:= 6;
    yysize 		:= 0.8;
    yypos 		:= 6;
    zzsize 		:= 1.3;
    zzpos 		:= 4;
    field 		:= (0.7,0.3,0.3);
    fieldpos 		:= 2;
    deltarad 		:= 0.45;
    deltapos 		:= 3;
    nalfad 		:= 1.0;
    theta 		:= 25;
    thetarad 		:= 0.3;
    thetapos 		:= 2;
    nalfapos 		:= 7;
    gammang 		:= 160;
    ngammapos 		:= 7;
    etalpharad 		:= 0.4;
    etalphapos 		:= 4;
    alphapos 		:= 6;
    gammapos 		:= 6;
    veccpos             := 3;

    xxcalabel   = btex $x_{\mbox{\protect \scriptsize layer}}$ etex;
    yycalabel   = btex $y_{\mbox{\protect \scriptsize layer}}$ etex;
    xxlabel 	= btex $x_{\mbox{\protect \scriptsize cone}}$ etex;
    yylabel 	= btex $y_{\mbox{\protect \scriptsize cone}}$ etex;
    zzlabel 	= btex $z_{\mbox{\protect \scriptsize cone}}$ etex;
    deltalabel 	= btex $\Delta$ etex;
    fieldlabel 	= btex $\mathbf{B}$ etex;
    thetalabel 	= btex $\theta$ etex;
    nalfalabel 	= btex $\mathbf{n}_\alpha$ etex;
    ngammalabel = btex $\mathbf{n}_\gamma$ etex;
    etalphalabel= btex $\eta_\alpha$ etex;
    alphalabel 	= btex $\alpha$ etex;
    gammalabel 	= btex $\gamma$ etex;
    vecclabel 	= btex $\mathbf{c}$ etex;
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
    
    for i=0 step 40 until 0:
       nazy 		:= radius*sind(gammang+i)*ncrossprod( nalfa, axxc );
       nax 		:= radius*cosd(gammang+i)*ncrossprod( ayyc, nalfa );
       ngamma[i] 		:= gammacnt+nax+nazy;
       drawarrow rp(orig)..rp(ngamma[i]);
    endfor;

    draw alfa dashed evenly;

    drawarrow rp(orig)..rp(axxc);
    drawarrow rp(orig)..rp(ayyc);
    drawarrow rp(orig)..rp(azzc);
    drawarrow rp(orig)..rp(field);
    drawarrow rp(orig)..rp(nalfa);
%    draw subpath (3.05,1.05) of alfa;
%    drawarrow subpath (4,3) of alfa;
%    draw subpath (0,1.3) of gamma;
%    draw subpath (7.5,8) of gamma;
%    drawarrow subpath (1,2) of gamma;
    anglinen( field, azzc, orig, deltarad, deltalabel, deltapos );	
    anglinen( azzc, nalfa, orig, thetarad, thetalabel, thetapos );	
    anglinen( ngamma0, field, orig, etalpharad, etalphalabel, etalphapos );
    labeln( xxlabel, rp(axxc), xxpos );
    labeln( yylabel, rp(ayyc), yypos );
    labeln( zzlabel, rp(azzc), zzpos );
    labeln( fieldlabel, rp(field), fieldpos );
    labeln( nalfalabel, rp(nalfa), nalfapos );
    labeln( ngammalabel, rp(ngamma0), ngammapos );
    labeln( alphalabel, point 1.8 of alfa, alphapos );
    labeln( gammalabel, point 8 of gamma, gammapos );
    labeln( yycalabel, rp(gammacnt+1.3*radius*(-1,0,0)), 0);
    labeln( xxcalabel, rp(gammacnt+1.3*radius*(0,1,-sind(theta))), 0);
    pickup pencircle scaled 1.5;
    draw rp(alfacnt);
    draw rp(gammacnt);
    pickup pencircle scaled 0.5;

    path arco, recta_um, recta_dois, recta_xlayer, path_ngamma;
    pickup pencircle scaled 0.5;
    recta_um   = rp(alfacnt)--rp(alfacnt+1.11*radius*(1,0,0));
    recta_dois = rp(alfacnt)--rp(alfacnt+1.11*radius*(0,1,0));
    draw recta_um;  
    draw recta_dois;

    pickup pencircle scaled 1;
    arco = (alfa cutbefore recta_dois) cutafter recta_um;
    draw arco;


    pickup pencircle scaled 0.5;
    drawarrow rp(gammacnt)--rp(gammacnt+1.3*radius*(-1,0,0));

    path_ngamma = rp(gammacnt)--rp(ngamma0)--rp((0,0,0));
    draw path_ngamma;
%%% c vector
    drawarrow rp(gammacnt)--rp(ngamma0);

    recta_xlayer = rp(gammacnt)--rp(gammacnt+1.3*radius*(0,1,-sind(theta)));
    drawarrow recta_xlayer;

    pickup pencircle scaled 1;

    arco := (gamma cutbefore recta_xlayer) cutafter path_ngamma; 
    draw gamma;
    undraw arco;
    
    pickup pencircle scaled 0.5;
    draw gamma dashed evenly;
    

%    draw subpath (0,0.9) of gamma;
%    draw subpath (7.6,8) of gamma;
%    draw subpath (1.5,8) of gamma;
    draw subpath (3.2,1.2) of alfa;
    labeln( vecclabel, rp(ngamma0), veccpos);
endfig;

verbatimtex \end{document} etex

end;
