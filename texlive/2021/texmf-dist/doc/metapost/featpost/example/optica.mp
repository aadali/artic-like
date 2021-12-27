% optica.mp
% P. J. Sebastiao 

input featpost3Dplus2D;

% bash: export TEX=latex
% tcsh: setenv TEX latex

verbatimtex \documentclass{article} \usepackage{beton} \usepackage{concmath} \usepackage{ccfonts} \begin{document} etex

beginfig(1);
    color orig, axxc, ayyc, azzc;
    color field, nalfa, ngamma[], alfacnt, gammacnt, nazy, nax, proj, align;
    path alfa, gamma;
    numeric nalfad, theta, radius, gammang, xxsize, yysize, zzsize;
    numeric xxpos, yypos, zzpos, fieldpos, deltarad, deltapos,angrecto;
    numeric thetarad, thetapos, nalfapos, ngammapos,veccpos,betarad,betapos;
    numeric etalpharad, etalphapos, alphapos, gammapos;
    picture xxlabel, yylabel, zzlabel, fieldlabel, deltalabel, xxcalabel,yycalabel,aMlabel;
    picture nalfalabel, ngammalabel, etalphalabel, alphalabel;
    picture gammalabel, thetalabel,vecclabel,betalabel;
    
    f 			:= (3,4,2);
    nbeta               := 1.3;
    beta                := 12;
    betarad 		:= 0.6;
    betapos 		:= 1;
    deltarad 		:= 0.8;
    deltapos 		:= 2;
    nalfad 		:= 1.0;
    theta 		:= 15;
    delta 		:= 35;
    thetarad 		:= 0.3;
    thetapos 		:= 2;
    nalfapos 		:= 1;
    gammang 		:= 60;
    ngammapos 		:= 7;
    etalpharad 		:= 0.6;
    etalphapos 		:= 2;
    alphapos 		:= 2;
    alpharad 		:= 0.9;
    gammapos 		:= 4;
    Spread              := 100;
    xxsize 		:= 0.8;
    xxpos 		:= 6;
    yysize 		:= 1.3;
    yypos 		:= 6;
    zzsize 		:= 1.3;
    zzpos 		:= 4;
%    angrecto            := 0.2;

    xxcalabel   = btex $x_{\mbox{\protect \scriptsize }}$ etex;
    yycalabel   = btex $y_{\mbox{\protect \scriptsize }}$ etex;
    xxlabel 	= btex $\mbox{\protect \scriptsize ligth}$ etex;
    yylabel 	= btex $\mbox{\protect \scriptsize analyser}$ etex;
    zzlabel 	= btex $\mbox{\protect \scriptsize polariser}$ etex;
    deltalabel 	= btex $\Delta$ etex;
    aMlabel 	= btex $\bar{\alpha}$ etex;
    fieldlabel 	= btex $\mathbf{B}$ etex;
    thetalabel 	= btex $\Delta$ etex;
    betalabel 	= btex $\xi$ etex;
    nalfalabel 	= btex $\mathbf{n}_a$ etex;
    ngammalabel = btex $\mathbf{n}_\gamma$ etex;
    etalphalabel= btex $\theta$ etex;
    alphalabel 	= btex $\alpha$ etex;
    gammalabel 	= btex $\gamma$ etex;
    vecclabel 	= btex $\mathbf{c}$ etex;
    veccpos             = 3;
    field 		= (0.7,0,0);
    fieldpos 		= 2;
    orig 		= (0,0,0);
    axxc 		= (xxsize,0,0);
    ayyc 		= (0,yysize,0);
    azzc 		= (0,0,zzsize);
    nalfa 		= (0,nbeta*cosd(beta)*sind(delta),nbeta*cosd(beta)*cosd(delta));
    alfacnt 		= nalfa;
    radius 		= nbeta*sind(beta);
    alfa 		= rigorouscircle( alfacnt, (0,0,1), radius );
    gammacnt 		= nalfa;
    gamma 		= rigorouscircle( gammacnt, nalfa, radius );
    
    for i=0 step 40 until 0:
       nazy 		:= radius*sind(gammang+i)*ncrossprod( nalfa, axxc );
       nax 		:= radius*cosd(gammang+i)*ncrossprod( ayyc, nalfa );
       ngamma[i] 		:= gammacnt+nax+nazy;
       pickup pencircle scaled 1.5;
       drawarrow rp(orig)..rp(ngamma[i]);
    endfor;
       pickup pencircle scaled 0.5;

%    draw alfa dashed evenly;
    align = 1*(0,+sind(10),cosd(10)); 
    draw rp(orig)..rp(align);
    drawarrow rp(orig)..rp(axxc);
    drawarrow rp(orig)..rp(ayyc);
%    draw rp(orig)..rp(-axxc);
    drawarrow rp(orig)..rp(azzc);
%    drawarrow rp(orig)..rp(field);
    draw rp(orig)..rp(nalfa);
%    draw subpath (3.05,1.05) of alfa;
%    drawarrow subpath (4,3) of alfa;
%    draw subpath (0,1.3) of gamma;
%    draw subpath (7.5,8) of gamma;
%    drawarrow subpath (1,2) of gamma;
%    anglinen( field, azzc, orig, deltarad, deltalabel, deltapos );	
    anglinen( ngamma0, nalfa, orig, betarad, betalabel, betapos );	
    anglinen( azzc, align, orig, deltarad, deltalabel, deltapos );	
    anglinen( align, nalfa, orig, deltarad, aMlabel, deltapos );	
%    anglinen( ngamma0, field, orig, etalpharad, etalphalabel, etalphapos );
    labeln( xxlabel, rp(axxc), xxpos );
    labeln( yylabel, rp(ayyc), yypos );
    labeln( zzlabel, rp(azzc), zzpos );
%    labeln( fieldlabel, rp(field), fieldpos );
%    labeln( nalfalabel, rp(nalfa), nalfapos );
    labeln( ngammalabel, rp(ngamma0), ngammapos );
%    labeln( alphalabel, point 1.8 of alfa, alphapos );
    labeln( gammalabel, point 2.1 of gamma, gammapos );
%    labeln( xxcalabel, rp(gammacnt+1.3*radius*(1,0,0)), 1);
%    labeln( yycalabel, rp(gammacnt+1.3*radius*(0,1,-sind(theta))), 2);
    pickup pencircle scaled 0.5;
    draw rp((0,0,0))..rp(alfacnt);
    draw rp(gammacnt);
    pickup pencircle scaled 0.5;

    path arco, recta_um, recta_quatro, recta_dois, recta_ylayer, recta_xlayer, path_ngamma, recta_tres;
    color tmp, tmpo;
    pickup pencircle scaled 0.5;
    tmp = gammacnt+radius*(0,-cosd(delta),sind(delta));
    tmpo = gammacnt+radius*(0,cosd(delta),-sind(delta));
    recta_um   = rp(tmpo)--rp(gammacnt)--rp(tmp)--rp((0,0,Z(tmp)));
    draw recta_um;

%    align = (0,Z(tmp)*sind(10),Z(tmp));
    drawarrow rp(orig)--rp((0,Z(tmp)*sind(10),Z(tmp)));
    labeln( btex $\mbox{\protect \scriptsize rubbing}$ etex, rp((0,Z(tmp)*sind(10),Z(tmp))), 1 );

    
    angrecto=nbeta*cosd(beta);
    recta_quatro = rp((0,angrecto*sind(delta),angrecto*cosd(delta)))--rp((0,angrecto*sind(delta),0));
%    recta_tres = rp((0,angrecto*sind(delta)+radius*sind(gammang)*cosd(delta),angrecto*cosd(delta)))--rp((0,angrecto*sind(delta)+radius*sind(gammang)*cosd(delta),0));
    recta_tres = rp((0,angrecto*sind(delta)+radius*sind(gammang)*cosd(delta),angrecto*cosd(delta)-radius*sind(gammang)*sind(delta)))--rp((0,angrecto*sind(delta)+radius*sind(gammang)*cosd(delta),0));
    
%    recta_dois = rp(angrecto*cosd(delta)*(0,0,1))--rp((0,angrecto*sind(delta)+radius*sind(gammang)*cosd(delta),angrecto*cosd(delta)));
%    recta_dois = rp(gammacnt)--rp((0,angrecto*sind(delta)+radius*sind(gammang)*cosd(delta),angrecto*cosd(delta)));
    recta_dois = rp(gammacnt)--rp((0,nbeta*cosd(beta)*sind(delta),0));
%    draw recta_dois;
    draw recta_tres;
%    draw recta_quatro;

    pickup pencircle scaled 1;
    arco = (alfa cutbefore recta_dois) cutafter recta_um;
%    draw arco;


    pickup pencircle scaled 0.5;
    drawarrow rp(gammacnt)--rp(gammacnt+1.3*radius*(1,0,0));

%    path_ngamma = rp(gammacnt)--rp(ngamma0);
    path_ngamma = rp(gammacnt)--rp(ngamma0)--rp((0,0,0));
    draw path_ngamma;
%%% c vector
%    drawarrow rp(gammacnt)--rp(ngamma0);

    recta_xlayer = rp(gammacnt)--rp(gammacnt+1.3*radius*(1,0,0));
%    recta_xlayer = rp(gammacnt)--rp(gammacnt+1.3*radius*(0,1,-sind(theta)));
    drawarrow recta_xlayer;


%    recta_ylayer = rp(gammacnt)--rp(gammacnt+1.3*radius*(0,1,-sind(theta)));
%    drawarrow recta_ylayer;

    proj = gammacnt+radius*sind(gammang)*(0,cosd(delta),-sind(delta));
    recta_ylayer = rp(gammacnt)--rp(proj);
    draw recta_ylayer;
    draw rp(ngamma0)..rp(proj);
    draw rp(orig)..rp(proj);
    draw rp(ngamma0)..rp(1.1*etalpharad*(1,0,0)+proj);
    draw rp(1.1*etalpharad*(1,0,0)+proj)..rp(1.1*etalpharad*(1,0,0));
%    draw rp(radius*cosd(gammang)*(1,0,0))..rp(ngamma0);
    anglinen( proj, nalfa, orig, alpharad, alphalabel, alphapos );	

    pickup pencircle scaled 1;

    arco := (gamma cutbefore recta_xlayer) cutafter path_ngamma; 
    draw gamma;
%    draw arco;
    undraw arco;
    arco := (gamma cutbefore recta_ylayer) cutafter path_ngamma; 
%    draw gamma;
%    draw arco;
    undraw arco;
    
    pickup pencircle scaled 0.5;
    draw gamma dashed evenly;

%    draw subpath (0,0.9) of gamma;
%    draw subpath (7.6,8) of gamma;
%    draw subpath (1.5,8) of gamma;
%    draw subpath (3.2,1.2) of alfa;
%    labeln( vecclabel, rp(ngamma0), veccpos);
    anglinen( ngamma0, field, orig, etalpharad, etalphalabel, etalphapos );
endfig;

verbatimtex \end{document} etex

end;
