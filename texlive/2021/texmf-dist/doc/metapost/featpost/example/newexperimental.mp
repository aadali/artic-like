% newexperimental.mp
% L. Nobre G.
% 2004

input featpost3Dplus2D;

%verbatimtex 
%\documentclass{article} 
%\usepackage{beton} 
%\usepackage{concmath} 
%\usepackage{ccfonts}
%\begin{document} 
%etex

verbatimtex 
\documentclass{article} 
\usepackage{palatino} 
\begin{document} 
etex

f := (4.5,-6,5);
%f := (-4.5,-6,5);

beginfig(1);
    numeric len[], ray[];
    color axe[], pos[], fix[], lab[];
    color reflabpos, reflabpoi, lincolor;
    picture txt[], dashpatt;
    numeric inclang, aberang, sampang, stepang;
    numeric squside, sqlevel, ul, outermagpos;
    numeric innerray, polesep, fracone, fracsam;
    numeric labelscale;
    path tmpath, aupath, copath;

    pickup pencircle scaled 1.2pt;
    dashpatt = evenly scaled 3;
    inclang = 150;
    aberang = 17;
    sampang = 30;
    stepang = 20;
    ul = 0.1;
    labelscale = 5;
    outermagpos = 30ul;
    squside = outermagpos; 
    sqlevel = -20ul;
    polesep = 15ul;
    fracone = 0.33;
    fracsam = 0.86;
    ray1 = 10ul;
    ray2 = 6ul;
    ray3 = 8ul;
    ray4 = 1.5ul;
    ray5 = 4ul;
    ray6 = 3ul;
    ray7 = 3.5ul;
    ray8 = 2ul;    
    len3 = 7ul;
    len4 = 20ul;
    len5 = 3ul;
    len7 = 3.5ul;
    len8 = 7ul;
    lincolor = 0.5background;

    reflabpos = ( -squside, 0, -sqlevel );     
    reflabpoi = ( -squside, squside, -sqlevel );     
    V1 := (squside,squside,-sqlevel);
    V2 := (-squside,squside,-sqlevel);
    V3 := (-squside,-squside,-sqlevel);
    V4 := (squside,-squside,-sqlevel);
%    makeface1(1,2,3,4);
%    draw_all_test( false );

    len1 = (1-fracone)*(outermagpos-ray5-polesep);
    pos1 = (-outermagpos,0,0);
    len2 = fracone*(outermagpos-ray5-polesep);
    pos2 = (len1-outermagpos,0,0);
    lab2 = reflabpos;
    txt2 = thelabel.top( btex magnet etex scaled labelscale, rp( lab2 ) );
    fix2 = (-outermagpos+0.5len1,0,ray1);
    axe3 = ( 0, cosd( inclang ), sind( inclang ) );
    pos3 = -( len3 + len4 + len5 )*axe3;
    lab3 = (1-fracone)*reflabpoi+fracone*V1;
    txt3 = thelabel.rt( btex rotor etex scaled labelscale, rp( lab3 ) );
    axe4 = axe3;
    pos4 = -( len4 + len5 )*axe4;
    fix3 = pos4+fracsam*ray3*(0,cosd( inclang -90 ),sind( inclang -90 ));
    axe5 = axe3;
    pos5 = -( len5 )*axe5;
    lab5 = (1-fracone)*reflabpoi+fracone*V2;
    txt5 = thelabel.rt( btex sample etex scaled labelscale, rp(lab5) );
    fix5 = fracsam*ray5*(0,cosd( inclang -90 ),sind( inclang -90 ));

    V5 := fracsam*ray5*( sind( sampang ),
	                 cosd( sampang )*cosd( inclang -90 ),
	                 cosd( sampang )*sind( inclang -90 ) );
    V6 := fracsam*ray5*( -sind( sampang ),
	                 cosd( sampang )*cosd( inclang -90 ),
	                 cosd( sampang )*sind( inclang -90 ) );
    V8 := fracsam*ray5*( sind( sampang ),
	                 -cosd( sampang )*cosd( inclang -90 ),
	                 -cosd( sampang )*sind( inclang -90 ) );
    V7 := fracsam*ray5*( -sind( sampang ),
	                 -cosd( sampang )*cosd( inclang -90 ),
	                 -cosd( sampang )*sind( inclang -90 ) );
    axe6 = ( 0, cosd( inclang+aberang ), sind( inclang+aberang ) );
    pos6 = ray1*axe6;
    len6 = len4;
    axe7 = ( 0, cosd( inclang-aberang ), sind( inclang-aberang ) );
    pos7 = ( ray1 + 0.5len6 )*axe7;
    lab6 = fracone*reflabpos+(1-fracone)*V3;
    txt6=thelabel.lft(btex camera etex scaled labelscale,rp(lab6));
    fix6 = pos6+len6*axe6;
    pos8 = pos7+len7*axe7;
    axe8 = axe7;
    lab8 = (1-fracone)*reflabpos+fracone*V3;
    txt8 = thelabel.lft( btex lamp etex scaled labelscale, rp( lab8 ) );
    fix8 = pos8+len8*axe8;
    len9 = len2;
    pos9 = (ray5+polesep,0,0);
    ray9 = ray2;
    ray10 = ray1;
    pos10 = pos9+len9*(1,0,0);
    len10 = len1;

    getready( "rigorousdisc( 0, true, pos1, ray1, len1*(1,0,0) )", pos1 );
    getready( "whatisthis(pos2,ray1*green,ray1*blue,conorm((-ray5-polesep,0,0)-pos2),ray2/ray1)", pos2 );
    getready( "rigorousdisc( 0, true, pos3, ray3, len3*axe3 )", pos3 );
    getready( "rigorousdisc( 0, true, pos4, ray4, len4*axe4 )", pos4 );
    getready( "rigorousdisc( 0, true, pos5, ray5, len5*axe5 )", pos5 );
    getready( "draw rp(V5)--rp(V6)--rp(V7)--rp(V8)--cycle", black );
    getready( "draw rp(pos6)--rp((0,0,0))--rp(pos7) dashed dashpatt",0.5*pos6);
    getready( "rigorousdisc( 0, true, pos6, ray6, len6*axe6 )", pos6 );
    getready( "whatisthis(pos7,ray7*red,-ray7*ncrossprod(red,axe7),conorm(pos8-pos7),ray8/ray7)", pos7 );
    getready( "rigorousdisc( 0, true, pos8, ray8, len8*axe8 );", pos8 );
    getready( "whatisthis(pos9,ray9*green,ray9*blue,conorm(pos10-pos9),ray10/ray9)", pos9 );
    getready( "rigorousdisc( 0, true, pos10, ray10, len10*(1,0,0) )", pos10 );
    doitnow;
    draw llcorner txt2--lrcorner txt2--rp(fix2) withcolor lincolor;
    draw lrcorner txt3--llcorner txt3--rp(fix3) withcolor lincolor;
    draw lrcorner txt5--llcorner txt5--rp(fix5) withcolor lincolor;
    draw llcorner txt6--lrcorner txt6--rp(fix6) withcolor lincolor;
    draw llcorner txt8--lrcorner txt8--rp(fix8) withcolor lincolor;
    draw txt2;
    draw txt3;
    draw txt5;
    draw txt6;
    draw txt8;

    produce_auto_scale;
endfig;

verbatimtex 
\end{document} 
etex

end.
