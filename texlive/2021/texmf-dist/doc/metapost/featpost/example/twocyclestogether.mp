% twocyclestogether.mp
% L. Nobre G.
% 2003

input featpost3Dplus2D;

f := (4.5,-6,5);

beginfig(1);
    numeric len[], ray[];
    color axe[], pos[], fix[], lab[], Ve[];
    color reflabpos, reflabpoi;
    picture txt[], dashpatt;
    numeric inclang, aberang, sampang;
    numeric squside, sqlevel, ul, outermagpos;
    numeric innerray, polesep, fracone, fracsam;
    numeric labelscale;
    path tmpath, aupath, copath;

    pickup pencircle scaled 1.2pt;
    dashpatt = evenly scaled 3;
    ul = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Change one of these at a time. %%%%%%%%
    inclang = 129;
    aberang = 19;
    outermagpos = 22ul;
    polesep = 1.7ul;
    len4 = 30ul;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sampang = 30;
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
    len5 = 3ul;
    len7 = 3.5ul;
    len8 = 7ul;

    len1 = (1-fracone)*(outermagpos-ray5-polesep);
    pos1 = (-outermagpos,0,0);
    rigorousdisc( 0, true, pos1, ray1, len1*(1,0,0) );

    len2 = fracone*(outermagpos-ray5-polesep);
    pos2 = (len1-outermagpos,0,0);
    tmpath := goodcirclepath( pos2, (1,0,0), ray1 ); 
    aupath := goodcirclepath( (-ray5-polesep,0,0), (1,0,0), ray2 );
    copath := twocyclestogether( tmpath, aupath );
    unfill copath;
    draw copath;
    draw aupath;
    
    Ve5 := fracsam*ray5*( sind( sampang ),
	                 cosd( sampang )*cosd( inclang -90 ),
	                 cosd( sampang )*sind( inclang -90 ) );
    Ve6 := fracsam*ray5*( -sind( sampang ),
	                 cosd( sampang )*cosd( inclang -90 ),
	                 cosd( sampang )*sind( inclang -90 ) );
    Ve8 := fracsam*ray5*( sind( sampang ),
	                 -cosd( sampang )*cosd( inclang -90 ),
	                 -cosd( sampang )*sind( inclang -90 ) );
    Ve7 := fracsam*ray5*( -sind( sampang ),
	                 -cosd( sampang )*cosd( inclang -90 ),
	                 -cosd( sampang )*sind( inclang -90 ) );
    tmpath := rp(Ve5)--rp(Ve6)--rp(Ve7)--rp(Ve8)--cycle;
    unfill tmpath;
    draw tmpath;

    axe6 = ( 0, cosd( inclang+aberang ), sind( inclang+aberang ) );
    pos6 = ray1*axe6;
    len6 = len4;
    axe7 = ( 0, cosd( inclang-aberang ), sind( inclang-aberang ) );
    pos7 = ( ray1 + 0.5len6 )*axe7;
    draw rp(pos6)--rp((0,0,0))--rp(pos7) dashed dashpatt;

    pos8 = pos7+len7*axe7;
    tmpath := goodcirclepath( pos7, axe7, ray7 ); 
    aupath := goodcirclepath( pos8, axe7, ray8 );
    copath := twocyclestogether( tmpath, aupath );
    unfill copath;
    draw copath;
    draw aupath;
    
    axe8 = axe7;
    rigorousdisc( 0, true, pos8, ray8, len8*axe8 );

    len9 = len2;
    pos9 = (ray5+polesep,0,0);
    ray9 = ray2;
    ray10 = ray1;
    tmpath := goodcirclepath( pos9, (1,0,0), ray9 ); 
    pos10 = pos9+len9*(1,0,0);
    aupath := goodcirclepath( pos10, (1,0,0), ray10 );
    copath := twocyclestogether( tmpath, aupath );
    unfill copath;
    draw copath;
    draw aupath;
    
    len10 = len1;
    rigorousdisc( 0, true, pos10, ray10, len10*(1,0,0) );

    produce_auto_scale;
endfig;

end.
