% unperfection.mp
% L. Nobre G.
% 2003

input featpost3Dplus2D;

Spread := 70;
f := (10,4,5);

beginfig(1);
	numeric i, len, polelen;
	numeric fws, NumLines, linex, inray, outay, axray;
	numeric wang, poleray;
	path aro[], oneline;
	color axe, pole, polec, aroc; 	
	picture floorwall;
	len := 0.6;
	wang := 60;
	axe := (0,cosd(90-wang),sind(90-wang));
	pole := (0,1,0);
	fws := 4;
	polelen := 0.8*fws;
	outay := 0.45*fws;
	inray := 0.7*outay;
	poleray := 0.07*inray;
	NumLines := 30;
	for i=0 upto NumLines:
	    linex := fws*(i/NumLines-0.5);
	    oneline := rp((linex,0.5*fws,-0.5*fws))
                     --rp((linex,-0.5*fws,-0.5*fws))
                     --rp((linex,-0.5*fws,0.5*fws));
	    draw oneline withcolor (0.7*background);
	endfor;
	aroc := outay*(0,cosd(wang),sind(wang))-0.5*(0,fws,fws);
	polec := (0,(outay-inray)*cosd(wang)-fws/2,
                    (outay+inray)*sind(wang)-poleray-fws/2);
	aro1 := goodcirclepath( aroc, axe, inray );
	aro2 := goodcirclepath( aroc+len*axe, axe, inray );
	aro3 := buildcycle( aro1, aro2 );
	floorwall := currentpicture;
	clip floorwall to aro3;
	aro4 := spatialhalfcircle( aroc+len*axe, axe, inray, false );
	aro5 := spatialhalfcircle( aroc+len*axe, axe, outay, false );
	aro6 := reverse aro5;
	aro7 := aro4--aro6--cycle;
	fill aro7 withcolor (0.6*background);
	rigorousdisc( 0, true, aroc, inray, axe*len );
	draw floorwall;
	draw aro3;
	rigorousdisc( 0, true, polec, poleray, pole*polelen );
	aro8 := spatialhalfcircle( aroc+len*axe, axe, inray, true );
	aro9 := spatialhalfcircle( aroc+len*axe, axe, outay, true );
	aro10 := reverse aro9;
	aro11 := aro8--aro10--cycle;
	fill aro11 withcolor (0.6*background);
	rigorousdisc( 0, false, aroc, outay, axe*len );
endfig;

end.
