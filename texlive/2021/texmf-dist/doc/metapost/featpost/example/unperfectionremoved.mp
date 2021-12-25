% unperfectionremoved.mp
% L. Nobre G.
% 2003

input featpost3Dplus2D;

Spread := 70;
f := (13,4,3);

beginfig(1);
	numeric i, len;
	numeric fws, NumLines, linex, inray, outay, axray;
	numeric wang;
	path oneline;
	color axe, aroc; 	
	len := 0.6;
	wang := 60;
	axe := (0,cosd(90-wang),sind(90-wang));
	fws := 4;
	outay := 0.45*fws;
	inray := 0.7*outay;
	NumLines := 30;
	for i=0 upto NumLines:
	    linex := fws*(i/NumLines-0.5);
	    oneline := rp((linex,0.5*fws,-0.5*fws))
                     --rp((linex,-0.5*fws,-0.5*fws))
                     --rp((linex,-0.5*fws,0.5*fws));
	    draw oneline withcolor (0.7*background);
	endfor;
	aroc := outay*(0,cosd(wang),sind(wang))-0.5*(0,fws,fws);
	rigorousdisc( inray, true, aroc, outay, axe*len );
endfig;

end.
