% stageforthree.mp
% L. Nobre G.
% 2003

input featpost3Dplus2D;

Spread := 35;
f := (13,7,3.5);

beginfig(1);
	numeric i, len, wang, reflen, frac, coordg;
	numeric fws, NumLines, inray, outay;
	path conepath, cira, cirb, ella, ellb, tuba, tubb, tubc;
	color axe, aroc, cubevertex, conecenter, conevertex; 	
	color allellaxe, ellaaxe, ellbaxe, pca, pea, pcb, peb;
	frac := 0.5;
	len := 0.6;
	wang := 60;
	axe := (0,cosd(90-wang),sind(90-wang));
	fws := 4;
	reflen := 0.2*fws;
	outay := 0.45*fws;
	inray := 0.7*outay;
	coordg := frac*fws;
	NumLines := 30;
	HoriZon := -0.5*fws;
	setthestage( 0.5*NumLines, 2*fws );
	cubevertex = (0.12*fws,-0.5*fws,-0.5*fws);
	kindofcube(false,true,cubevertex,180,0,0,0.65*fws,0.2*fws,fws);
	aroc := outay*(0,cosd(wang),sind(wang))-0.5*(0,fws,fws);
	rigorousdisc( inray, true, aroc, outay, axe*len );
	allellaxe := reflen*( 0.707, 0.707, 0 );
	ellaaxe := reflen*( 0.707, -0.707, 1.0 );
	ellbaxe := reflen*( -0.707, 0.707, 1.0 );
	conecenter = ( coordg, coordg, -0.5*fws );
	pca := ( coordg, -coordg, -0.5*fws );
	pcb := ( -coordg, coordg, -0.5*fws );
	pea := ( coordg, -coordg, 0.9*fws );
	peb := ( -coordg, coordg, 0.9*fws );
	cira := goodcirclepath( pca, blue, reflen );
	cirb := goodcirclepath( pcb, blue, reflen );
	ella := ellipticpath( pea, allellaxe, ellaaxe );
	ellb := ellipticpath( peb, allellaxe, ellbaxe );
	tuba := twocyclestogether( cira, ella );
	tubb := twocyclestogether( cirb, ellb );
	tubc := twocyclestogether( ella, ellb );
	unfill tubb; draw tubb;
	unfill tubc; draw tubc;
	unfill tuba; draw tuba;
	conevertex = conecenter + ( -3.5*reflen, 0, 0.8*fws );
	verygoodcone(false,conecenter,blue,reflen,conevertex);
endfig;

Spread := 23;
f := (6,13,3);
ShadowOn := true;

beginfig(2);
	numeric reflen, frac, coordg;
	numeric fws, NumLines;
	path ella, ellb;
	color axe, cubevertex; 	
	color allellaxe, ellaaxe, pca, pcb;
	frac := 0.5;
	wang := 60;
	axe := (0,cosd(90-wang),sind(90-wang));
	fws := 4;
	reflen := 0.35*fws;
	coordg := frac*fws;
	NumLines := 52;
	HoriZon := -0.5*fws;
	setthestage( 0.5*NumLines, 3.6*fws );
	cubevertex = (0.3*fws,-0.5*fws,-0.5*fws);
	kindofcube(false,true,cubevertex,180,0,0,0.6*fws,0.2*fws,fws);
	allellaxe := reflen*( 0.707, 0.707, 0 );
	ellaaxe := reflen*( 0.5, -0.5, 0 );
	pca := ( coordg, coordg, -0.5*fws );
	pcb := ( -coordg, coordg, -0.5*fws );

	whatisthis( pca, 0.5*ellaaxe, allellaxe, reflen, 0.3 );
	whatisthis( pcb, ellaaxe*sqrt(2), allellaxe/sqrt(2), 2*reflen, 0.6 );
	
endfig;

Spread := 23;
f := (13,6,4.5);
ShadowOn := true;
LightSource := 10*(4,-3,6);
beginfig(3);
	numeric reflen, frac, coordg;
	numeric fws, NumLines;
	path ella, ellb;
	color axe, cubevertex, conecenter, conevertex; 	
	color allellaxe, ellaaxe, pca, pcb;
	frac := 0.5;
	wang := 60;
	axe := (0,cosd(90-wang),sind(90-wang));
	fws := 4;
	reflen := 0.35*fws;
	coordg := frac*fws;
	NumLines := 45;
	HoriZon := -0.5*fws;
	setthestage( 0.5*NumLines, 3.3*fws );
	cubevertex = (0.3*fws,-0.5*fws,0);
	tropicalglobe( 7, cubevertex, 0.5*fws, axe );
	allellaxe := reflen*( 0.707, 0.707, 0 );
	ellaaxe := reflen*( 0.5, -0.5, 0 );
	pcb := ( -coordg, coordg, 0 );
	rigorousdisc( 0, true, pcb, 0.5*fws, -ellaaxe );
	conecenter = ( coordg, coordg, -0.5*fws );
	conevertex = conecenter + ( 0, 0, 0.9*fws );
	verygoodcone(false,conecenter,blue,reflen,conevertex);
endfig;

end.
