% pathernon.mp
% L. Nobre G.
% 2004

input featpost3Dplus2D;

beginfig(0);
        f := (4.8,3,3.65);
	viewcentr := black;
	Spread := 50;
	horiZon := 0;
	ShadowOn := false;
	LightSource := (1,-1,7);
	SphericalDistortion := true;
	OverRidePolyhedricColor := true;
	numeric side, coluheight, blocheight, roofheight, coluray, fact;
	numeric acroside, acronum, auxa, raya, rayb , poshemi, lena, lenb;
	color roofcenter, ba, bb, bc, bd, coluvec, colucen, auxb, auxc, auxd;
	side = 1;
	roofheight = 0.1;
	coluheight = 0.7;
	acroside = 4;
	acronum = 20;
	coluray = 0.06;
	blocheight = 0.1;
	poshemi = 0.25*side;
	raya = 0.08;
	rayb = 0.12;
	lena = 0.13;
	lenb = 0.26;
	fact = 1.3;
	auxa = 0.5*side-coluray;
	coluvec = blue*coluheight;
	colucen = 0.5*coluvec;
	auxc = (0,-poshemi,0);
	auxd = (0,-poshemi,raya+rayb);
	setthestage( acronum, acroside ); 
	V1 := (0.5*side,-0.5*side,coluheight+blocheight);
	V2 := (0.5*side,0.5*side,coluheight+blocheight);
	V3 := (-0.5*side,0.5*side,coluheight+blocheight);
	V4 := (-0.5*side,-0.5*side,coluheight+blocheight);
	V5 := (0,-0.5*side,coluheight+roofheight+blocheight);
	V6 := (0,0.5*side,coluheight+roofheight+blocheight);
	makeface1(1,4,5);
	makeface2(1,2,6,5);
	makeface3(2,3,6);
	makeface4(3,4,5,6);
	ba = (auxa,-auxa,0);
	bb = (auxa,auxa,0);
	bc = (-auxa,auxa,0);
	bd = (-auxa,-auxa,0);
	roofcenter = (V1+V2+V3+V4+V5+V6)/6;
	auxb = roofcenter-blue*blocheight;
	V7 := (-0.5*side,-0.5*side,coluheight);
	getready( "draw_invisible( false, false, white, black )", roofcenter );
	getready( "rigorousdisc( 0, true, ba, coluray, coluvec )", ba+colucen);
	getready( "rigorousdisc( 0, true, bb, coluray, coluvec )", bb+colucen);
	getready( "rigorousdisc( 0, true, bc, coluray, coluvec )", bc+colucen);
	getready( "rigorousdisc( 0, true, bd, coluray, coluvec )", bd+colucen);
	getready("kindofcube(false,true,V7,0,0,0,side,side,blocheight)",auxb);
	getready( "whatisthis(auxc,lenb*red,lena*green,raya,fact)", auxc );
	getready( "spatialhalfsfear( auxd, blue, rayb )", auxd );
	doitnow;
	produce_vga_border;
endfig;

end;
