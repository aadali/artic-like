% cellevolve.mp
% L. Nobre G.
% IYP (2005)

def projectyz( expr theta, phi ) = 
	(cosd(theta)*sind(phi),sind(theta))
enddef;

def rotatangl( expr ori, rang ) =
	( (xpart ori)*cosd(rang)-(ypart ori)*sind(rang),
	  (xpart ori)*sind(rang)+(ypart ori)*cosd(rang))
enddef;

numeric vertn, horin, grids, thex, they, varx, xcomp;
numeric ymax, xmax, i, j, size, u, ycomp, theta, phi;
numeric power, xlim, ylim, frac, vary, allang, shift;
numeric perang, ampper, npers, velang, displamp, hdk;
pair actpos, direct, one, two, displvec;
path cirma, VGAborder;
color thedark;

u := 0.5mm;
allang := 15;
ampper := 20;
npers := 4;
velang := 5;
hdk := 0.6;
displamp := 8u;
frac := 1.5;
power := 0.8;
vertn := 20;
horin := 41;
size := 3u; 
grids := 8u;
ymax := floor(0.5*(vertn*frac-1));
xmax := floor(0.5*(horin*frac-1));
xlim := grids*0.5*horin;
ylim := grids*0.5*vertn;
varx := 180/(0.5*(horin-1));
vary := 90/(0.5*(vertn-1));
cirma := fullcircle scaled (grids);

	    VGAborder := (-xlim,-ylim)--         %
			 ( xlim,-ylim)--         %
			 ( xlim, ylim)--         %
			 (-xlim, ylim)--cycle;   %

beginfig(0);				
	pickup pencircle scaled u;
	for i=-xmax upto xmax:	
	    for j=-ymax upto ymax:
		actpos:=grids*rotatangl((i,j),allang);
		thex:=i*varx;
		they:=j*vary;
		theta:=90*cosd(they)*cosd(thex);
		phi:=60*sind(2*they)*sind(thex);
		perang := allang+ampper*sind(npers*they);
		direct:=rotatangl(projectyz(theta,phi),perang);
		displvec:=dir(velang)*displamp*cosd(npers*they);
		actpos := actpos+displvec; 
		one:=actpos-size*direct;
		two:=actpos+size*direct;
		thedark := (1-hdk*(ypart direct))*white;
		fill cirma shifted actpos withcolor thedark;
		draw one--two;
	    endfor; 
	endfor;
	clip currentpicture to VGAborder;
	undraw VGAborder;
endfig;

end;
