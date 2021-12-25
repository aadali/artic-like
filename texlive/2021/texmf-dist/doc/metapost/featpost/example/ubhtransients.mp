% ubhtransients.mp
% L. Nobre G.
% 2003

numeric SingleContinuosParameter;
SingleContinuosParameter = 0.95;

numeric vertn, horin, grids, thex, they, varx, xcomp;
numeric ymax, xmax, i, j, size, u, ycomp, theta, phi;
numeric power, xlim, ylim, frax, fray, allang, shift;
numeric perang, ampper, velang, displamp, hdk;
numeric velfunc, auxfunc, homang, vpower;
pair actpos, direct, one, two, displvec;
path cirma;
color thedark;
path VGAborder;

u := 0.5mm;
size := 1.14u; 
grids := 4u;
hdk := 0.47;

% This is for the sudden rotation
%allang := SingleContinuosParameter[0,80];
%homang := 0;

% This is for the UH evolution
%         allang := 80;
%	  homang := SingleContinuosParameter[0,80]; 
%	  ampper := 0;
%	  displamp := 0;
%	  frax := 1.2;
%	  fray := 1.6;
%	  velang := 0;
%	  power := 1;
%	  vpower := 1;

% The next homang parameter when in UBH evolution
% has two parts:
% Part1 - SingleContinuosParameter from 0 to 1 [0,45]
% Part2 - SingleContinuosParameter from 1 to 0 [85,45]
%	  homang := SingleContinuosParameter[0,45]; 
%	  allang := 85;
%	  ampper := SingleContinuosParameter[0,33];
%	  velang := 15;
%	  displamp := SingleContinuosParameter[0,5u];
%	  frax := 1.2;
%	  fray := 2.0;
%	  power := 1;
%	  vpower := 1;

% This is for the UB evolution
          allang := 90;
	  homang := 0; 
	  ampper := SingleContinuosParameter[0,84];
	  displamp := SingleContinuosParameter[0,9u];
	  frax := 0.8;
	  fray := 1.5;
	  velang := 0;
	  power := 1;
	  vpower := 1;

%power := SingleContinuosParameter[1,0.47];
%vpower := 1.0/power;

xlim := 115.00;
ylim := 86.025;
vertn := 2.0*ylim/grids; 
horin := 2.0*xlim/grids;
ymax := floor(0.5*(vertn*fray-1));
xmax := floor(0.5*(horin*frax-1));
%varx := 450/(0.5*(vertn-1));% in UH and UBH
varx := 350/(0.5*(vertn-1));% in UB

	    VGAborder := (-xlim,-ylim)--         %
			 ( xlim,-ylim)--         %
			 ( xlim, ylim)--         %
			 (-xlim, ylim)--cycle;   %
			 
    def mypower( expr base, expo ) =
        begingroup
	    save aux;
	    numeric aux;
	    if base=0:
	        aux = 1;
	    else:
		aux = abs(base)/base;
	    fi;
	    ( aux*(abs(base)**expo) )
	endgroup
    enddef;

    def produce_vga_border =
	begingroup
	    draw VGAborder withcolor background withpen pencircle scaled 0;
	    clip currentpicture to VGAborder
	endgroup
    enddef;
    
    beginfig(0);				
	pickup pencircle scaled u;
	cirma := fullcircle scaled (1.43*grids);
	for i=-xmax upto xmax:	
	    for j=-ymax upto ymax:
		actpos:=(grids*(i,j)) rotated (-allang);
		thex:=i*varx;
		they:=j*varx;
		auxfunc:=dir(velang) dotprod (thex,they);
		perang :=ampper*mypower(sind(auxfunc),power);
		theta:=180-allang+perang+homang;
		direct:=dir(theta);
		velfunc:=-displamp*mypower(cosd(auxfunc),vpower);
		displvec:=dir(velang)*velfunc;
		actpos := actpos+displvec; 
		thedark := (1-hdk*(ypart direct))*white;
		fill cirma shifted actpos withcolor thedark;
	    endfor; 
	endfor;
	for i=-xmax upto xmax:	
	    for j=-ymax upto ymax:
		actpos:=(grids*(i,j)) rotated (-allang);
		thex:=i*varx;
		they:=j*varx;
		auxfunc:=dir(velang) dotprod (thex,they);
		perang :=ampper*mypower(sind(auxfunc),power);
		theta:=180-allang+perang+homang;
		direct:=dir(theta);
		velfunc:=-displamp*mypower(cosd(auxfunc),vpower);
		displvec:=dir(velang)*velfunc;
		actpos := actpos+displvec; 
		one:=actpos-size*direct;
		two:=actpos+size*direct;
		draw one--two;
	    endfor; 
	endfor;
	produce_vga_border;
endfig;

end;
