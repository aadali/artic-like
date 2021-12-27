% lcurvature.mp
% L. Nobre G.
% 2002

input featpost3Dplus2D;

verbatimtex 
\documentclass{article} 
\usepackage{beton} 
\usepackage{concmath} 
\usepackage{ccfonts} 
\begin{document} 
etex

beginfig(1);
	f := (5,3,4);
	Spread := 70;
	numeric gridstep, anglstep, cylength, cyradius, numcylis;
	numeric i, j, curcoord, vcx, vcy, vcz, curangle, curnorma;
	numeric major, othcoord, othangle;
	color base, lenvec;
	cylength := 0.1;
	cyradius := 0.3;
	gridstep := cyradius*2.2;
	anglstep := 10;
	numcylis := 5;
	major := 0.35*gridstep*numcylis;
	for i=1 upto numcylis:			% bend
		curnorma := i-0.5*(numcylis+1);
		curangle := curnorma*anglstep;
		vcx := cylength*cosd(curangle);
		vcy := 0;	
		vcz := cylength*sind(curangle);	
		base := major*(sind(curangle),-1,-cosd(curangle));
		lenvec := (vcx,vcy,vcz);
		rigorousdisc( 0, true,base,cyradius,lenvec);
	endfor;
	for i=1 upto numcylis:			% twist
		curnorma := i-0.5*(numcylis+1);
		curcoord := curnorma*gridstep;
		curangle := curnorma*anglstep;
		vcx := 0;
		vcy := cylength*sind(curangle);	
		vcz := cylength*cosd(curangle);	
		base := (curcoord,-vcy/2-major,-vcz/2+major);
		lenvec := (vcx,vcy,vcz);
		rigorousdisc( 0, true,base,cyradius,lenvec);
	endfor;
	for i=1 upto numcylis:			% splay
		curnorma := i-0.5*(numcylis+1);
		curcoord := curnorma*gridstep;
		curangle := curnorma*anglstep;
		vcx := cylength*sind(curangle);
		vcy := 0;	
		vcz := -cylength*cosd(curangle);	
		lenvec := (vcx,vcy,vcz);
		for j=1 upto numcylis:
			othcoord := (j-0.5*(numcylis+1))*gridstep;
			base := major*(0,1,-cosd(curangle));
			base := base +(curcoord,othcoord,0);
			rigorousdisc( 0, true,base,cyradius,lenvec);
		endfor;
	endfor;
	for i=1 upto numcylis:			% saddle
		curnorma := i-0.5*(numcylis+1);
		curcoord := curnorma*gridstep;
		curangle := curnorma*anglstep;
		vcx := cylength*sind(curangle);
		for j=1 upto numcylis:
			othcoord := (j-0.5*(numcylis+1))*gridstep;
			othangle := (j-0.5*(numcylis+1))*anglstep;
			vcy := -cylength*sind(othangle);	
			vcz := -cylength*cosd(curangle)*cosd(othangle);	
			base := major*(0,1,1-cosd(curangle)+cosd(othangle));
			base := base + (curcoord,othcoord,0);
			lenvec := (vcx,vcy,vcz);
			rigorousdisc( 0, true,base,cyradius,lenvec);
		endfor;
	endfor;
endfig;

beginfig(2);
	f := (5,3,4);
	Spread := 70;
	numeric gridstep, cylength, cyradius, numcylis;
	numeric i, j, curcoord, curangle, curnorma;
	numeric othcoord, othangle, zcoord, coef, othnorma;
	color base, lenvec;
	cylength := 0.1;
	cyradius := 0.3;
	gridstep := cyradius*2.2;
	numcylis := 5;
	coef := 0.32;
	for i=1 upto numcylis:			% splay
	    curnorma := i-0.5*(numcylis+1);
	    curcoord := curnorma*gridstep;
	    curangle := -curcoord*2*coef;
	    for j=1 upto numcylis:
		othnorma := j-0.5*(numcylis+1);
		othcoord := othnorma*gridstep;
		othangle := -othcoord*2*coef;
		zcoord := -coef*(othcoord**2+curcoord**2);
		base := (curcoord,othcoord,zcoord);
		lenvec := cylength*N((-curangle,-othangle,1));
		rigorousdisc( 0, true,base,cyradius,lenvec);
	    endfor;
	endfor;
endfig;

verbatimtex 
\end{document} 
etex

end;
