% nemsaddl.mp and nemsplay.mp and nemtwist.mp and nembends.mp
% L. Nobre G.
% 2002

input featpost3Dplus2D;

 

beginfig(1);
	f := 2*(5,3,4);
	Spread := 70;

	numeric hstep, hmax, hsdev, basestep, basemax, basesdev;
	numeric i, j, k, angsdev, cylength, cyradius, basen, hn;
	numeric vcx, vcy, vcz, hcurr, xcurr, ycurr, hbase, aone;
	numeric atwo, zcurr, angstep;
	color lenvec, currpos;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% splay, sddle and twist
	angstep    := 10; 
	cylength   := 0.6;
	cyradius   := 0.1;
	basen	   := 3;
	hn	   := 1;
	basestep   := cyradius*4.2;
	hstep	   := cylength*1.4;
	basesdev   := cyradius*0.003;
	hsdev	   := hstep*0.005;
	hbase	   := -0.8;
	angsdev	   := 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bend
%	angstep    := 15; 
%	cylength   := 0.6;
%	cyradius   := 0.1;
%	basen	   := 3;
%	hn	   := 1;
%	basestep   := cylength*1.6;
%	hstep	   := cyradius*3.4;
%	basesdev   := basestep*0.01;
%	hsdev	   := hstep*0.005;
%	hbase	   := -0.8;
%	angsdev	   := 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	basemax	   := basen*basestep;
	hmax	   := hn*hstep;
	for i=-basen upto basen:
	    for j=-basen upto basen:
		hcurr := hbase;
		forever:
		    exitunless hcurr < hmax+hbase;
		    hcurr   := hcurr + hstep;
		    zcurr   :=                hcurr + hsdev*normaldeviate;
		    xcurr   :=        i*basestep + basesdev*normaldeviate;
		    ycurr   :=        j*basestep + basesdev*normaldeviate;
		    aone    := angle((xcurr,ycurr))+angsdev*normaldeviate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%		    atwo    := ( i ++ j )*angstep + angsdev*normaldeviate;
%		    vcz     := cosd( atwo );
%		    vcy     :=  sind( atwo )*sind( aone );  %   pure-splay
%		    vcx     :=  sind( atwo )*cosd( aone );  %   pure-splay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%		    atwo    := ( i ++ j )*angstep + angsdev*normaldeviate;
%		    vcz     := cosd( atwo );
%		    vcy     :=  sind( atwo )*sind( aone );  % saddle-splay
%		    vcx     := -sind( atwo )*cosd( aone );  % saddle-splay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		    atwo    := ( i ++ j )*angstep + angsdev*normaldeviate;
		    vcz     := cosd( atwo );
		    vcy     := sind( atwo )*sind( aone+90 );%   pure-twist
		    vcx     := sind( atwo )*cosd( aone+90 );%   pure-twist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%		    atwo    := j *angstep + angsdev*normaldeviate;
%		    vcz     := sind( atwo );
%		    vcy     := cosd( atwo );                %   pure-bend
%		    vcx     := 0;                           %   pure-bend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		    currpos := ( xcurr, ycurr, zcurr );
		    lenvec  := cylength*( vcx, vcy, vcz ); 
		    rigorousdisc( 0, true, currpos-0.5*lenvec, cyradius, lenvec);
		endfor;
	    endfor;
	endfor;
endfig;

end;
