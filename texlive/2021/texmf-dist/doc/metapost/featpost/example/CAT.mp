beginfig(1);
	path pat[];
	numeric i, a, sidel;
	sidel = 2cm;
	pickup pencircle scaled 3pt;
	pair pa, pb, pc, pd;
	for i=4 upto 10:
		a := 360/i;
		pa := origin+(i-4)*right*1.5*sidel;
		pb := pa+right*sidel;
		pd := pa+dir(a)*sidel;
		pc := pb+dir(a)*sidel;
		pat[i] = pa--pb--pc--pd--cycle;
	endfor;
	z0 = (0,2.5*sidel);
	z1 = z0+2*right*sidel;
	z2 = z0+2*dir(60)*sidel;
	pat0 = z0--z1--z2--cycle;
	z3 = z1+1.5*right*sidel;
	z4 = z3+right*sidel;
	z5 = z4+dir(60)*sidel;
	z6 = z5+dir(120)*sidel;
	z8 = z3+dir(120)*sidel;
	z7 = z8+dir(60)*sidel;
	pat1 = z3--z4--z5--z6--z7--z8--cycle;
	z9 = z4+1.5*right*sidel;
	z10= z9+right*sidel;
	z11= z10+dir(60)*sidel;
	z12= z11+dir(120)*sidel;
	z14= z9+dir(60)*sidel;
	z13= z14+dir(120)*sidel;
	pat2 = z9--z10--z11--z12--z13--z14--cycle;
	z15= z10+1.5*right*sidel;
	z16= z15+right*sidel;
	z17= z16+dir(60)*sidel;
	z18= z15+2*dir(60)*sidel;
	pat3 = z15--z16--z17--z18--cycle;
	background := 0.85white;
	for i = 0 upto 10:
		unfill pat[i];
		draw pat[i];
	endfor;
endfig;

end.