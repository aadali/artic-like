% used many times in LaTeX documents (\usebox{})
% L. Nobre G. 
% 2001

beginfig(0); 						% grayrectanglebox
	numeric u, Nu;
	u = 0.5mm;
	Nu = 20;
	transform T;
	T = identity scaled (Nu*u) xscaled 1.33;
	pickup pensquare scaled u;
	draw (unitsquare transformed T) withcolor 0.7white;
endfig;

beginfig(1); 						% bigrectanglebox
	numeric u, Nu;
	u = 0.5mm;
	Nu = 70;
	transform T;
	T = identity scaled (Nu*u) xscaled 1.33;
	pickup pensquare scaled u;
	draw (unitsquare transformed T) withcolor 0.7white;
endfig;

end;
