% halfcirclesfear 
% L. Nobre G. 
% 2001

input featpost3Dplus2D;

% input la3Dmacros;

verbatimtex 
\documentclass{article} 
\usepackage{beton} 
\usepackage{concmath} 
\usepackage{ccfonts} 
\begin{document} 
etex

beginfig(1);
	color centre, axpin, axpon, axpun;
	numeric size, radii;
	f := (4,1.5,3.5);
	size = 0.5;
	radii = 1;
	centre = (0,0,0);
	axpin = (0,0,1);
	axpon = (0,1,0);
	axpun = (1,0,0);
	cartaxes( 1.2,1.2,0.8 );
	pickup pencircle scaled 2;
	draw rigorouscircle( centre, axpin, radii ) withcolor red;
	draw goodcirclepath( centre, axpin, radii );
	draw spatialhalfcircle( centre, axpon, size, true );
	draw spatialhalfcircle( centre, axpon, size, false ) dashed evenly;
	draw spatialhalfcircle( centre, axpun, size, true );
	draw spatialhalfcircle( centre, axpun, size, false ) dashed evenly;
	draw rigorousfearpath( centre, size ); 
endfig;


verbatimtex 
\end{document} 
etex

end;
