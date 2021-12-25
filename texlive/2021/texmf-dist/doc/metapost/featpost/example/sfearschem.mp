% sfearchem.mp
% L. Nobre G.
% 2003

input featpost3Dplus2D;

verbatimtex 
\documentclass{article} 
\usepackage{beton} 
\usepackage{concmath} 
\usepackage{ccfonts} 
\begin{document} 
etex

defaultfont:="cmssu30";

beginfig(1);
	f := (4,2,1.5);
	cartaxes( 1.5,1.5,0.75 );
	drawoptions( withcolor red withpen pencircle scaled 1.3pt );
	draw goodcirclepath( (0,0,0), (0,0,1), 1 );
	drawoptions();
	draw rigorouscircle( (0,0,0), (0,0,1), 1 );
	draw spatialhalfcircle( (0,0,0), (0,0,1), 0.5, true );
	draw spatialhalfcircle( (0,0,0), (0,0,1), 0.5, false ) dashed evenly;
	draw rigorousfearpath( (0,0,0), 0.5 );
	draw rp((0,0,0.5)) withpen pencircle scaled 1.3pt;
endfig;

verbatimtex 
\end{document} 
etex

end;
