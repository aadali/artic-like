verbatimtex 
\documentclass{article} 
\usepackage{palatino} 
\usepackage{concmath} 
\usepackage{ccfonts} 
\begin{document} 
etex

beginfig(1);
  numeric u, dist, desv, rayf, aux;
  u = 1cm;
  dist = 2u;
  desv = 1.5u;
  rayf = 1.2;
  ahangle := 30;
  aux = (rayf-1.0)*u;
  z1 = (dist,desv);
  z2 = (dist,-0.5u);
  z3 = (dist-aux,0.33*desv);
  z4 = (dist-aux,0.67*desv);
  drawdblarrow origin--(dist,0);
  drawarrow subpath (1,2) of fullcircle scaled (rayf*u);
  drawarrow z3--z4;
  pickup pencircle scaled 2pt;
  draw fullcircle scaled u;
  pickup pencircle scaled 5pt;
  draw z2--z1 withcolor 0.5white;
  label.urt(btex $i_2$ etex, (0,0.5*rayf*u) );
  label.llft(btex $i_1$ etex, z4 );
  label.top(btex $l$ etex, (0.5*dist,0) );
endfig;

verbatimtex 
\end{document} 
etex

end.