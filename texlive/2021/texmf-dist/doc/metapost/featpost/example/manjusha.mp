
numeric u;
u = 1cm;
pen thin, thick;
thin = pencircle scaled 1pt;
thick= pencircle scaled 3pt;

verbatimtex
\documentclass{article}
\usepackage{palatino}
\begin{document}
etex

beginfig(1);
  pair a, b, c, d;
  a = origin;
  b = (1u,2u);
  c = (4u,2u);
  d = (6u,0);
  path lines, curve;
  lines = a--b--c--d;
  curve = a..controls b and c..d; %%%%%%%%%%% This is quite powerfull.
  draw lines withpen thin withcolor red+blue;
  draw curve withpen thin withcolor blue;
  drawoptions( withpen thick );
  draw a;
  draw b;
  draw c;
  draw d;
endfig;

beginfig(4);
  pair a, b, c, d, e;
  a = origin;
  b = (2.5u,0.5u);
  c = (4u,2u);
  d = (6u,0);
  e = (4u,-0.5u);
  drawoptions( withcolor blue );
  picture thecaption;
  thecaption = btex Degree 3 etex;
  label( thecaption, e );
  path lines, curve;
  lines = a--b--c--d;
  curve = a..controls b and c..d;
  draw lines withpen thin withcolor red+blue;
  draw curve withpen thin withcolor blue;
  drawoptions( withpen thick withcolor red+green );
  draw a;
  draw b;
  draw c;
  draw d;
endfig;

beginfig(6);
  numeric n, i;
  n = 100;
  pair p[];
  for i=1 upto n:
    p[i] = (uniformdeviate(4u),uniformdeviate(3u));
    draw p[i] withpen thick withcolor green;
  endfor;
  path everywhere, convexhull;
  everywhere = p1 for i=2 upto n: --p[i] endfor --cycle;
  convexhull = makepath makepen everywhere; %%%%%%%%%%%%%% Wow!
  draw convexhull withpen thin withcolor red+green;
  for i=1 upto length convexhull:
    draw point i of convexhull withpen thick withcolor blue;
  endfor;
endfig;
    
verbatimtex
\end{document}
etex

end.
