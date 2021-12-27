% xcmplaca2buraquads.mp
% L. Nobre G.
% 2008

verbatimtex
\documentclass{article}
\usepackage{palatino,concmath,ccfonts}
\begin{document}
etex

beginfig(1);
  numeric fy, tf, fl, tk, fu, fd;
  fy = 0.25+uniformdeviate( 0.35 );
  tf = 0.4*fy/sqrt(2.0);
  fl = 0.3*tf+uniformdeviate( 0.7*tf );
  tk = fl/sqrt(2.0);
  fu = tk+uniformdeviate( 0.24-tk );
  fd = 0.26+uniformdeviate( 0.24-tk );
  numeric u,h,w;
  u = 5mm;
  w = 10u;
  h = 11*u*fy;
  drawarrow origin--(0,h);  
  drawarrow origin--(w,0);
  label.lft( btex $y$ etex, (0,h) );  
  label.bot( btex $x$ etex, (w,0) );
  numeric marg, tl;
  tl = 0.05*h;
  marg = 0.1*w;
  numeric xzero, yzero, xum, xdois, yum, ydois, lado, ang;
  xzero = w-marg;
  yzero = h-marg;
  xum = fu * xzero;
  xdois = fd * xzero;
  yum = 0.7 * yzero;
  ydois = 0.3 * yzero;
  lado =  fl * xzero;
  ang = 10+uniformdeviate( 26 );
  path placa, quadrad, quadrum, quadois;
  placa = origin--(xzero,0)--(xzero,yzero)--(0,yzero)--cycle;
  fill placa withcolor 0.5white;
  draw placa;
  quadrad = (origin--(lado,0)--(lado,lado)--(0,lado)--cycle) 
            shifted (-0.5*lado*(1,1)) rotated ang;
  quadrum = quadrad shifted (xum,yum);
  quadois = quadrad shifted (xdois,ydois);
  unfill quadrum;
  draw quadrum;
  unfill quadois;
  draw quadois;
  z1 = (xum,0);
  z2 = (xdois,0);
  z3 = (xzero,0);
  z4 = (0,yzero);
  draw z1--(z1+tl*up);
  draw z2--(z2+tl*up);
  draw z3--(z3+tl*up);
  draw z4--(z4+tl*right);
  label.lft( btex $y_0$ etex, z4 );  
  label.bot( btex $x_1$ etex, z1 );
  label.bot( btex $x_2$ etex, z2 );  
  label.bot( btex $x_0$ etex, z3 );
endfig;

verbatimtex
\end{document}
etex

end.
