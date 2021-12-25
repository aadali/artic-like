% simplecarparam.mp
% L. Nobre G. 
% 2011

input featpost3Dplus2D;

verbatimtex 
\documentclass{article} 
\begin{document} 
etex

%f := 0.25*(2,2.15,15.3);
f := (3.5,8,1.1);
Spread := 220;

beginfig(1);
  color lencol, froncol, rearcol;
  color a, b, c, d, ca, cb, cc, cd;
  numeric lx, ly, lz, fx, fy, fz, hei;
  lx = 0.8;
  ly = 0.35;
  lz = 0.18;
  fx = 0.1;
  fy = 0.2;
  fz = 0.132;
  hei = 0.035;
  lencol = (lx, ly, lz);
  show lencol;
  froncol = (fx, fy, fz);
  show froncol;
  rearcol = (0.6*fx, 0.3*fy, 0.76*fz );
  show rearcol;
  a = (0.5*lx,0.5*ly,-0.5*lz);
  b = (0.5*lx-fx,0.5*ly,-0.5*lz);
  c = (0.5*lx-fx,0.5*ly+fy,-0.5*lz);
  d = (0.5*lx-fx-fz,0.5*ly+fy,-0.5*lz);
  ca= 0.5*(lx,ly,-lz);
  cb= 0.5*(lx,-ly,-lz);
  cc= 0.5*(-lx,ly,-lz);
  cd= 0.5*(-lx,-ly,-lz);
  drawoptions( withpen pencircle scaled 1.7pt withcolor 0.65blue );
  simplecar( black, black, lencol, froncol, rearcol );
  drawoptions( withpen pencircle scaled 1.3pt );
  draw rp(cd)--rp(cb) dashed evenly;
  draw rp(cd)--rp(cc) dashed evenly;
  draw rp(cd)--rp(-ca) dashed evenly;
  draw rp(ca)--rp(-ca) dashed evenly;
  draw rp(cb)--rp(-cb) dashed evenly;
  draw rp(cc)--rp(-cc) dashed evenly;
  draw rp(cd)--rp(-cd) dashed evenly;
  drawdblarrow rp(a)--rp(b);
  drawdblarrow rp(b)--rp(c);
  drawdblarrow rp(c)--rp(d);
  label.top( btex $\vec{o}$ etex, rp(black));
  drawoptions( withcolor 0.87red );
  labelinspace( true, 0.3[a,b], -red, hei*blue, btex Xf etex );
  labelinspace( true, 0.3[b,c], green, hei*blue, btex Yf etex );
  labelinspace( true, 0.3[c,d], -red, hei*blue, btex Zf etex );
endfig;

verbatimtex 
\end{document} 
etex

end;





