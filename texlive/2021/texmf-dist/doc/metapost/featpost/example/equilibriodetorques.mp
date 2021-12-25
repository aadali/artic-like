% equilibriodetorques.mp
% L. Nobre G.
% 2012

verbatimtex \documentclass{article}\begin{document} etex

beginfig(1);
  numeric u, longarm, smallarm, senodenom, senonumer, seno, ang, aux;
  u = 5mm;
  smallarm = 3u;
  longarm = 5u;
  senodenom = 10;
  senonumer = 8;
  seno = senonumer/senodenom;
  ang = angle(1+-+seno,seno);
  pen mainpen;
  mainpen = pencircle scaled (0.15u);
  path barpath, angpath, lenpath;
  barpath = (longarm,0.25u)--(-smallarm,0.25u)--
            (-smallarm,-0.25u)--(longarm,-0.25u)--cycle;
  filldraw barpath withpen mainpen withcolor (0.5white);
  drawarrow (-smallarm-u,0)--(longarm+u,0) dashed evenly;
  drawarrow (0,-2u)--(0,3u) dashed evenly;
  draw origin withpen mainpen;
  pair a, b, fa, fb;
  a = (-smallarm,0);
  b = (longarm,0);
  angpath = (halfcircle cutafter (origin--dir(ang))) scaled (2.5u) shifted a;
  draw angpath;
  fa = a+2.5u*dir(ang);
  fb = b+up*abs(fa)*seno*smallarm/longarm;
  drawarrow a--fa withpen mainpen;
  drawarrow b--fb withpen mainpen;
  aux = 0.45*(ypart fb);
  lenpath = (0,aux)--(longarm,aux);
  drawdblarrow lenpath;
  label.rt(btex $\theta$ etex, point (0.5*(length angpath)) of angpath );
  label.lft(btex $\vec{F}_1$ etex, fa);
  label.rt(btex $\vec{F}_2$ etex, fb);
  label.llft(btex $Y$ etex, (0,3u));
  label.llft(btex $X$ etex, (longarm+u,0));
  label.top(btex $r_2$ etex, point 0.5 of lenpath );
endfig;
  
verbatimtex \end{document} etex

end.