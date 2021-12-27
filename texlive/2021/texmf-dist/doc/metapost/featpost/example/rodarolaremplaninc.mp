% rodarolaremplaninc.mp
% L. Nobre G.
% 2008

verbatimtex
\documentclass{article}
\usepackage{palatino,concmath,ccfonts}
\begin{document}
etex

beginfig(1);
  numeric Anginc;
  Anginc = 10+uniformdeviate(30);
  numeric ray; ray = 7mm;
  pickup pencircle scaled 1pt;
  numeric dh; dh = ray/sind( Anginc );
  numeric roll; roll = ray*(270-Anginc)*1.74533/100.0;
  numeric i, N; 
  N = 10+floor(uniformdeviate(10));
  z0 = dir(-Anginc);
  z1 = (-dh,0);
  z2 = z1+z0*(dh*cosd(Anginc)+roll+ray);
  z3 = z1+z0*dh*cosd(Anginc);
  z4 = z0*roll+(ray,0);
  path circ, arcangi, arca, spoke;
  spoke = (0.07[origin,z3])--z3;
  draw spoke; 
  picture roda;
  circ = fullcircle scaled (2*ray);
  draw circ;
  roda = currentpicture;
  for i = 1 upto N-1:
    draw spoke rotated ((Anginc-270)*i/N) shifted (z0*i*roll/N) withcolor
    0.5white;
  endfor;
  draw z1--z2--(-dh, ypart z2)--cycle;
  draw roda rotated (Anginc-270) shifted (z0*roll);
  arcangi = halfcircle scaled roll cutbefore (origin--roll*dir(180-Anginc));
  arca = arcangi shifted z2;
  z5 = point (0.5*(length arca)) of arca;
  pickup pencircle scaled 0.7pt;
  draw arca;
  label.lft( btex $\alpha$ etex, z5 ); 
  label.llft( btex $P_1$ etex, z3 ); 
  label.rt( btex $P_2$ etex, z4 ); 
endfig;

verbatimtex
\end{document}
etex

end.
