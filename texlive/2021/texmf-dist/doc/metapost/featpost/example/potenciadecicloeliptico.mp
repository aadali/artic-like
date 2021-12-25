% potenciadecicloeliptico.mp
% L. Nobre G.
% 2008

verbatimtex
\documentclass{article}
\usepackage{palatino,concmath,ccfonts}
\begin{document}
etex

beginfig(1);
  numeric u,w,h;
  u = 5mm;
  w = 11u;
  h = 6u;
  drawarrow origin--(0,h);  
  drawarrow origin--(w,0);
  label.lft( btex $p$ etex, (0,h) );  
  label.bot( btex $V$ etex, (w,0) );
  numeric wm, hm, tl;
  tl = 0.05*h;
  wm = 0.1+uniformdeviate(0.07); 
  hm = 0.1+uniformdeviate(0.07);
  path elipse;
  elipse = fullcircle xscaled (w-2*w*wm) yscaled (h-2*h*hm) 
                      shifted (w/2,h/2);
  draw elipse;
  z1 = (w*wm,0);
  z3 = (w*(1-wm),0);
  z2 = (0,h*(1-hm));
  z4 = (0,h*hm);
  draw z1--(z1+tl*up);
  draw z3--(z3+tl*up);
  draw z2--(z2+tl*right);
  draw z4--(z4+tl*right);
  label.lft( btex $p_b$ etex, z2 );  
  label.bot( btex $V_a$ etex, z1 );
  label.lft( btex $p_d$ etex, z4 );  
  label.bot( btex $V_c$ etex, z3 );
endfig;

verbatimtex
\end{document}
etex

end.
