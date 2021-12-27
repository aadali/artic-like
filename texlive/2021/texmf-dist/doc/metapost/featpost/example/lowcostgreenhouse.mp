% lowcostgreenhouse.mp
% L. Nobre G.
% 2014

prologues := 3;

verbatimtex \documentclass{article}\begin{document} etex;

beginfig(1);
  numeric u;
  u = 7mm;
  numeric numcenhexagons, marg;
  numcenhexagons = 2;
  marg = 0.2u;
  numeric wid, len, ray, basewid;
  basewid = 11*sqrt(3)/3;
  wid = basewid*u;
  len = 2*u*(3*numcenhexagons+basewid);
  draw (-marg,-marg)--(marg+wid,-marg)--
    (marg+wid,marg+len)--(-marg,marg+len)--cycle;
  path botcircle, topcircle;
  z1 = (wid/2,0);
  botcircle = halfcircle scaled wid shifted z1;
  z2 = (0,len);
  topcircle = botcircle yscaled -1 shifted z2;
  draw botcircle;
  draw topcircle;
  z3 = (0,len-wid/2);
  z4 = (wid,len-wid/2);      
  draw z3--z4 withcolor red;
  z5 = (0,wid/2);
  z6 = (wid,wid/2);      
  draw z5--z6 withcolor red;
  z7 = (0,-2marg);
  z8 = (wid,-2marg);
  drawdblarrow z7--z8;
  label.bot(btex $w$ etex, 0.5[z7,z8]);
  z9 = (-2marg, wid);
  z10 = (-2marg, 0);
  drawdblarrow z9--z10;
  label.lft(btex $w=\frac{11}{3}\sqrt{3}$ etex, 0.5[z9,z10]);
  numeric j, k;
  k = 11;
  for j=wid*(1-0.5*sqrt(3)) step 6u until len-sqrt(3)*wid:
    z[k] = (0,j);
    z[k+1] = (wid,j+sqrt(3)*wid);
    draw z[k]--z[k+1];
    z[k+2] = (wid,j);
    z[k+3] = (0,j+sqrt(3)*wid);
    draw z[k+2]--z[k+3];
    z[k+4] = (0,j+0.5*sqrt(3)*wid-1.5u);
    z[k+5] = (wid,j+0.5*sqrt(3)*wid-1.5u);
    draw z[k+4]--z[k+5] withcolor (red+blue);
    z[k+6] = (0,j+0.5*sqrt(3)*wid+1.5u);
    z[k+7] = (wid,j+0.5*sqrt(3)*wid+1.5u);
    draw z[k+6]--z[k+7] withcolor (red+blue);
    k := k+8;
  endfor;
  z100 = z16+2*marg*right;
  z101 = z18+2*marg*right;
  drawdblarrow z101--z100;
  label.rt(btex $3$ etex, 0.5[z100,z101]);
  z200 = z19-2*marg*right;
  z201 = z17-2*marg*right;
  drawdblarrow z201--z200;
  label.lft(btex $1$ etex, 0.5[z200,z201]);
endfig;

verbatimtex \end{document} etex

end.