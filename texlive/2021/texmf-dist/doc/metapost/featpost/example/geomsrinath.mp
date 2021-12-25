% geomsrinath.mp
% L. Nobre G. 
% 2013

input featpost3Dplus2D;

verbatimtex
  \documentclass{article}
  \usepackage{palatino,mathpazo}
  \begin{document}
etex

beginfig(1);
  f := 3*(1.7,1.2,0.8);
  Spread := 125;
  numeric thetan, phiang, faca, phiray, cray, mar;
  color origr, nothr, mfxxp, mfxyp, mfyyp, mfzzp, vmagf, mfxyext, ac, bc;
  thetan = 55;
  phiang = 65;
  faca = 1.15;
  phiray = 0.485;
  cray = 0.067;
  origr = (0,0,0);
  mar = 0.15;
  vmagf = faca*(cosd(phiang)*cosd(thetan),sind(phiang)*cosd(thetan),sind(thetan));
  nothr = (0,0,1);
  bc = (faca-2*mar)*N(vmagf);
  ac = mar*N(vmagf);
    mfxxp = (X(vmagf),0,0);
    mfxyp = (X(vmagf),Y(vmagf),0);
    mfyyp = (0,Y(vmagf),0);
    mfzzp = (0,0,Z(vmagf));
    cartaxes( 0.8, 0.8, 1 );
    rigorousdisc(0,true,ac,cray,bc);
    draw rp(origr)--rp(mfxyp)--rp(vmagf)--rp(mfzzp) dashed evenly;
    draw rp(mfxxp)--rp(mfxyp)--rp(mfyyp) dashed evenly;
    if conorm(mfxyp)<phiray:
      mfxyext = phiray*N(mfxyp);
      draw rp(mfxyext)--rp(mfxyp) dashed evenly;
    fi;
    angline(mfxxp,mfxyp,origr,phiray,btex $\phi$ etex,bot);
    angline(vmagf,nothr,origr,0.5,btex $\theta$ etex,top);
    pickup pencircle scaled 1.15pt;
    drawarrow rp(origr)..rp(vmagf);
    label.rt(btex $\vec{n}$ etex,rp(nothr));
    label.rt(btex $\hat{a}$ etex,rp(vmagf));
    pickup pencircle scaled 2.6pt;
    draw rp(ac);
    draw rp(ac+bc);
endfig;

verbatimtex 
\end{document} 
etex

end;

