% hexagonalgrid.mp
% L. Nobre G.
% 2014

verbatimtex \documentclass{article}\begin{document} etex

prologues := 3;

input graph;

beginfig(1);
  numeric u, vertlim, horilim, dotsiz, inistep, mininistep, maxinistep;
  numeric securitymargin;
  u = 1.01cm;
  vertlim = 30u;
  horilim = vertlim/sqrt(2);
  dotsiz = 0.3u;
  mininistep = 4.867u;
  maxinistep = 4.868u;
  inistep = 0.5*(mininistep+maxinistep);
  securitymargin = 4u;
  pair upperleft, horizrigh;
  upperleft = (-15u,45u)+45u*dir(135);
  horizrigh = (35u,vertlim);
  pen defpen;
  defpen = pencircle scaled dotsiz;
  numeric i, j, Num;
  Num = 150;
  j = 0;
  for i=mininistep step (maxinistep-mininistep)/Num until maxinistep:
    pair Dh[]v[];  
    Dh0v0 = (0,vertlim); 
    Dh0v1 = (0,vertlim-i);
    numeric indh;
    indh = 0;
    forever:
      Dh[indh+1]v1 = whatever[upperleft,Dh[indh]v0];
      Dh[indh+1]v1 = whatever[horizrigh,Dh[indh]v1];
      xpart Dh[indh+1]v0 = xpart Dh[indh+1]v1;
      Dh[indh+1]v0 = whatever[horizrigh,Dh[indh]v0];
      exitif xpart Dh[indh+1]v0 > horilim;
      indh := indh + 1;
    endfor;
    z[j] = (i,xpart Dh[indh+1]v0 - horilim)/u;
    j := j + 1;
  endfor;
  path thecurve;
  picture theplus;
  theplus = btex $+$ etex;
  theplus := theplus shifted -(center theplus);
  thecurve = z0 for i=1 upto j-1: --z[i] endfor;
  picture thegraph;
  thegraph = begingraph(9u,6u);
    glabel.bot( btex inistep/u etex, OUT );
    glabel.lft( btex deviation/u etex, OUT );
    gdraw thecurve withpen defpen withcolor blue;
    gdraw thecurve plot theplus;
    autogrid(grid.bot,grid.lft) withcolor red dashed evenly;
  endgraph;
  draw thegraph;
endfig;

beginfig(2);  
  pair Dh[]v[];  
  Dh0v0 = (0,vertlim); 
  Dh0v1 = (0,vertlim-inistep);
  drawoptions( withpen defpen );
  draw Dh0v0;
  draw Dh0v1;
  numeric indh, indv, maxindv[];
  indh = 0;
  indv = 2;
  forever:
    Dh[indh+1]v[indv-1] = whatever[upperleft,Dh[indh]v[indv-2]];
    Dh[indh+1]v[indv-1] = whatever[horizrigh,Dh[indh]v[indv-1]];
    draw Dh[indh+1]v[indv-1];
    forever:
      xpart Dh[indh+1]v[indv-1] = xpart Dh[indh+1]v[indv];
      Dh[indh+1]v[indv] = whatever[upperleft,Dh[indh]v[indv-1]];
      draw Dh[indh+1]v[indv]; 
      xpart Dh[indh]v[indv] = xpart Dh[indh]v[indv-1];
      Dh[indh]v[indv] = whatever[horizrigh,Dh[indh+1]v[indv]];
      draw Dh[indh]v[indv];
      exitif ypart Dh[indh]v[indv] < -securitymargin;
      indv := indv + 1;
    endfor;
    maxindv[indh] = indv; 
    xpart Dh[indh+1]v0 = xpart Dh[indh+1]v1;
    Dh[indh+1]v0 = whatever[horizrigh,Dh[indh]v0];
    draw Dh[indh+1]v0;
    show xpart Dh[indh+1]v0 - xpart Dh[indh]v0;
    exitif xpart Dh[indh+1]v0 > horilim+securitymargin;
    indh := indh + 1;
    indv := 2;
    Dh[indh+1]v[indv-1] = whatever[upperleft,Dh[indh]v[indv-2]];
    Dh[indh+1]v[indv-1] = whatever[horizrigh,Dh[indh]v[indv-1]];
    draw Dh[indh+1]v[indv-1];
    xpart Dh[indh+1]v0 = xpart Dh[indh+1]v1;
    show xpart Dh[indh+1]v0 - xpart Dh[indh]v0;
    Dh[indh+1]v0 = whatever[horizrigh,Dh[indh]v0];
    draw Dh[indh+1]v0;
    indh := indh + 1;
  endfor;
  drawoptions();
  for i=2 step 2 until indh+1:
    for j=0 upto maxindv[i-2]-1:
      draw Dh[i]v[j+1]--Dh[i-1]v[j]--Dh[i]v[j]--cycle;
      draw Dh[i-1]v[j+1]--Dh[i-2]v[j]--Dh[i-1]v[j]--cycle;
    endfor;
  endfor;
  path borde, bordf;
  borde = (-dotsiz,-dotsiz)
    --(-dotsiz,vertlim+dotsiz)
    --(horilim+dotsiz,vertlim+dotsiz)
    --(horilim+dotsiz,-dotsiz)--cycle;
  bordf = origin--(0,vertlim)--(horilim,vertlim)--(horilim,0)--cycle;
  clip currentpicture to bordf;
  draw bordf;
  picture lasercutprofile;
  lasercutprofile = currentpicture;
  draw borde;
  %draw upperleft withpen defpen;
  %draw horizrigh withpen defpen;
endfig;

beginfig(3);
  draw lasercutprofile withpen pencircle scaled 0 withcolor red;
endfig;
  
verbatimtex \end{document} etex

end.