% minimal.mp
verbatimtex \documentclass{beamer}\begin{document} etex
prologues := 3;
outputformat := "eps";
outputtemplate := "%j-%c.mps";
beginfig(1);
  numeric u, h, v;
  picture title;
  path arm;
  u = 4mm;
  h = 9u;
  v = 2u;
  z1 = (0,-v);
  z2 = (h,-v);
  z3 = (h,v);
  z4 = (0,v);
  arm = z1---z2..z3---z4..cycle;
  title = thelabel( btex minimal etex scaled 2, origin );
  unfill arm;
  draw arm withpen pencircle scaled 3bp;
  draw title rotated 7 shifted (0.5h,0);
endfig;
verbatimtex \end{document} etex
end
