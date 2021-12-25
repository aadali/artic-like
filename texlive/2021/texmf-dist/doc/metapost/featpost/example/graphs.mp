% graphs.mp
% L. Nobre G.
% 2007

verbatimtex
\documentclass{article}
\usepackage{color,SIunits}
\begin{document}
etex

input graph;

beginfig(1);
  picture thegraph, thedash, thestar;
  path thecurve;
  pen thepen;
  numeric theu;
  theu = 1.7cm;
  z1 = (1,2);
  z2 = (3,8);
  z3 = (4,8);
  z4 = (2,8);
  z5 = (5,1);
  z6 = (9,1);
  z7 = (8,2);
  z8 = (7,3);
  z9 = (6,4);
  thestar = btex $\star$ etex scaled (3.5theu/in);
  thestar := thestar shifted -(center thestar); 
  thepen = pencircle scaled 2mm;
  thedash = dashpattern( off 3mm on 1mm off 4mm on 2mm );
  thecurve = z1..z2..z3..z4..z5..z6..z7..z8..z9;
  thegraph = begingraph(6theu,4theu);
    glabel.bot( btex $\unit{t}{/\;\second}$ etex, OUT );
    glabel.lft( btex $\unit{v}{/\;\metre\per\second}$ etex, OUT );
    gdraw thecurve withpen thepen withcolor blue dashed thedash;
    gdraw thecurve plot thestar;
    autogrid(grid.bot,grid.lft) withcolor red dashed evenly;
  endgraph;
  draw thegraph;  
endfig;

beginfig(2);
  draw begingraph(7cm,5cm);
    pickup pencircle scaled 1pt;
    path p;
    for x = -20 step 0.2 until -0.2:
      augment.p (x, sind(x*180/3.14159)/x);
    endfor
    augment.p (0, 1);
    for x = 0.2 step 0.2 until 20:
      augment.p (x, sind(x*180/3.14159)/x);
    endfor
    glabel.lrt (btex $\frac{\sin(x)}{x}$ etex, (-20,1));
    gdraw p;
  endgraph;
endfig;
  
verbatimtex
\end{document}
etex

end.

  
