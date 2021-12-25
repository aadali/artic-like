% ultraeye.mp
% L. Nobre G.
% 2007

prologues := 1; 
beginfig(1);
  LaTeX_unitlength := 14.22636pt; % from the emp LaTeX package
  w := 15*LaTeX_unitlength;
  h := 10*LaTeX_unitlength;
  for i=h step -h/20 until h/20:
    fill fullcircle scaled i xscaled (w/h) withcolor blue;
    fill fullcircle scaled i withcolor red;
  endfor;
endfig;
end.
endinput;
