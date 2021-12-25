% ropepatterns.mp
% L. Nobre G.
% 2011

input featpost3Dplus2D;

prologues := 2;

RopeColors := TableColors;                

beginfig(1);
  numeric i, u;
  u = 1cm;
  path auxpath, ropath;
  pair auxpair;
  for i=0 upto TableColors:
    RopeColorSeq[i] := i;
    auxpair := u*(i,TableColors);
    auxpath := fullcircle scaled u shifted auxpair;
    fill auxpath withcolor TableC[i];
    draw auxpath;
    label( decimal(i), auxpair );
  endfor;
  ropath = origin...(u*(0.5*TableColors,TableColors-2))...(u*TableColors*(1,0));
  ropepattern(ropath, u, 10);
  RopeColors := floor(0.5*TableColors);                
  ropepattern(ropath yscaled (-1) shifted (u*up*(TableColors-2)), u, 7);  
endfig;

end.