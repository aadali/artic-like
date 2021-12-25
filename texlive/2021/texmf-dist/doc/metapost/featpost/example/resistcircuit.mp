% resistcircuit.mp
% L. Nobre G. 
% IYP (2005)

input featpost3Dplus2D;

prologues := 1;

beginfig(1); 						
  boolean kind[];
  numeric largmolas, lengmolas, numolas, numturns, u, i, fraction, stepx;
  numeric totalx, totaly, actuax, actuay, numparal, numserie, stepy;
  path refpath;
  pair refpoint, startpoint, ini, fin;
  pickup pencircle scaled 2pt;
  fraction = 0.6;
  numolas = 9;
  u = 2mm;
  largmolas = 2u;
  lengmolas = 12u;
  numturns = 10;
  actuax = 1;
  actuay = 1;
  numparal = 0;
  numserie = 1;
  stepy = largmolas/fraction;
  stepx = 0.5*(1-fraction)*lengmolas;
  draw (left*stepx)--origin;
  for i=2 upto numolas:
    if uniformdeviate(2) > 1:
      kind[i] = false;
      numserie := incr( numserie );
    else:
      kind[i] = true;
      numparal := incr( numparal );
    fi;
  endfor;
  totalx = numserie*lengmolas;
  totaly = numparal*stepy;
  draw (right*totalx)--(right*(totalx+stepx));
  startpoint = (up*0.5*totaly);
  refpoint = startpoint;
  draw (down*0.5*totaly)--refpoint;
  refpath = springpath(origin,right*lengmolas,numturns,largmolas,fraction);
  draw refpath shifted refpoint;
  refpoint := refpoint + right*lengmolas;
  for i=2 upto numolas:
    if kind[i]:
      ini := startpoint+down*stepy*actuay;
      fin := ini + 0.5*(actuax-1)*lengmolas*right;
      draw ini--fin;
      draw refpath shifted fin;
      ini := fin + lengmolas*right;
      fin := ini + 0.5*(actuax-1)*lengmolas*right;
      draw ini--fin--refpoint;
      refpoint := refpoint + down*stepy*0.5;
      actuay := incr( actuay );
    else:
      draw refpath shifted refpoint;
      refpoint := refpoint + right*lengmolas;
      actuax := incr( actuax );
    fi;
  endfor;
endfig;

end.
