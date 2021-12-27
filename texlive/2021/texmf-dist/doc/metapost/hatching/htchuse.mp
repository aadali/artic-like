%%% addto hatchfill
%%% addto beginfig
%%% enddef endfig
%%% length withcolor dashed off on
% This is a sample file demonstrating how one can use the macros defined
% in hatching.mp.

input hatching.mp;

vardef star(expr n) =
 for i_:=0 upto 2n-1:
  if odd i_: 1/2 fi (right rotated (180*(i_/n))) --
 endfor cycle
enddef;

beginfig(100);
 path p; p:=unitsquare xscaled 50mm yscaled 25mm;
 hatchfill p withcolor red withcolor (45,2mm,-.5bp) withcolor (-45,2mm,-.5bp);

endfig;

beginfig(101);
 interim hatch_match:=0;
 path p; p:=star(10) xscaled 35mm yscaled 20mm rotated 20;
 hatchfill p withcolor (0,1,.5);
 draw image(hatchfill p
  withcolor (45,3bp,-.5bp) withcolor (-45,3bp,-.5bp);
 ) withcolor red dashed evenly; % |evenly=dashpattern(on 3bp off 3bp)|
endfig;

beginfig(102);
 path p; p:=star(10) xscaled 35mm yscaled 20mm rotated 20;
 interim hatch_match:=0;
 hatchoptions(withcolor blue dashed evenly scaled 2);
 hatchfill p withcolor .75white withcolor (20,6bp,-.5bp);
 hatchoptions(withcolor (blue+green) dashed evenly shifted (3/2bp,0));
 hatchfill p withcolor (110,6bp,-.5bp);
endfig;

beginfig(103);
 path p; p:=star(10) xscaled 35mm yscaled 35mm;
 save draw_hatched_band;
 vardef draw_hatched_band(expr za,zb,a,l,d) =
  save n_; n_:=length(za-zb)/l;
  for i_:=0 upto ceiling n_:
   fill star(10) xscaled 2/5l yscaled 2/5l
     shifted (i_/n_)[za,zb] withcolor (i_/n_)[green,blue];
  endfor
 enddef;
 hatchfill p withcolor (red+green) withcolor (45,10bp,-1bp);
endfig;

end.
%%\end

