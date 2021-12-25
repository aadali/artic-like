% tdcircarrow.mp
% L. Nobre G.
% 2012

input featpost3Dplus2D;

f := (3,5,2);
Spread := 200;
TDAtiplen := 0.06; 
TDAhalftipbase := 0.04; 
TDAhalfthick := 0.025;   

beginfig(1);
  setthearena(11,1);
  numeric i, begsum, amp, margi;
  margi = 6;
  begsum = 0;
  amp = 36;
  for i=1 upto 4:
    tdcircarrow( black, blue, 0.6, begsum, amp-margi );
    begsum := begsum+amp;
    amp := 36*(i+1);
  endfor;
  tdarrow(black,0.6*ncrossprod( -f, blue ));
  tdarrow(black,0.6*ncrossprod( ccrossprod( -f, blue ), blue ));
endfig;

end.
  