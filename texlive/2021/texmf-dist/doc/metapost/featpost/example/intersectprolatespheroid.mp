% intersectprolatespheroid.mp
% L. Nobre G. 
% 2013

input featpost3Dplus2D;

beginfig(1);
  numeric ray, i, Num, zcoord;
  color startpoi, lindir, outpoi, auxpoi;
  ray = 0.6;
  zcoord = 0.75;
  startpoi = (1.5,1.5,zcoord);
  lindir = (4,3,0);
  Num = 12;
  spheroid( black, blue, ray );
  drawoptions( dashed (withdots scaled 0.25) );
  draw rigorouscircle( black, blue, ray );
  draw rigorouscircle( (0,0,zcoord), blue, ray*(1 +-+ zcoord) );
  drawoptions();
  for i=0 upto Num:
    auxpoi := startpoi+red*i/Num;
    outpoi := intersectprolatespheroid( black, blue, ray, auxpoi, lindir );
    draw rp(auxpoi)--rp(outpoi) withpen pencircle scaled 1pt;
  endfor;
endfig;

end.