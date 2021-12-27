% intersectorus.mp
% L. Nobre G. 
% 2014

input featpost3Dplus2D;

f := (1,-5,2);

beginfig(1);
  numeric ray, i, Num, zcoord;
  color startpoi, lindir, outpoi, auxpoi;
  ray = 0.4;
  zcoord = 0.5*ray;
  startpoi = (-1-0.5*ray,-1,zcoord);
  lindir = (0,3,-2.5);
  Num = 50;
  smoothtorus( black, blue, 1, ray );
  for i=0 upto Num:
    auxpoi := startpoi+2*(1+0.5*ray)*red*i/Num;
    outpoi := intersectorus( black, blue, 1, ray, auxpoi, lindir );
    if pointinsidetorus( auxpoi, black, blue, 1, ray ):
      draw rp(auxpoi)--rp(outpoi) withpen pencircle scaled 2pt withcolor red;
    else:
      draw rp(auxpoi)--rp(outpoi) withpen pencircle scaled 2pt;
    fi;
  endfor;
  draw rp(startpoi)--rp(auxpoi);
endfig;

end.