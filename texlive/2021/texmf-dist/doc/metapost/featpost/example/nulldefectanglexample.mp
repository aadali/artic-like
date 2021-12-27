% nulldefectanglexample.mp
% L. Nobre G.
% 2012

input featpost3Dplus2D;

prologues := 1;

beginfig(1);
  color po[];
  numeric ta, tb;
  pen stroke;
  stroke = pencircle scaled 3pt;
  ta = sind(36)/cosd(36);
  tb = sind(18)/cosd(18);
  po1 = (1,ta,0);
  po2 = (-1,ta,0);
  po3 = (-1,tb,0);
  po4 = (-1,-ta,0);
  po5 = (0,-ta,0);
  po6 = (1,-ta,0);
  po7 = (1,0,0);
  path frame;
  frame = rp(po1)--rp(po2)--rp(po4)--rp(po6)--cycle;
  unfill frame;
  drawoptions( withpen stroke );
  draw frame;
  drawoptions( dashed evenly );
  draw rp(black)--rp(po1);
  draw rp(black)--rp(po3);
  draw rp(black)--rp(po5);
  draw rp(black)--rp(po7);
endfig;

beginfig(2);
  numeric parang, refdista, refdistb, auxa, auxb;
  parang = 55;
  refdista = conorm(po1-po3);
  refdistb = 1 ++ ta;
  auxa = (ta-tb)*sind(18);
  auxb = (ta-tb)*cosd(18);
  po7 := (cosd(parang),0,-sind(parang));
  po6 := po7+po5;
  po1 := improvertex( black, refdistb, po7, ta, po3, refdista, po1 );
  po8 = po3+N(po3)*auxa;
  po2 := po8+auxb*ncrossprod(po3,ccrossprod(po1,po3));
  drawoptions( withpen stroke );
  path a, b, c, d;
  a = rp(black)--rp(po3)--rp(po4)--rp(po5)--cycle;
  b = rp(black)--rp(po5)--rp(po6)--rp(po7)--cycle;
  c = rp(black)--rp(po7)--rp(po1)--cycle;
  d = rp(black)--rp(po1)--rp(po2)--rp(po3)--cycle;
  unfill a; draw a;
  unfill b; draw b;
  unfill c; draw c;
  unfill d; draw d;
endfig;

beginfig(3);
  numeric a, b, c, i, j;
  a = 0.5;
  b = 0.5*sqrt(3);
  V1 := (b,-a,a);
  V2 := (2*b,-a,0);
  V3 := (b,a,a);
  V4 := (2*b,a,0);
  V5 := (b,a,0);
  V6 := (0,1,0);
  c = 6;
  for j=1 upto 2: 
    for i=1 upto 6:
      c := incr(c);
      V[c] := rotvecaroundanother( j*120, V[i], blue );
    endfor;
  endfor;
  makeface1(1,3,7,9,13,15);
  makeface2(1,2,4,3);
  makeface3(3,4,5);
  makeface4(3,5,6,7);
  makeface5(7,6,8);
  makeface6(7,8,10,9);
  makeface7(9,10,11);
  makeface8(9,11,12,13);
  makeface9(13,12,14);
  makeface10(13,14,16,15);
  makeface11(15,16,17);
  makeface12(15,17,18,1);
  makeface13(1,18,2);
  f := (5,0.9,9);
  draw_invisible( false, false, white, black );
endfig;

beginfig(4);
  numeric bsang, dista, distb, distc, distd, verdi, hordi, si, co;
  color cauxa, cauxb, cauxd;
  bsang = 15;
  cauxa = (-1-sind(22),-1-cosd(22),1);
  cauxb = (-1-cosd(22),-1-sind(22),1);
  cauxd = (1+cosd(15),sind(15)-1,1);
  V1 := (-1,1,1);
  V2 := (-1,-1,1);
  distc = conorm(V2-cauxb);
  V3 := (1,-1,1);
  dista = conorm(cauxa-V3);
  V4 := (1,1,1);
  distd = conorm(V4-cauxd);
  V6 := (-1-cosd(bsang)*cosd(22),1+sind(22),1-sind(bsang)*cosd(22));
  V7 := (-1-cosd(bsang)*cosd(22),-1-sind(22),1-sind(bsang)*cosd(22));
  distb = conorm(cauxa-cauxb);
  V8 := improvertex( V3, dista, V7, distb, V2, distc, cauxb-(0,0,0.5) );
  V5 := (X(V8),-Y(V8),Z(V8));
  verdi = 1-Z(V8);
  hordi = Y(V8)+1;
  si = verdi/(verdi++hordi);
  co = hordi/(verdi++hordi);
  V9 := (1+0.5*sqrt(2),-1+0.5*sqrt(2)*co,1-0.5*sqrt(2)*si);
  V10 := improvertex( V9, 1, V3, 1, V4, distd, (1,-1,0) );
  V11 := (X(V10),-Y(V10),Z(V10));
  V12 := (X(V9),-Y(V9),Z(V9));
  NF := 0;
  makeface1(1,2,3,4);
  makeface2(1,5,6);
  makeface3(1,2,7,6);
  makeface4(2,7,8);
  makeface5(2,3,9,8);
  makeface6(3,9,10);
  makeface7(3,4,11,10);
  makeface8(4,11,12);
  makeface9(1,5,12,4);
  f := (9,2,5);
  draw_invisible( false, false, white, black );  
endfig;

end.