% durgayantra.mp
% L. Nobre G.
% 2010

prologues := 1;

numeric u, numthings, bang, tang, hei, rb, rt, ro;
numthings = 8;
bang = 75;
tang = 75;
rb = u+hei;
rt = rb+hei;
ro = rt+2.5mm;
numeric adang, reflen, marlen, sqrlen;
adang = 19.896;
reflen = 0.4*ro;
marlen = ro*sqrt(2)*(sind(45+adang)-0.5*sqrt(2));
hei = marlen/6;
ro+marlen+hei = 18mm; % defines u and hei
color fundo, risco;
fundo = (15*16+15,14*16+11,12*16+13)/255;
risco = (11*16-15,3*16-14,15*16-15)/255;

vardef makedurgayantra =
  path outercircle, firstcircle, seconcircle;
  draw fullcircle scaled (2*u);
  draw fullcircle scaled (2*rb);
  draw fullcircle scaled (2*rt);
  path pthingb, pthingbr;
  pthingb = (rt,0){dir(90-bang)}..
                  {dir(180/numthings+90-tang)}(ro*dir(180/numthings));
  pthingbr = pthingb yscaled (-1);
  numeric i;
  for i=1 upto numthings:
    draw (reverse pthingbr) rotated ((i+0.5)*360/numthings);
    draw pthingb rotated ((i+0.5)*360/numthings);
  endfor;
  z200 = (ro+marlen,0);
  z201 = (ro+marlen,1.5*reflen);
  z202 = (ro+0.5*marlen,1.5*reflen);
  z203 = (ro+0.5*marlen,0.5*reflen);
  z204 = (ro,0.5*reflen);
  z205 = (ro,ro);
  z100 = z200-(hei,0);
  z101 = z201-(hei,hei);
  z102 = z202+(hei,-hei);
  z103 = z203+(hei,-hei);
  z104 = z204-(hei,hei);
  z105 = z205-(hei,hei);
  z300 = z200+(hei,0);
  z301 = z201+(hei,hei);
  z302 = z202-(hei,-hei);
  z303 = z203-(hei,-hei);
  z304 = z204+(hei,hei);
  z305 = z205+(hei,hei);
  path octp[], octpr[];
  numeric j;
  octp[1] = z200--z201--z202--z203--z204--z205;
  octp[0] = z100--z101--z102--z103--z104--z105;
  octp[2] = z300--z301--z302--z303--z304--z305;
  for j=0 upto 2:
    octpr[j] = (reverse octp[j]) yscaled (-1) ;
    for i=1 upto 4:
      draw octpr[j] rotated (i*90);
      draw octp[j] rotated (i*90);
    endfor;
  endfor;
  z1 = (0,-u);
  z2 = u*0.5*(sqrt(3),1);
  z3 = u*0.5*(-sqrt(3),1);
  draw z1--z2--z3--cycle;
  numeric refang, xr, dr;
  xr = sqrt(15);
  dr = (xr-sqrt(3))*u*0.25;
  refang = angle(xr,1); show refang;
  z11 = u*0.25*(xr,-1);     
  z12 = 0.5[z2,z3]+dr*dir(60);   
  z13 = 0.5[z2,z3]+dr*dir(120);   
  z14 = u*0.25*(-xr,-1);  
  z15 = 0.5[z1,z3]+dr*dir(-120);
  z16 = 0.5[z2,z1]+dr*dir(-60);   
%  z11 = u*dir(-refang);     
%  z12 = u*dir(90-refang);   
%  z13 = u*dir(90+refang);   
%  z14 = u*dir(refang-180);  
%  z15 = u*dir(3*refang-180);
%  z16 = u*dir(-3*refang);   
  draw z11--z12--z15--z16--z13--z14--cycle;
enddef;

beginfig(1);
  numeric hs, al, do, m, diamfur, mfur;
  hs = 20mm;
  al = 16cm;
%  do = 3mm;
  do = 0;
  m = 3.8mm;
%  diamfur = 3mm;
  mfur = 6mm;
  z21 = (hs,m-hs);
  z22 = (hs,al-hs+do-m);
  z23 = (hs-m,al-hs+do);
  z24 = (m-hs,al-hs+do);
  z25 = (-hs,al-hs+do-m);
  z26 = (-hs,m-hs);
  z27 = (m-hs,-hs);
  z28 = (hs-m,-hs);
  fill z21---z22...z23---z24...z25---z26...z27---z28...cycle withcolor fundo;
%  draw (hs,al-hs)--(-hs,al-hs) dashed evenly;
%  draw fullcircle scaled diamfur shifted (mfur-hs,al-hs-mfur);
%  draw fullcircle scaled diamfur shifted (hs-mfur,al-hs-mfur);
  
  drawoptions( withpen pencircle scaled 0.3mm withcolor risco );
  makedurgayantra;
  drawoptions( withpen pencircle scaled 0.1mm withcolor risco );
  numeric lettersiz, lr, ls, lh;
  lettersiz = 1.1mm;
  lr = 1+0.5*sqrt(2);
  ls = 0.9mm; 
  lh = lettersiz*lr;
  path lepa[];
  pair av, auxa, lc, auxd;
  numeric lph[], auxb, auxc;
  auxc = 0.5*lettersiz;
  lc = (hs-m,0.5*al-hs);
  av = (auxc,lh);
  lepa0 = origin--av--(lettersiz,0);
  lph0 = 0;
  auxa = (lh-lettersiz)*unitvector(av);
  lepa7 = auxa--(lettersiz-xpart auxa, ypart auxa);
  lph7 = 0;
  lepa1 = origin--(0,lh)---av...(auxc,lh-lettersiz)---(0,lh-lettersiz);
  lph1 := lettersiz+ls;
  lepa2 = origin--(0,lh)--(auxc,lh-auxc)--(lettersiz,lh)--(lettersiz,0);
  lph2 = 2*(lettersiz+ls);
  auxb = 0.25*lettersiz*sqrt(2);
  lepa3 = (auxc-auxb,0.5*lh+auxb)--(auxc+auxb,0.5*lh-auxb);
  lepa13 = (auxc-auxb,0.5*lh-auxb)--(auxc+auxb,0.5*lh+auxb);
  lph3 = 3*(lettersiz+ls);
  lph13 = lph3;
  lepa4 = (0,lh)--origin--(lettersiz,0);
  lph4 = 4*(lettersiz+ls);
  lepa5 = origin--(0,lh)--(lettersiz,0)--(lettersiz,lh);
  lph5 = 5*(lettersiz+ls);
  lepa6 = (lettersiz,lh-auxc){up}...av...(0,lh-auxc)---(0,auxc)...(auxc,0)...(lettersiz,auxc)--(auxc,auxc);
  lph6 = 6*(lettersiz+ls);
  lepa8 = (lettersiz,auxc)--(lettersiz,0);
  lph8 = lph6;
  auxd = sqrt(2)*auxc*(1,1);
  lepa9 = (0,lh-auxc){up}...av...(lettersiz,lh-auxc)...auxd---origin--(lettersiz,0);
  lph9 = 8*(lettersiz+ls);
  lepa10 = (0,auxc)---(0,lh-auxc)...av...(lettersiz,lh-auxc)---(lettersiz,auxc)...(auxc,0)...cycle;
  lph10 = lph9 + lettersiz + ls;
  lepa11 = (0,lh-ypart auxd)--(xpart auxd,lh)--(xpart auxd,0);
  lph11 = lph10 + lettersiz + ls;
  lepa12 = lepa10;
  lph12 = lph11 + ls + xpart auxd; 
  numeric i;
  for i=0 upto 13:
    draw lepa[i] rotated 90 shifted (lc+(0,lph[i]));
  endfor;
  picture finalpicture;
  finalpicture = currentpicture;
  currentpicture := nullpicture;
  drawoptions();
  draw finalpicture rotated (-45);
endfig;

beginfig(2);
%  drawoptions( withcolor risco );
  makedurgayantra;
endfig;

end.
