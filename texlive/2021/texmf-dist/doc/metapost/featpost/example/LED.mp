% LED.mp
% L. Nobre G.
% 2009

prologues := 1;

beginfig(1);
  numeric u, hthick, vthick, lspace, espace, dspace, h, lbase, dbase;
  numeric dinray, halfh, tp, scal, marg, dcorrection;
  u = 1.5mm;
  hthick = 20u;
  vthick = 10u;
%  lspace = 15u;
  lspace = 8u;
  espace = 8u;
  dspace = 8u;
  h = 3*vthick+2*espace;
%  lbase = 38u;
  lbase = 45u;
%  dbase = 16u;
  dbase = 19u;
  halfh = h/2;
  dinray = halfh-vthick;
  tp = 1.15;
  dcorrection = 1.3u;
  scal = (h-vthick-espace)/h;
  marg = 4u;
  z200 = (marg,vthick+espace-marg);
  path pathl, pathe, pathd, pathu;
  z0 = (0,0);
  z1 = (lbase,0);
  z2 = (lbase,vthick);
  z3 = (hthick,vthick);
  z4 = (hthick,h);
  z5 = (0,h);
  pathl = z0--z1--z2--z3--z4--z5--cycle;
  z6 = (lbase+lspace,0);
  z7 = z6 + (h,0);
  z8 = z7 + (0,vthick);
  z9 = z6 + (vthick,vthick);
  z10= z9 + (0,espace);
  z11= z10+ (vthick,0);
  z12= z10+ (vthick,vthick);
  z13= z6 + (0,2*vthick+espace);
  pathe = z6--z7--z8--z9--z10--z11--z12--z13--cycle;
  z14= z7 + (dspace,0);
  z15= z14+ (dbase,0);
  z16= z15+ (halfh,halfh+dcorrection);
  z17= z15+ (dcorrection,h);
  z18= z14+ (0,h);
  z19= z18- (0,vthick);
  z20= z19+ (dbase,0);
  z21= z20+ (dinray,-dinray);
  z22= z15+ (0,vthick);
  z23= z14+ (0,vthick); 
  pathd = z14---z15..z16..tension tp..z17---z18--z19---z20..tension tp..z21..z22---z23--cycle;
  z24= z6 + (0,h);
  z25= z24+ (h,0);
  z26= z25- (0,vthick);
  z27= z24- (0,vthick);
  pathu = z24--z25--z26--z27--cycle;
%   fill pathl withcolor 0.78red;
%   fill pathe;
%   fill pathd;
%   fill pathu;
%  dotlabels( 14,15,16,17,18,19,20,21,22,23 );
  z100 = (0,-2h);
  drawoptions( shifted z100 );
  fill pathl withcolor 0.78red;
  fill pathe shifted (-z6) scaled scal shifted z200; 
  fill pathd shifted (-z6) scaled scal shifted z200; 
  fill pathu shifted (-z6) scaled scal shifted z200;
endfig;

end.
  
