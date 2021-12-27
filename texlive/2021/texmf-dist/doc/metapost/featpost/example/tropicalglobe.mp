% tropicalglobe.mp
% L. Nobre G. 
% 2012

input featpost3Dplus2D;

verbatimtex \documentclass{article}\begin{document} etex

beginfig(1);
  numeric radius, numc;
  pen thickp;
  
  f 		:= (3,0,2);
  viewcentr     := 0.5*(1,1,1);
  Spread 	:= 20;
  radius	= 1;
  numc		= 18;
  thickp	= pencircle scaled 1mm;
  
  drawoptions( withpen thickp );
  tropicalglobe( numc, black, radius, blue );
  produce_auto_scale;
endfig;

beginfig(2);
  color gammacnt, earthaxis;
  numeric newradius, aux, i, foc, lc;
  path conepath, latpath;
  picture theglobe;
    
  f 		:= (5,2,4);
  viewcentr     := black;
  Spread 	:= 70;
  radius	:= 1;
  numc		:= 9;
  thickp	:= pencircle scaled 0.5mm;
  
  earthaxis = eulerrotation( 0, 0, 23, blue);
  pickup thickp;
  tropicalglobe( numc, black, radius, earthaxis );
  for i=1 upto numc:				% longitudes
    aux		:= (i-1)*180/numc;
    gammacnt 	:= eulerrotation( 0, 0, 23, (cosd(aux),sind(aux),0) );
    draw spatialhalfcircle(black,gammacnt,radius,true);
  endfor;
  theglobe = currentpicture;
endfig;

beginfig(3);
  draw theglobe;
  TDAtiplen := 0.15;              % tdarrow and tdcircarrow
  TDAhalftipbase := 0.06;         % Three-Dimensional (Circular)
  TDAhalfthick := 0.03;           % Arrow
  path planet, twilight, shadowplanet, fullplanet;
  numeric endtime, begintime, cutmargin, scalevecs;
  color planetcolor, sunposition, earthvelocity, refpoi, surfacevel;
  color ua, ub, uc, ud, ue;
  pair va, vb, vc, vd, ve;
  cutmargin = 0.01;
  scalevecs = 1.3;
  planetcolor = 0.3*(255,152,55)/255;
  sunposition = red;
  earthvelocity = green;
  surfacevel = eulerrotation( 0, 0, 23, earthvelocity );
  refpoi = 1.1*sunposition;
  planet = rigorousfearpath( black, radius );
  twilight = spatialhalfcircle( black, sunposition, radius+cutmargin, true );
  endtime = xpart ( planet intersectiontimes twilight );
  begintime = xpart ( reverse planet intersectiontimes twilight );
  shadowplanet = subpath (-begintime-cutmargin,endtime+cutmargin) of planet;
  fullplanet = buildcycle( shadowplanet, reverse twilight );
  fill fullplanet withcolor planetcolor;
  tdcircarrow( black, earthaxis, 1.1, 0, 90 );
%  tdcircarrow( black, earthaxis, 1.1, 90, -90 );
  ua = refpoi+scalevecs*earthaxis;
  ub = refpoi+scalevecs*surfacevel;
  uc = refpoi-scalevecs*earthaxis;
  ud = refpoi-scalevecs*surfacevel;
  ue = 2.2*earthvelocity;
  tdarrow( 1.05*earthvelocity, ue );
  tdarrow( refpoi, uc );
  tdarrow( refpoi, ud );
  tdarrow( refpoi, ub );
  tdarrow( refpoi, ua );
  va = rp( ua );
  vb = rp( ub );
  vc = rp( uc );
  vd = rp( ud );
  ve = rp( ue );
  label.top( btex $\vec{A}$ etex, va );
  label.rt( btex $\vec{B}$ etex, vb );
  label.bot( btex $\vec{C}$ etex, vc );
  label.lft( btex $\vec{D}$ etex, vd );
  label.rt( btex $\vec{E}$ etex, ve );
endfig;
  
verbatimtex \end{document} etex

end.