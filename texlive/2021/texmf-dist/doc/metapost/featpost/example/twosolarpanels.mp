% twosolarpanels.mp
% L. Nobre G.
% 2014

input featpost3Dplus2D;

f := (2.7,5,2);
Spread := 170;  
color poi[], colthr, coltwo;
pen penthr, pentwo;
numeric projx, projy, projz;
penthr = pencircle scaled 2.4pt;
pentwo = pencircle scaled 1.0pt;
colthr = blue;
coltwo = black;
projx = -0.5;
projy = -1;
projz = -0.5;

def drawonesegmentmanytimes( expr poia, poib ) =
  begingroup
    draw rp((projx,Y(poia),Z(poia)))--rp((projx,Y(poib),Z(poib))) withpen pentwo withcolor coltwo;
    draw rp((X(poia),projy,Z(poia)))--rp((X(poib),projy,Z(poib))) withpen pentwo withcolor coltwo;
    draw rp((X(poia),Y(poia),projz))--rp((X(poib),Y(poib),projz)) withpen pentwo withcolor coltwo;
    draw rp(poia)--rp(poib) withpen penthr withcolor colthr
  endgroup
enddef;

vardef segments( text t ) =
  color tmpcol[];
  numeric tmpnum, i;
  tmpnum := 0;
  forsuffixes $=t:
    tmpnum := tmpnum+1;
    tmpcol[tmpnum] = poi$;
  endfor;
  for i=1 upto tmpnum-1:
    drawonesegmentmanytimes(tmpcol[i],tmpcol[i+1]);
  endfor
enddef;

beginfig(1);
  numeric fac;
  fac = 0.5*sqrt(2);
  numeric panelangle, roofangle, panelwid, panellen;
  panelangle = 61;
  roofangle = 22;
  panelwid = 0.5;
  panellen = 0.8;
  numeric panelz, panelx;
  panelz = panellen*sind(panelangle);
  panelx = panellen*cosd(panelangle);
  numeric forwarddescent, backwardascent, middledescent;
  forwarddescent = panelwid*fac*sind(roofangle)/cosd(roofangle);
  backwardascent = panelx*fac*sind(roofangle)/cosd(roofangle);
  middledescent = (panelwid-panelx)*fac*sind(roofangle)/cosd(roofangle);
  poi1 = (panelwid,0,panelz);
  poi2 = (panelwid,panelwid,panelz);
  poi3 = (0,panelwid,panelz);
  poi4 = (panelwid,0,backwardascent);
  poi5 = (0,panelwid,backwardascent);
  poi6 = (panelwid+panelx,0,0);
  poi7 = (panelwid+panelx,panelwid,0);
  poi8 = (panelwid,panelwid,0);
  poi9 = (panelwid,panelwid+panelx,0);
  poi10= (0,panelwid+panelx,0);
  poi11= (panelwid,panelwid,-middledescent);
  poi12= (panelwid+panelx,panelwid,-forwarddescent);
  poi13= (panelwid,panelx+panelwid,-forwarddescent);
  segments( 6, 12, 11, 13, 10 );
  segments( 12, 7, 8, 9, 13 );
  segments( 7, 6, 1, 2, 7 );
  segments( 9, 10, 3, 2, 9 );
  segments( 11, 8, 2 );
  segments( 6, 4, 1 );
  segments( 10, 5, 3 );
  segments( 1, 3 );
  segments( 12, 4, 5, 13 );
  segments( 6, 10 );
  segments( 4, 11, 5 );
  unfill rp(poi7)--rp(poi6)--rp(poi1)--rp(poi2)--cycle;
  unfill rp(poi9)--rp(poi10)--rp(poi3)--rp(poi2)--cycle;
endfig;

end.