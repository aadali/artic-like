% fieldlinesnormapaper.mp
% L. Nobre G.
% 2009

input featpost3Dplus2D;

Spread := 25;
f := (20,0,0);

color posa, posb, perpd;
numeric limy, limz, bord;
posa = (0,-3,3);  %%%%%%%% allways make Y(posa) < Y(posb) %%%%%%%%%%%
posb = (0,3,-2);  %%%%%%%%%%%% never make Y(posa) = Y(posb) %%%%%%%%%
perpd = (1,0,0);
limy = 7;
limz = 5;
bord = 2.5;

def vecfunc( expr posit )=
  begingroup
    save field;
    color field, horiz, verti, vech, vecv;
    numeric modula, modulb, ypos, zpos, hd, hf, vd, vf, zref, yref;
    modula = conorm( posit-posa );
    modulb = conorm( posit-posb );
    field = N( ( posit-posa )/(modula**3)-( posit-posb )/(modulb**3) );
    ypos = Y( posit );
    zpos = Z( posit );
    if Z(posa) <> Z(posb):
      (ypos,zref) = whatever[(Y(posa),Z(posa)),(Y(posb),Z(posb))];
      (yref,zpos) = whatever[(Y(posa),Z(posa)),(Y(posb),Z(posb))];
      if ( abs(zpos) > limz-bord ):
	hd = abs(zpos)-limz+bord;
	if ypos < yref:
	  if zpos > 0:
	    vech = -green;
	  else:
	    vech = green;
	  fi;
	else:
	  if zpos > 0:
	    vech = green;
	  else:
	    vech = -green;
	  fi;
	fi;
      else:
	hd = 0;
	vech = black;
      fi;
    else:
      zref = Z(posa);
      vech = green;
      if ( abs(zpos) > limz-bord ):
	hd = abs(zpos)-limz+bord;
      else:
	hd = 0;
      fi;
    fi;
    if ( abs(ypos) > limy-bord ):
      vd = abs(ypos)-limy+bord;
      if zpos > zref:
	if ypos < 0:
	  vecv = blue;
	else:
	  vecv = -blue;
	fi;
      else:
	if ypos < 0:
	  vecv = -blue;
	else:
	  vecv = blue;
	fi;
      fi;
    else:
      vd = 0;
      vecv = black;
    fi;
    hf = (hd/bord)**2;
    vf = (vd/bord)**2;
    ( N(field*(0.975-hf)*(0.975-vf)+hf*vech+vf*vecv) )
  endgroup
enddef;

def potfunc( expr posit )=
  begingroup
    color field;
    field = vecfunc( posit );
    ( ncrossprod( field, perpd ) )
  endgroup
enddef;
  
beginfig(1);
  numeric i, lena, lenb, numa, numb, sa, sb, diffstep, ray;
  numeric j, gridstep, fac;
  color stp, grdp, locv;
  path oneline;
  pen grossa, fina;
  grossa = pencircle scaled 3pt; 
  fina = pencircle scaled 2pt; 
  lena = 60;
  lenb = 40;
  numa = 15;
  numb = 30;
  diffstep = 0.1;
  ray = 0.5;
  sa = 40;
  sb = 5;
  fac = 0.275;
  gridstep = 0.3;
  for i=-limy+0.5*gridstep step gridstep until limy-0.5*gridstep:
    for j=-limz+0.5*gridstep step gridstep until limz-0.5*gridstep:
      grdp := (0,i,j);
      locv := 0.5*fac*vecfunc( grdp );
      drawarrow rp(grdp-locv)--rp(grdp+locv);
    endfor;
  endfor;
  for i=sa step numa until (360-sa):
    stp := posa+( 0, ray*cosd(i), ray*sind(i) );
    oneline := fieldlinepath( lena, stp, diffstep, vecfunc );
    draw oneline withpen grossa;
    draw oneline withpen fina withcolor background;
  endfor;
  for i=sb step numb until (360-sb):
    stp := posb+( 0, ray*cosd(i), ray*sind(i) );
    oneline := fieldlinepath( lenb, stp, -diffstep, vecfunc );
    draw oneline withpen grossa;
    draw oneline withpen fina withcolor background;
  endfor;
  drawoptions(withpen grossa);
  draw rp(posa);
  draw rp(posb);
  draw rp((0,-limy,-limz))--rp((0, limy,-limz))--
       rp((0, limy, limz))--rp((0,-limy, limz))--cycle;
  drawoptions();
endfig;

beginfig(2);
  numeric i, lena, numa, diffstep, fac, gridstep, margnum;
  color stp, grdp, locv;
  path oneline, twoline, bordeline;
  pen grossa, fina;
  grossa = pencircle scaled 3pt; 
  fina = pencircle scaled 2pt; 
  lena = 85;
  numa = 13;
  margnum = 2;
  diffstep = 0.2;
  fac = 0.275;
  gridstep = 0.3;
  for i=-limy+0.5*gridstep step gridstep until limy-0.5*gridstep:
    for j=-limz+0.5*gridstep step gridstep until limz-0.5*gridstep:
      grdp := (0,i,j);
      locv := 0.5*fac*potfunc( grdp );
      drawarrow rp(grdp-locv)--rp(grdp+locv);
    endfor;
  endfor;
  for i=margnum upto numa-margnum:
    stp := (i/numa)[posa,posb];
    oneline := fieldlinepath( lena, stp, diffstep, potfunc );
    twoline := fieldlinepath( lena, stp, -diffstep, potfunc );
    draw oneline withpen grossa;
    draw twoline withpen grossa;
  endfor;
  drawoptions(withpen grossa);
  draw rp(posa);
  draw rp(posb);
  bordeline = rp((0,-limy,-limz))--rp((0, limy,-limz))--
  rp((0, limy, limz))--rp((0,-limy, limz))--cycle;
  clip currentpicture to bordeline;
  draw bordeline;
  drawoptions(withpen fina);
  numeric side;
  path aux, electri;
  side = 10mm;
  aux = origin--(side,0)--(side,side)--cycle;
  z0 = 0.667[origin,(side,0.5*side)];
  electri = aux shifted (-z0) rotated (45+angle(rp(posb)-rp(posa)));
  fill electri shifted rp(posb);
  fill electri rotated 180 shifted rp(posa);
endfig;

end.
