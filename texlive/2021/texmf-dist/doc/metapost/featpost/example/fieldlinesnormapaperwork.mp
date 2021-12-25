% fieldlinesnormapaperwork.mp
% L. Nobre G.
% 2013

input featpost3Dplus2D;

prologues := 1;

Spread := 25;
f := (20,0,0);

color posa, posb, perpd;
numeric limy, limz, bord, reduct;
posa = (0,-4,0);  %%%%%%%% allways make Y(posa) < Y(posb) %%%%%%%%%%%
posb = (0,4,0);  %%%%%%%%%%%% never make Y(posa) = Y(posb) %%%%%%%%%
perpd = (1,0,0);
limy = 7;
limz = 5;
bord = 2.5;
reduct = 0.18;
  
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
    ( N(field*(0.975-hf)*(0.975-vf)+(hf*vech+vf*vecv)*reduct) )
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
  numeric i, lena, numa, diffstep, ray;
  color stp;
  path bordeline;
  pen grossa, fina;
  pair auxp;
  grossa = pencircle scaled 3pt; 
  fina = pencircle scaled 1pt; 
  lena = 80;
  numa = 13;
  diffstep = 0.12;
  ray = 0.05;
  % drawoptions(withpen grossa);
  % for i=1 upto numa-1:
  %   auxp := (i/numa)[rp(posa),rp(posb)];
  %   draw auxp;
  % endfor;
  drawoptions(withpen fina);
  for i=2 step 2 until numa-2:
    stp := (i/numa)[posa,posb];
    draw fieldlinepath( lena, stp, diffstep, potfunc );
    draw fieldlinepath( lena, stp, -diffstep, potfunc );
  endfor;
%  draw fieldlinepath( lena, (0,Y(posa),Z(posa)+ray), diffstep, vecfunc );
  bordeline = rp((0,-limy,-limz))--rp((0, limy,-limz))--
  rp((0, limy, limz))--rp((0,-limy, limz))--cycle;
  clip currentpicture to bordeline;
  draw bordeline;
  drawoptions(withpen fina);
  % numeric side;
  % path aux, electri;
  % side = 10mm;
  % aux = origin--(side,0)--(side,side)--cycle;
  % z0 = 0.667[origin,(side,0.5*side)];
  % electri = aux shifted (-z0) rotated (45+angle(rp(posb)-rp(posa)));
  % draw electri shifted rp(posb);
  % draw electri rotated 180 shifted rp(posa);
  path electri;
  electri = fullcircle scaled 7mm;
  fill electri shifted rp(posb) withcolor 0.5white;
  fill electri shifted rp(posa) withcolor 0.5white;
endfig;

end.
