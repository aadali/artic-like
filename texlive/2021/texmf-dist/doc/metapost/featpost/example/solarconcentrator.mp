% solarconcentrator.mp
% L. Nobre G.
% 2013

prologues := 1;

input featpost3Dplus2D;

color actpoi, actdir, actperp;
boolean finish;

%def findintersection( expr SourcePoi, RayDir ) =
%  begingroup
%    
%  endgroup
%enddef;

vardef reflectrayr(expr sourcepoi, raydir )=
  save J, Perp, newdir, cropro, tandir, theangle;
  color J, Perp, newdir, cropro, tandir;
  numeric theangle;
  findintersection( sourcepoi, raydir );
  if not finish:
    draw rp(sourcepoi)--rp(actpoi);
    theangle = 180-getangle( raydir, actperp );
    if theangle<0.5: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DANGER!
      cropro = ncrossprod( raydir, actperp );
      tandir = ncrossprod( actperp, cropro );
      newdir = planarrotation( actperp, tandir, theangle );
    else:
      newdir = actperp;
    fi;
    reflectrayr( actpoi, newdir );
  else:
    J=sourcepoi+raydir;
    drawarrow rp(sourcepoi)--rp(J);
    finish := false;
  fi;
enddef;

beginfig(1);
  f := (6,4,5);
  numeric iz, jang, numcircs, zcoord;
  numcircs = 10;
  draw rp((1,0,1))--rp(red)--rp(-red)--rp((-1,0,1));
  for iz=1 upto numcircs:
    zcoord := iz/numcircs;
    draw
      rp((1-zcoord,zcoord,zcoord))--rp((zcoord-1,zcoord,zcoord))..
      for jang=5 step 5 until 175:
      rp((zcoord-1-zcoord*sind(jang),zcoord*cosd(jang),zcoord))..
      endfor
      rp((zcoord-1,-zcoord,zcoord))--rp((1-zcoord,-zcoord,zcoord))..
      for jang=5 step 5 until 175:
      rp((1-zcoord+zcoord*sind(jang),-zcoord*cosd(jang),zcoord))..
      endfor
      cycle;
  endfor;
  draw rp(red)--rp((0,1,1))--rp(-red)--rp((0,-1,1))--cycle;    
endfig;



end.