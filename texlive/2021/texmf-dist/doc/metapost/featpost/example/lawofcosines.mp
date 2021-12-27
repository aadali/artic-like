% lawofcosines.mp
% L. Nobre G.
% 2007

prologues := 1;

beginfig(1);
  numeric u, limit, major, minor, abvcos, blwcos;
  color gridc;
  pen gridpmaj, gridpmin, dotp;
  u=0.7mm;
  limit=64;
  major=20;
  minor=2;
  abvcos=0.5;
  blwcos=1/3;
  gridc=0.5white;
  gridpmaj=pencircle scaled 2pt;
  gridpmin=pencircle scaled 1.3pt;
  dotp=pencircle scaled 2.7pt;
%  labeloffset := 3u;
  numeric i,j, abvsin, blwsin, stpmaj, stpmin;
  abvsin=1 +-+ abvcos;
  blwsin=1 +-+ blwcos;
  stpmaj=major*u;
  stpmin=minor*u;
  pair abvdir, blwdir, thedir;
  abvdir=(abvcos,abvsin);
  blwdir=(blwcos,-blwsin);
  boolean ismajor, isaxe;
  drawoptions( withpen gridpmin withcolor gridc );
  for i=-limit step minor until limit:
    ismajor := not ((i mod major)<>0);
    if ismajor:
      isaxe := not (i <> 0);
      if isaxe:
	drawoptions( withpen gridpmaj );
      else:
	drawoptions( withpen gridpmaj withcolor gridc );
      fi;
    fi;
    draw ((i*u,0)+limit*u*blwdir)--(i*u,0)--((i*u,0)+limit*u*abvdir);    
    if i<0:
      thedir := -blwdir;
    else:
      thedir := abvdir;
    fi;
    draw ((-limit*u,0)+i*u*thedir)--((limit*u,0)+i*u*thedir);
    if ismajor:
      drawoptions( );
      label.bot( decimal( abs( i ) ), (i*u,0)+limit*u*blwdir );
      if i<0:
	label.llft( decimal( abs( i ) ), (-limit*u,0)+i*u*thedir );
      elseif not (i<>0):
	label.lft( decimal( abs( i ) ), (-limit*u,0)+i*u*thedir );
      else:
	label.ulft( decimal( abs( i ) ), (-limit*u,0)+i*u*thedir );
      fi;
      label.top( decimal( abs( i ) ), (i*u,0)+limit*u*abvdir );
      drawoptions( withpen gridpmin withcolor gridc );
    fi;    
  endfor;
  drawoptions( withpen dotp withcolor red );
  numeric signal;
  for i=-limit upto limit:
    if i <> 0:
      if i<0:
	signal := -1;
      else:
	signal := 1;
      fi;
      for j=-limit upto -2:
	if not ( ( sqrt((i**2)+(j**2)+signal*2*abs(i*j)/3) mod 1 ) <> 0 ):
	  draw (i*u,0)-j*u*blwdir;
	fi;
      endfor;
      for j=limit downto 2:
	if not ( ( sqrt((i**2)+(j**2)+signal*2*abs(i*j)*abvcos) mod 1 ) <> 0 ):
	  draw (i*u,0)+j*u*abvdir;
	fi;
      endfor;
    fi;
  endfor;
  
endfig;

beginfig(2);
  numeric u, limit, major, minor, abvcos, blwcos;
  color gridc;
  pen gridpmaj, gridpmin, dotp;
  u=0.7mm;
  limit=24;
  major=20;
  minor=2;
  abvcos=0.8;
  blwcos=0.9;
  gridc=0.5white;
  gridpmaj=pencircle scaled 1.3pt;
  gridpmin=pencircle scaled 0.7pt;
  dotp=pencircle scaled 1.4pt;
%  labeloffset := 3u;
  numeric i,j, abvsin, blwsin, stpmaj, stpmin;
  abvsin=1 +-+ abvcos;
  blwsin=1 +-+ blwcos;
  stpmaj=major*u;
  stpmin=minor*u;
  pair abvdir, blwdir, thedir;
  abvdir=(abvcos,abvsin);
  blwdir=(blwcos,-blwsin);
  boolean ismajor, isaxe;
  drawoptions( withpen gridpmin withcolor gridc );
  for i=-limit step minor until limit:
    ismajor := not ((i mod major)<>0);
    if ismajor:
      isaxe := not (i <> 0);
      if isaxe:
	drawoptions( withpen gridpmaj );
      else:
	drawoptions( withpen gridpmaj withcolor gridc );
      fi;
    fi;
    draw ((i*u,0)+limit*u*blwdir)--(i*u,0)--((i*u,0)+limit*u*abvdir);    
    if i<0:
      thedir := -blwdir;
    else:
      thedir := abvdir;
    fi;
    draw ((-limit*u,0)+i*u*thedir)--((limit*u,0)+i*u*thedir);
    if ismajor:
      drawoptions( );
      label.bot( decimal( abs( i ) ), (i*u,0)+limit*u*blwdir );
      if i<0:
	label.llft( decimal( abs( i ) ), (-limit*u,0)+i*u*thedir );
      elseif not (i<>0):
	label.lft( decimal( abs( i ) ), (-limit*u,0)+i*u*thedir );
      else:
	label.ulft( decimal( abs( i ) ), (-limit*u,0)+i*u*thedir );
      fi;
      label.top( decimal( abs( i ) ), (i*u,0)+limit*u*abvdir );
      drawoptions( withpen gridpmin withcolor gridc );
    fi;    
  endfor;
  drawoptions( withpen dotp withcolor red );
  numeric signal;
  for i=-limit upto limit:
    if i <> 0:
      if i<0:
	signal := -1;
      else:
	signal := 1;
      fi;
      for j=-limit upto -2:
	if not ( ( sqrt((i**2)+(j**2)+signal*18*abs(i*j)/10) mod 1 ) <> 0 ):
	  draw (i*u,0)-j*u*blwdir;
	fi;
      endfor;
      for j=limit downto 2:
	if not ( ( sqrt((i**2)+(j**2)+signal*8*abs(i*j)/5) mod 1 ) <> 0 ):
	  draw (i*u,0)+j*u*abvdir;
	fi;
      endfor;
    fi;
  endfor;
  draw origin withcolor blue withpen dotp;
endfig;

end.