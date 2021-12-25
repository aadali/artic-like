% electricpotential.mp
% L. Nobre G.
% 2004
% PAY ATTENTION: THESE FIGURES ARE INDEPENDANT. CONDUCTOR SURFACE POTENTIAL
% CHANGES UNDER THE INFLUENCE OF EXTERNAL CHARGES.

input featpost3Dplus2D;

f := 16*(4,1,1);
Spread := 35;
LightSource := 10*(4,-3,4);
ActuC := 2;

numeric potlim, cac, cbc, ccc, cdc;
potlim = 2;
cac = 2;
cbc = -2;
ccc = 4;
cdc = 1.5;
pair cap, cbp, ccp;
cap = (-2,-4);
ccp = (-2,4);
cbp = (2,0);
numeric np, ssize;
np = 56;
ssize = 18;

def zi( expr xc, yc ) =
  0
enddef;

def za( expr xc, yc ) =
  begingroup
    save out;
    numeric dist, out;
    dist = abs( (xc,yc)-cap );
    if abs(dist) < 0.005:
      if cac<0:
	out = -potlim;
      else:
	out = potlim;
      fi;
    else:
      out = cac/dist;
    fi;
    if out > potlim:
      out := potlim;
    elseif out < -potlim:
      out := -potlim;
    fi;
    (out)
  endgroup
enddef;

def zb( expr xc, yc ) =
  begingroup
    save out;
    numeric dist, out;
    dist = abs( (xc,yc)-cbp );
    if abs(dist) < 0.005:
      if cbc<0:
	out = -potlim;
      else:
	out = potlim;
      fi;
    else:
      out = cbc/dist + za( xc, yc );
    fi;
    if out > potlim:
      out := potlim;
    elseif out < -potlim:
      out := -potlim;
    fi;
    (out)
  endgroup
enddef;

def zc( expr xc, yc ) =
  begingroup
    save out;
    numeric dist, out;
    dist = abs( (xc,yc)-ccp );
    if abs(dist) < 0.005:
      if ccc<0:
	out = -potlim;
      else:
	out = potlim;
      fi;
    else:
      out = ccc/dist + zb( xc, yc );
    fi;
    if out > potlim:
      out := potlim;
    elseif out < -potlim:
      out := -potlim;
    fi;
    (out)
  endgroup
enddef;

def zd( expr xc, yc ) =
  begingroup
    save out;
    numeric dist, out, i, tmp;
    pair poi;
    tmp = 0;
    for i=0 upto np:
      poi := (i/np)[cap,ccp];
      dist := abs( (xc,yc)-poi );
      if abs(dist) < 0.005:
	if cdc<0:
	  out := -potlim;
	else:
	  out := potlim;
	fi;
      else:
	out := cdc/dist/np;
      fi;
      tmp := tmp+out;
    endfor;
    tmp := tmp + zc( xc, yc );
    if tmp > potlim:
      tmp := potlim;
    elseif tmp < -potlim:
      tmp := -potlim;
    fi;
    (tmp)
  endgroup
enddef;

beginfig(0); 
  hexagonaltrimesh( false, np, ssize, zi );
endfig;

beginfig(1); 
  hexagonaltrimesh( false, np, ssize, za );
endfig;

beginfig(2); 
  hexagonaltrimesh( false, np, ssize, zb );
endfig;

beginfig(3); 
  hexagonaltrimesh( false, np, ssize, zc );
endfig;

% Going into "MetaPost capacity exceeded!"
%
%beginfig(4); 
%  hexagonaltrimesh( false, np, ssize, zd );
%endfig;

end;
