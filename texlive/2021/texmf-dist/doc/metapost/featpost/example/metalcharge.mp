% metalcharge.mp
% L. Nobre G.
% IYP (2005)

input featpost3Dplus2D;

def parabolline( expr ParC, FracofLen ) =
  begingroup
    numeric Xlengh, Xcontx, Xconty, Xpos, Xa, Xb, Xyfin;
    Xlengh = X( ParC );
    Xcontx = Y( ParC );
    Xconty = Z( ParC );
    Xpos = Xlengh*FracofLen;
    if FracofLen < Xcontx/Xlengh :
      Xa = Xconty*(2/Xcontx - 1/(Xcontx-Xlengh));
      Xb = Xconty*(Xcontx/(Xcontx-Xlengh)-1)/(Xcontx**2);
      Xyfin = Xa*Xpos+Xb*(Xpos**2);
    else:
      Xyfin = Xconty*((Xpos-Xcontx)/(Xcontx-Xlengh)+1);
    fi;
    ( Xyfin )
  endgroup
enddef;
    
beginfig(1);
  numeric i, j, k, ang, lengh, contx, conty, resol, pot[], lim, resola, u;
  numeric ringcharg[], rinray[], lstep, coun, sumc, dist, sump, fact, sca;
  color actpoi, scanpoi, dirv;
  pair shiv;
  path form, chargedensity, simetr, revolenvel;
  lengh = 1.6;
  contx = 0.5;
%  contx = 0.9;
  conty = 0.75;
  fact = 0.3;
  lim = 30;
  resol = 30;
  sca = 0.17;
  u = 6cm;
  shiv = (4cm, 12cm);
  resola = 12; %make it even
  viewcentr := (0,contx,0);
  lstep = lengh/resol;
  coun = 0;
  for i = 1 upto resol:
    rinray[i] = parabolline( (lengh, contx, conty), (i-0.5)/resol );
    ringcharg[i] = 1;
  endfor;
  forever:
    coun := incr( coun );
    sump := 0;
    for k=1 upto resol:
      actpoi := green*lengh*(k-0.5)/resol + blue*rinray[k];
      pot[k] := 0;
      for i=1 upto resol:
	for j=1 upto resola:
	  ang := 360*(j-0.5)/resola-180;
	  scanpoi:=green*lengh*(i-0.5)/resol+rinray[i]*(sind(ang),0,cosd(ang));
	  dirv := actpoi - scanpoi;
	  dist := 50*conorm( dirv );
	  pot[k] := pot[k] + ringcharg[i]/dist;
	endfor;
      endfor;
      sump := sump + pot[k];
    endfor;
    sump := sump/resol;
    sumc := 0;
    for k=1 upto resol:
      ringcharg[k] := ringcharg[k]*(1+fact*(sump-pot[k])/sump);
      sumc := sumc + ringcharg[k];
    endfor;
    for k=1 upto resol:
      ringcharg[k] := ringcharg[k]*resol/sumc;
    endfor;
    exitif coun>lim;
    show ringcharg[resol];
  endfor;
  for i=1 upto resol:
    ringcharg[i] := ringcharg[i]/rinray[i];
  endfor;
  chargedensity = ((lstep*u,sca*ringcharg1*u)
         for i=2 upto resol:
	   --(i*lstep*u,sca*ringcharg[i]*u)
	 endfor) shifted shiv;
  draw chargedensity;
  form = ((lstep*u,rinray1*u)
                  for i=2 upto resol:
		    ..(i*lstep*u,rinray[i]*u)
		  endfor) shifted shiv;
  simetr = ((lstep*u,-rinray1*u)
                  for i=2 upto resol:
		    ..(i*lstep*u,-rinray[i]*u)
		  endfor) shifted shiv;
  revolenvel = form...(reverse simetr)...cycle; 
  draw revolenvel withcolor red;
endfig;

end.