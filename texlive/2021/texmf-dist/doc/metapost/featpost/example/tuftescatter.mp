% tuftescatter.mp
% L. Nobre G.
% IYP (2005)

input featpost3Dplus2D;

def linscatter = uniformdeviate( 2 ) - 1 enddef;

def fullcalculator =
  begingroup
    color aux;
    numeric cx, cy;
    forever:
      cx := linscatter;
      cy := linscatter;
      exitif cx ++ cy < 0.8;               % Saddle ray 
    endfor;
    aux = ( cx, cy, cx*cx-cy*cy );
    ( eulerrotation( -25, -15, -5, aux ) )  % Rotation angles
  endgroup
enddef;

beginfig(1);
  numeric totalnum, classnum, sparksiz, pointray;
  totalnum = 250;                          % Number of data points
  classnum = 30;                           % Number of sparkline classes
  sparksiz = 0.34;                         % Maximum sparkline height
  pointray = 0.035;                        % Size of the data points
  f := ( 7, 6, 5 );
  Spread := 120;
  numeric i, j, k, l, hx[], hy[], hz[], classtep;
  classtep = 2/classnum;
  color tmp;
  for j=0 upto classnum:
    hx[j] = 0;
    hy[j] = 0;
    hz[j] = 0;
  endfor;
  for l=1 upto totalnum:
    tmp := fullcalculator;
    for k=1 upto classnum:
      if ( X(tmp) <= k*classtep-1 ) and ( X(tmp) > (k-1)*classtep-1 ):
	hx[k] := incr( hx[k] );
      fi;
      if ( Y(tmp) <= k*classtep-1 ) and ( Y(tmp) > (k-1)*classtep-1 ):
	hy[k] := incr( hy[k] );
      fi;
      if ( Z(tmp) <= k*classtep-1 ) and ( Z(tmp) > (k-1)*classtep-1 ):
	hz[k] := incr( hz[k] );
      fi;
    endfor;
    draw rigorouscircle( (-1,Y(tmp),Z(tmp)), red, pointray );
    draw rigorouscircle( (X(tmp),-1,Z(tmp)), green, pointray );
    draw rigorouscircle( (X(tmp),Y(tmp),-1), blue, pointray );
    getready("fill rigorousfearpath("
	& cstr( tmp )
	& ",pointray) withcolor blue; draw rigorousfearpath("
	& cstr( tmp )
	& ",pointray); ", tmp);
  endfor;
  numeric maxn, sparkfac;
  maxn = 0;
  for i=1 upto classnum:
    if hx[i]>maxn: maxn := hx[i]; fi;
    if hy[i]>maxn: maxn := hy[i]; fi;
    if hz[i]>maxn: maxn := hz[i]; fi;
  endfor;
  sparkfac = sparksiz/maxn;
  for j=1 upto classnum:
    draw rp(((j-1)*classtep-1,1+sparkfac*hx[j-1],-1))
           --rp(((j-1)*classtep-1,1+sparkfac*hx[j],-1))
           --rp((j*classtep-1,1+sparkfac*hx[j],-1));
    draw rp((-1,(j-1)*classtep-1,1+sparkfac*hy[j-1]))
           --rp((-1,(j-1)*classtep-1,1+sparkfac*hy[j]))
           --rp((-1,j*classtep-1,1+sparkfac*hy[j]));
    draw rp((-1,1+sparkfac*hz[j-1],(j-1)*classtep-1))
           --rp((-1,1+sparkfac*hz[j],(j-1)*classtep-1))
           --rp((-1,1+sparkfac*hz[j],j*classtep-1));
  endfor;
  doitnow;
endfig;

end.
