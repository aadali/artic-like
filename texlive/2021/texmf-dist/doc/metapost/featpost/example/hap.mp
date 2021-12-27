% hap.mp
% L. Nobre G.
% IYP (2005)

numeric HalfSideLength;
numeric ThicknesFactor, BorderThFactor;

HalfSideLength = 6mm;
ThicknesFactor = 0.5;
BorderThFactor = 1.85;        % Borders over TraceThickness

def PickInteger( expr MaxNum ) = floor( uniformdeviate( MaxNum+1 ) ) enddef; 

def PickColor =
  begingroup
    numeric c[], i, min, max, ref;
    ref = 0.76;
    for i=1 upto 3:
      c[i] = uniformdeviate(1);
    endfor;
    min = c1;
    max = c1;
    for i=2 upto 3:
      if c[i] < min: min := c[i]; fi;
      if c[i] > max: max := c[i]; fi;
    endfor;
    for i=1 upto 3:
      c[i] := ref+(0.94-ref)*(c[i]-min)/(max-min);
    endfor;
    ( ( c1 , c2, c3 ) )
  endgroup
enddef;

vardef LocateIndex[] =
  dotlabel.urt( str @, z@ )
enddef;

beginfig(1);
  numeric tracethickness, i, j, hexconst, angaux, cou, nupersi;
  numeric ThisInd, ThatInd, dist[], ang[];
  color TheColor;
  pen tracepen, bordepen;
  path hexpath, ActualPath, FinalCut;
  pair corn[], nco, aco, vec, Dir[], urCorner, ulCorner, llCorner, lrCorner;
  boolean Filled[];
  picture Tile[];
  hexconst = sqrt(3)/2;
  nupersi = PickInteger( 2 );
  tracethickness = ThicknesFactor*
                   (0.5+uniformdeviate(0.5))*
                   HalfSideLength/(nupersi+1);
  TheColor = PickColor;
  tracepen = pencircle scaled tracethickness;
  bordepen = tracepen scaled BorderThFactor;
  for i=0 upto 6:
      corn[i] = 2*HalfSideLength*dir(i*60);
  endfor;
  hexpath = corn0--corn1--corn2--corn3--corn4--corn5--cycle;
  dist0 = 0;
  ang0 = 40+uniformdeviate( 100 );
  for i=1 upto nupersi:
    dist[i] = dist[i-1]+(0.7+uniformdeviate(0.45))*HalfSideLength/(nupersi+1);
    angaux := 40+20*i/nupersi;
    ang[i] = angaux + uniformdeviate( 180-2*angaux );
  endfor;
  for i=0 upto 5:
    aco := corn[i];
    nco := corn[1+i];
    vec := nco-aco;
    z[i] = aco+0.5*vec;
    Dir[i] = dir(angle(vec)+ang0);
    Filled[i] = false;
  endfor;
  cou = 5;
  for i=0 upto 5:
    aco := corn[i];
    nco := corn[i+1];
    vec := unitvector(nco-aco);
    for j=1 upto nupersi:
      cou := incr( cou );
      z[cou] = 0.5[aco,nco]+dist[j]*vec;
      Dir[cou] = dir(angle(vec)+ang[j]);
      Filled[cou] = false;
      cou := incr( cou );
      z[cou] = 0.5[aco,nco]-dist[j]*vec;
      Dir[cou] = dir(angle(vec)+ang[j]);
      Filled[cou] = false;
    endfor;
  endfor;
%  for i=0 upto cou:
%    LocateIndex[i];
%  endfor;
  for i=1 step 2 until cou:    
    forever:
      ThisInd := PickInteger( cou );
      exitunless Filled[ThisInd];
    endfor;
    Filled[ThisInd] := true;
    forever:
      ThatInd := PickInteger( cou );
      exitunless Filled[ThatInd];
    endfor;
    Filled[ThatInd] := true;
    ActualPath := z[ThisInd]{Dir[ThisInd]}
      ...z[ThatInd]{-Dir[ThatInd]};
    draw ActualPath withcolor background withpen bordepen;
    draw ActualPath withcolor TheColor withpen tracepen;
  endfor;
  clip currentpicture to hexpath;
  urCorner = ( 8.13in, 11.533in );
  llCorner = ( 0.118in, 35.5bp );
  ulCorner = ( xpart llCorner, ypart urCorner );
  lrCorner = ( xpart urCorner, ypart llCorner );
  FinalCut = urCorner--ulCorner--llCorner--lrCorner--cycle;
  for i=0 upto 5:
    Tile[i] = currentpicture rotated (i*60);
  endfor;
  currentpicture := nullpicture;
  numeric j, horizo, vertic, vertigo;
  horizo = ( xpart urCorner ) + 2HalfSideLength; 
  vertic = ( ypart urCorner ) + 2HalfSideLength;
  vertigo = 0.5*horizo;
  pair paird, pairv;
  pairv = up*4*HalfSideLength*hexconst;
  paird = dir(30)*4*HalfSideLength*hexconst;
  numeric numd, numv;
  numv = (vertic+2*vertigo)/(4*HalfSideLength*hexconst);
  numd = horizo/(3*HalfSideLength);
  for i = 0 upto numd:
    for j = 0 upto numv:
      draw Tile[PickInteger( 5 )] shifted ( (0,-vertigo) + i*paird + j*pairv );
    endfor;
  endfor;
  clip currentpicture to FinalCut;
endfig;

end.
