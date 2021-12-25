% qap.mp
% L. Nobre G.
% IYP (2005)

def PickColor =
  ( uniformdeviate(1) , uniformdeviate(1), uniformdeviate(1) )
enddef;

numeric HalfSideLength, ThinnessFactor;
numeric ThicknesFactor, BorderThFactor;

HalfSideLength = 10mm;
ThicknesFactor = 0.17;         % Maximum TraceThickness over Half Side
ThinnessFactor = 0.12;         % Minimum TraceThickness over Half Side
BorderThFactor = 1.35;        % Borders over TraceThickness

def RandomStep =
  begingroup
    numeric margin;
    margin = ThicknesFactor-ThinnessFactor;
    ThinnessFactor*HalfSideLength + uniformdeviate( margin*HalfSideLength )
  endgroup
enddef;

def PickInteger( expr MaxNum ) =
  begingroup
    save aux;
    numeric aux;
    aux = floor( uniformdeviate( MaxNum+0.995 ) );
    ( aux )
  endgroup
enddef; 

vardef LocateIndex[] =
  dotlabel.urt( str @, z@ )
enddef;

beginfig(1);
  numeric QNum;
  QNum = 5 + PickInteger( 5 );
  numeric tracethickness, i, j, k, angadd, aa, ax, ay;
  tracethickness = RandomStep;
  pen tracepen, bordepen;
  color TheColor;
  tracepen = pencircle scaled tracethickness;
  bordepen = tracepen scaled BorderThFactor;
  path qpath;
  aa = 360/QNum;
  ax = HalfSideLength*2*cosd(0.5*aa);
  ay = HalfSideLength*2*sind(0.5*aa);
  z0 = (ax,0);
  z1 = (0,ay);
  z2 = (-ax,0);
  z3 = (0,-ay);
  z4 = z0;
  qpath = z0--z1--z2--z3--cycle;
  angadd = uniformdeviate( 20 );
  TheColor = PickColor;  
  pair Dir[], pos[];
  boolean Filled[];
  for i = 0 upto 11:
    j := floor( i/3 );
    Dir[i] = dir(angle(z[j+1]-z[j])+80+angadd);
    Filled[i] = false;
    pos[i] = (0.25*(i-3*j+1))[z[j],z[j+1]];
  endfor;
%  for i=0 upto 5:
%    LocateIndex[i];
%  endfor;
  numeric ThisInd, ThatInd;
  path ActualPath;
  for i=1 step 2 until 11:    
    forever:
      ThisInd := PickInteger( 11 );
      exitunless Filled[ThisInd];
    endfor;
    Filled[ThisInd] := true;
    forever:
      ThatInd := PickInteger( 11 );
      exitunless Filled[ThatInd];
    endfor;
    Filled[ThatInd] := true;
    ActualPath := pos[ThisInd]{Dir[ThisInd]}...pos[ThatInd]{-Dir[ThatInd]};
%    draw ActualPath withcolor background withpen bordepen;
    draw ActualPath withcolor TheColor withpen tracepen;
  endfor;
%  draw qpath;
%   clip currentpicture to qpath;
  pair urCorner, ulCorner, llCorner, lrCorner, cent;
  path FinalCut;
  urCorner = ( 8.13in, 11.533in );
  llCorner = ( 0.118in, 35.5bp );
  ulCorner = ( xpart llCorner, ypart urCorner );
  lrCorner = ( xpart urCorner, ypart llCorner );
  FinalCut = urCorner--ulCorner--llCorner--lrCorner--cycle;
  cent = 0.25*(urCorner+ulCorner+llCorner+lrCorner);
  picture Tile[];
  for i=0 upto 1:
    Tile[i] = currentpicture rotated (i*180);
  endfor;
  currentpicture := nullpicture;
  pair pairh, pairv;
  for k = 1 upto QNum:
    pairh := HalfSideLength*2*dir((k-1)*aa);
    pairv := HalfSideLength*2*dir(k*aa);
    for i = 0 upto 11:
      for j = 0 upto 11:
 	draw Tile[PickInteger( 1 )]
 	rotated ((k-0.5)*aa)
 	shifted ( cent + ax*dir((k-0.5)*aa) + i*pairh + j*pairv );
      endfor;
    endfor;
  endfor;
  clip currentpicture to FinalCut;
endfig;

end.
