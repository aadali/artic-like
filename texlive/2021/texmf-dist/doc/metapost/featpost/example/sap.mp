% sap.mp
% L. Nobre G.
% 2002

numeric NumberOfColors, HalfSideLength;
numeric ThicknesFactor, TolerancFactor;
NumberOfColors = 7;
HalfSideLength = 15mm;
ThicknesFactor = 0.2;         % Maximum TraceThickness over Half Side
ThinnessFactor = 0.1;         % Minimum TraceThickness over Half Side
TolerancFactor = 0.5 + 0.02;  % Side Gaps over Thickness

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
    aux = ceiling( uniformdeviate( MaxNum ) );
    if aux = 0:
      aux := 1;
    fi;
    ( aux )
  endgroup
enddef; 

def PickColor =
  ( uniformdeviate(1) , uniformdeviate(1), uniformdeviate(1) )
enddef;

vardef LocateIndex[] =
  dotlabel.urt( str @, z@ )
enddef;

def SpreadAllOver =
  begingroup
    pair urCorner, ulCorner, llCorner, lrCorner;
    path FinalCut;
    urCorner = ( 8.13in, 11.533in );
    llCorner = ( 0.118in, 35.5bp );
    ulCorner = ( xpart llCorner, ypart urCorner );
    lrCorner = ( xpart urCorner, ypart llCorner );
    FinalCut = urCorner--ulCorner--llCorner--lrCorner--cycle;
    
    picture Tile[];
    Tile[1] = currentpicture;
    Tile[2] = currentpicture rotated 90;
    Tile[3] = currentpicture rotated 180;
    Tile[4] = currentpicture rotated 270;
    currentpicture := nullpicture;
    numeric i, j, horizo, vertic;
    horizo = ( xpart urCorner ) + HalfSideLength; 
    vertic = ( ypart urCorner ) + HalfSideLength; 
    for i = 0 step 2HalfSideLength until horizo:
      for j = 0 step 2HalfSideLength until vertic:
	draw Tile[PickInteger( 4 )] shifted ( i, j );
      endfor;
    endfor;
    clip currentpicture to FinalCut
  endgroup
enddef;

beginfig(1);
  linecap := butt;
  
  numeric tracethickness;
  tracethickness = RandomStep;
  currentpen := pensquare rotated 45 scaled tracethickness;
 
  numeric i, currentcoord, maincoords[], numaincoords, toleranc;
  i = 0;
  currentcoord = 2TolerancFactor*tracethickness;
  toleranc = HalfSideLength - currentcoord;
  forever:
    i := incr( i );
    maincoords[i] = currentcoord;
    currentcoord := currentcoord + 4TolerancFactor*tracethickness;
    currentcoord := currentcoord +
                    uniformdeviate( ThicknesFactor*HalfSideLength );
    exitif currentcoord > toleranc;
  endfor;
  numaincoords = i;

  numeric j;
  pair Dir[];
  for i = 1 upto numaincoords:
    j := i;
    z[j] = ( HalfSideLength, maincoords[i] );
    Dir[j] = left;
    j := i + numaincoords;
    z[j] = ( maincoords[i], HalfSideLength );
    Dir[j] = down;
    j := i + 2numaincoords;
    z[j] = ( -maincoords[i], HalfSideLength );
    Dir[j] = down;
    j := i + 3numaincoords;
    z[j] = ( -HalfSideLength, maincoords[i] );
    Dir[j] = right;
    j := i + 4numaincoords;
    z[j] = ( -HalfSideLength, -maincoords[i] );
    Dir[j] = right;
    j := i + 5numaincoords;
    z[j] = ( -maincoords[i], -HalfSideLength );
    Dir[j] = up;
    j := i + 6numaincoords;
    z[j] = ( maincoords[i], -HalfSideLength );
    Dir[j] = up;
    j := i + 7numaincoords;
    z[j] = ( HalfSideLength, -maincoords[i] );
    Dir[j] = left;
  endfor;

%  for i=1 upto 8numaincoords:
%    LocateIndex[i];
%  endfor;

  boolean Filled[];
  for i=1 upto 8numaincoords:
    Filled[i] = false;
  endfor;
  numeric Chosen[];
  for i=0 upto numaincoords-1:
    Chosen[i] = 0;
  endfor;
  color LookUpColor[];
  for i=1 upto NumberOfColors:
    LookUpColor[i] = PickColor;
  endfor;
  
  numeric ThisInd, ThatInd, ThisCoo, ThatCoo;
  color ThisColor, ThatColor;
  path ActualPath;
  for j=2 step 2 until 8numaincoords:
    
    forever:
      ThisInd := PickInteger( 8numaincoords );
      exitunless Filled[ThisInd];
    endfor;
    Filled[ThisInd] := true;
    ThisCoo := ThisInd mod numaincoords;
    if Chosen[ThisCoo] = 0:
      i := PickInteger( NumberOfColors );
      ThisColor := LookUpColor[i];
      Chosen[ThisCoo] := i;
    else:
      ThisColor := LookUpColor[Chosen[ThisCoo]];
    fi;
    
    forever:
      ThatInd := PickInteger( 8numaincoords );
      exitunless Filled[ThatInd];
    endfor;
    Filled[ThatInd] := true;
    ThatCoo := ThatInd mod numaincoords; 
    if Chosen[ThatCoo] = 0:
      i := PickInteger( NumberOfColors );
      ThatColor := LookUpColor[i];
      Chosen[ThatCoo] := i;
    else:
      ThatColor := LookUpColor[Chosen[ThatCoo]];
    fi;

    ActualPath := z[ThisInd]{Dir[ThisInd]}
               ...z[ThatInd]{-Dir[ThatInd]};
    draw subpath (0,0.6) of ActualPath withcolor ThisColor;
    draw subpath (0.4,1) of ActualPath withcolor ThatColor;
    
  endfor;
  
  SpreadAllOver;
endfig;

end.
