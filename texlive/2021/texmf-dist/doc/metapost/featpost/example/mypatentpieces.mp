% mypatentpieces.mp
% L. Nobre G.
% IYP (2005)

input featpost3Dplus2D;

prologues := 1;

vardef P@# = rp(V@#) enddef;

%f := 2*(17,7,9);
f := 2*(21,5,8);
Spread := 36;
background := white;

numeric u, i;
u = 0.1;
V1 := (0,0,0);
V2 := (0,0,4u);
V3 := (0,16u,4u);
V4 := (4u,16u,4u);
V5 := (4u,16u,8u);
V6 := (0,16u,8u);
V7 := (0,24u,16u);
V8 := (4u,24u,16u);
V9 := (4u,35u,4u);
V10:= (4u,35u,0);
V11:= (16u,23u,0);
V12:= (16u,23u,4u);
V13:= (16u,0,4u);
V14:= (16u,0,0);
V15:= (0,43u,8u);
V16:= (4u,43u,8u);
V17:= (8u,8u,0);
V18:= (4u,43u,16u);
V19:= (0,43u,16u);

V22:= (-62.5u+16u,0,0);
V23:= (62.5u-16u,0,0);
V24:= (-62.5u,0,0);
V25:= (62.5u,0,0);
V26:= (-62.5u,0,4u);
V27:= (62.5u,0,4u);
V28:= (-62.5u+16u,0,4u);
V29:= (62.5u-16u,0,4u);

V40:= (95.5u,-43u,-8u); % zero
V42:= (-62.5u+21u,-43u,8u); % 200
V43:= (-62.5u+17u,-8u,-43u); %300
V44:= (95.5u,8u,-43u); % 400
V55:= (-32u,-4u,-35u); % 500
V56:= (82u,-35.25u,-4u); % 600 %%%
V57:= (50u,0,-35u); % 700
V58:= (0,-35.13u,0); % 800 %%%

def indicp( expr Ind ) =
  ( rp(V[Ind])--rp( (X(V[Ind]),Y(V[Ind]),40u) ) )
enddef;

def indicl( expr Ind ) =
  ( rp((X(V[Ind]),Y(V[Ind]),40u)) )
enddef;

% 100
for i=1 upto 19:
  V[i+100] := (X(V[i]),-Y(V[i]),Z(V[i]));
endfor;

% 200
for i=101 upto 119:
  V[i+100] := (X(V[i]),-Y(V[i]),-Z(V[i]));
endfor;

% 300
for i=101 upto 119:
  V[i+200] := (X(V[i]),Z(V[i]),-Y(V[i]));
endfor;

% 400
for i=301 upto 319:
  V[i+100] := (-X(V[i]),-Y(V[i]),Z(V[i]));
endfor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% zero
for i=1 upto 19:
  V[i] := V[i]+V40;
endfor;

% 200
for i=201 upto 219:
  V[i] := V[i]+V42;
endfor;

% 300
for i=301 upto 319:
  V[i] := V[i]+V43;
endfor;

% 400
for i=401 upto 419:
  V[i] := V[i]+V44;
endfor;

% 500
for i=22 upto 29:
  V[i+500] := (X(V[i]),Z(V[i]),-Y(V[i]))+V55;
endfor;

% 600
for i=22 upto 29:
  V[i+600] := V[i]+V56;
endfor;

% 700
for i=22 upto 29:
  V[i+700] := (X(V[i]),Z(V[i]),-Y(V[i]))+V57;
endfor;

% 800
for i=22 upto 29:
  V[i+800] := V[i]+V58;
endfor;

numeric dhip, dcat, dhei, dang;
dhip = X(V627)-X(V824);
dcat = X(V25)-X(V24);
dhei = dcat +-+ (0.5*dhip);
dang = angle(0.5*dhip,dhei);

% 1000
for i=822 upto 829:
  V[i+100] := (X(V[i])*cosd(dang)-Y(V[i])*sind(dang),
    X(V[i])*sind(dang)+Y(V[i])*cosd(dang),Z(V[i]));
endfor;
for i=922 upto 929:
  V[i+100] := V[i]+V627-V925;
endfor;

% 900
for i=622 upto 629:
  V[i-500] := (X(V[i])*cosd(dang)+Y(V[i])*sind(dang),
    -X(V[i])*sind(dang)+Y(V[i])*cosd(dang),Z(V[i]));
endfor;
for i=122 upto 129:
  V[i+800] := V[i]+V824-V126;
endfor;

beginfig(1);
  pen penedge;
  penedge= pencircle scaled 1.5pt;
  pickup penedge;
  path refillpat[], retouchpat[];

  path pat[];
  pat1 = P309--P312--P313--P314--P311--P310--P309--P304--P305
    --P308--P307;
  pat2 = P311--P312;
  pat11 = P313--P302--P303;
  pat3 = pat11 cutafter pat1;
  pat10 = P301-- for i=-90 step 5 until 85:
    rp(V315+(0,8u*sind(i),8u*cosd(i)))..
    endfor P319--P307;
  pat4 = P310-- for i=-90 step 5 until 85:
    rp(V316+(0,8u*sind(i),8u*cosd(i)))..
    endfor P318--P308;
  pair auxp,auxP;
  auxp = point 0 of spatialhalfcircle(V315,red,8.05u,false);
  auxP = point 0 of spatialhalfcircle(V316,red,8.05u,false);
  pat9 = P315--auxp;  
  pat5 = pat10 cutbefore pat9;
  pat6 = auxp--auxP;
  for i=1 upto 6:
    draw pat[i];
  endfor;
  fakehole(V315,red*4u, 2u );
  fakehole(V317,green*4u, 2u );
  V388 := 0.7[V308,V318];
  V379 := 0.7[V307,V319];
  refillpat2 = P309--P304--P305--P308--P307--P379--P388--cycle;
  retouchpat3 = P309--P304--P305--P308--P307--P379;
  retouchpat4 = P388--P308;

  
  path pat[];
  pat1 = for i=-90 step 2 until 90:
    rp(V522+(-2u*cosd(i),0,2u*sind(i)))...
    endfor
    for i=-90 step 2 until 90:
    rp(V523+(2u*cosd(i),0,2u*sind(-i)))...
    endfor
    cycle;
  pat2 = for i=-90 step 2 until 90:
    rp(V528+(-2u*cosd(i),0,2u*sind(i)))...
    endfor
    for i=-90 step 2 until 90:
    rp(V529+(2u*cosd(i),0,2u*sind(-i)))...
    endfor
    cycle;
  pat3 = for i=-90 step 2 until 90:
    rp(V524+(-8u*cosd(i),0,8u*sind(i)))...
    endfor
    for i=-90 step 2 until 90:
    rp(V525+(8u*cosd(i),0,8u*sind(-i)))...
    endfor
    cycle;
  pat4 = for i=-90 step 2 until 90:
    rp(V526+(-8u*cosd(i),0,8u*sind(i)))...
    endfor
    for i=-90 step 2 until 90:
    rp(V527+(8u*cosd(i),0,8u*sind(-i)))...
    endfor
    cycle;
  pair auxp, auxP, auxm, auxM;
  auxp = point 0 of spatialhalfcircle(V524,green,8.05u,true);
  auxP = point 0 of spatialhalfcircle(V526,green,8.05u,true);
  auxm = point 0 of spatialhalfcircle(V525,green,8.05u,false);
  auxM = point 0 of spatialhalfcircle(V527,green,8.05u,false);
  pat5 = P524--auxp;
  pat6 = P525--auxm;
  pat7 = pat3 cutbefore pat5 cutafter pat6;
  pat8 = auxp--auxP;
  pat9 = auxm--auxM;

  unfill pat3;
%  draw pat1;
  draw pat2;
  draw pat7;
  draw pat4;
  draw pat8;
  draw pat9;
  fakehole(V524,green*4u, 2u );
  fakehole(V525,green*4u, 2u );


  path pat[];
  pat1 = for i=-90 step 2 until 90:
    rp(V722+(-2u*cosd(i),0,2u*sind(i)))...
    endfor
    for i=-90 step 2 until 90:
    rp(V723+(2u*cosd(i),0,2u*sind(-i)))...
    endfor
    cycle;
  pat2 = for i=-90 step 2 until 90:
    rp(V728+(-2u*cosd(i),0,2u*sind(i)))...
    endfor
    for i=-90 step 2 until 90:
    rp(V729+(2u*cosd(i),0,2u*sind(-i)))...
    endfor
    cycle;
  pat3 = for i=-90 step 2 until 90:
    rp(V724+(-8u*cosd(i),0,8u*sind(i)))...
    endfor
    for i=-90 step 2 until 90:
    rp(V725+(8u*cosd(i),0,8u*sind(-i)))...
    endfor
    cycle;
  pat4 = for i=-90 step 2 until 90:
    rp(V726+(-8u*cosd(i),0,8u*sind(i)))...
    endfor
    for i=-90 step 2 until 90:
    rp(V727+(8u*cosd(i),0,8u*sind(-i)))...
    endfor
    cycle;
  pair auxp, auxP, auxm, auxM;
  auxp = point 0 of spatialhalfcircle(V724,green,8.05u,true);
  auxP = point 0 of spatialhalfcircle(V726,green,8.05u,true);
  auxm = point 0 of spatialhalfcircle(V725,green,8.05u,false);
  auxM = point 0 of spatialhalfcircle(V727,green,8.05u,false);
  pat5 = P724--auxp;
  pat6 = P725--auxm;
  pat7 = pat3 cutbefore pat5 cutafter pat6;
  pat8 = auxp--auxP;
  pat9 = auxm--auxM;
  pat10= pat4 cutbefore P727--auxM cutafter P726--auxP;
  unfill pat7--pat9--pat10--(reverse pat8)--cycle;
%  draw pat1;
  draw pat2;
  draw pat7;
  draw pat4;
  draw pat8;
  draw pat9;
  fakehole(V724,green*4u, 2u );
  fakehole(V725,green*4u, 2u );

  path pat[];
  pat1 = P401-- for i=90 step -5 until -85:
    rp(V415+(0,8u*sind(i),8u*cosd(i)))..
    endfor P419--P407--P406--P403--P402--cycle;
  pat10 = P401--P414--P411--P410-- for i=90 step -5 until -85:
    rp(V416+(0,8u*sind(i),8u*cosd(i)))..
    endfor P418;
  pair auxp,auxP;
  auxp = point 0 of spatialhalfcircle(V415,red,8.05u,false);
  auxP = point 0 of spatialhalfcircle(V416,red,8.05u,false);
  pat9 = P416--auxP;  
  pat2 = pat10 cutafter pat9;
  pat3 = auxp--auxP;
  unfill pat1;
  unfill pat10--P403--cycle;
  for i=1 upto 3:
    draw pat[i];
  endfor;
  fakehole(V415,-red*4u, 2u );
  fakehole(V417,-green*4u, 2u );


  path pat[];
  pat1 = P3--P2--P13--P14--P11--P10--P9--P12--P13;
  pat2 = P11--P12;
  pat3 = P9--P4--P5--P6--P7--P8--P5;
  numeric i;
  pat8 = P1-- for i=-90 step 5 until 85:
    rp(V15+(0,8u*cosd(i),8u*sind(i)))..
    endfor P19--P7;
  pat5 = P10-- for i=-90 step 5 until 85:
    rp(V16+(0,8u*cosd(i),8u*sind(i)))..
    endfor P18--P8;
  pat6 = pat1 cutbefore pat3;
  pair auxp,auxP;
  auxp = point 0 of spatialhalfcircle(V15,red,8.05u,true);
  auxP = point 0 of spatialhalfcircle(V16,red,8.05u,true);
  pat9 = rp(V15)--auxp;  
  pat4 = pat8 cutbefore pat9;
  pat7 = auxp--auxP;
  unfill buildcycle( pat1, pat5, pat8, reverse pat3 );
  for i=2 upto 7:
    draw pat[i];
  endfor;
  fakehole(V15,red*4u, 2u );
  fakehole(V17,blue*4u, 2u );
  refillpat1 = P9--P4--P5--P6--P7--P19--P18--cycle;
  retouchpat1 = P9--P4--P5--P6--P7--P19;
  retouchpat2 = P18--P8--P5;
  retouchpat5 = P8--P7;
  
   path pat[];
   pat2 = for i=-90+dang step 2 until 90+dang:
     rp(V928+(-2u*cosd(i),2u*sind(i),0))...
     endfor
     for i=-90+dang step 2 until 90+dang:
     rp(V929+(2u*cosd(i),2u*sind(-i),0))...
     endfor
     cycle;
   pat3 = for i=-90+dang step 2 until 90+dang:
     rp(V924+(-8u*cosd(i),8u*sind(i),0))...
     endfor
     for i=-90+dang step 2 until 90+dang:
     rp(V925+(8u*cosd(i),8u*sind(-i),0))...
     endfor
     cycle;
   pat4 = for i=-90+dang step 2 until 90+dang:
     rp(V926+(-8u*cosd(i),8u*sind(i),0))...
     endfor
     for i=-90+dang step 2 until 90+dang:
     rp(V927+(8u*cosd(i),8u*sind(-i),0))...
     endfor
     cycle;
   pair auxp, auxP, auxm, auxM;
   auxp = point 0 of spatialhalfcircle(V924,blue,8.05u,false);
   auxP = point 0 of spatialhalfcircle(V926,blue,8.05u,false);
   auxm = point 0 of spatialhalfcircle(V925,blue,8.05u,true);
   auxM = point 0 of spatialhalfcircle(V927,blue,8.05u,true);
   pat5 = P924--auxp;
   pat6 = P925--auxm;
   pat7 = pat3 cutbefore pat5 cutafter pat6;
   pat8 = auxp--auxP;
   pat9 = auxm--auxM;
   pat10= pat4 cutbefore P927--auxM cutafter P926--auxP;
   unfill pat7--pat8--pat10--(reverse pat9)--cycle;

   draw pat2;
   draw pat7;
   draw pat4;
   draw pat8;
   draw pat9;
   fakehole(V924,blue*4u, 2u );
   fakehole(V925,blue*4u, 2u );
  
  path pat[];
  pat1 = for i=-90 step 2 until 90:
    rp(V622+(-2u*cosd(i),2u*sind(i),0))...
    endfor
    for i=-90 step 2 until 90:
    rp(V623+(2u*cosd(i),2u*sind(-i),0))...
    endfor
    cycle;
  pat2 = for i=-90 step 2 until 90:
    rp(V628+(-2u*cosd(i),2u*sind(i),0))...
    endfor
    for i=-90 step 2 until 90:
    rp(V629+(2u*cosd(i),2u*sind(-i),0))...
    endfor
    cycle;
  pat3 = for i=-90 step 2 until 90:
    rp(V624+(-8u*cosd(i),8u*sind(i),0))...
    endfor
    for i=-90 step 2 until 90:
    rp(V625+(8u*cosd(i),8u*sind(-i),0))...
    endfor
    cycle;
  pat4 = for i=-90 step 2 until 90:
    rp(V626+(-8u*cosd(i),8u*sind(i),0))...
    endfor
    for i=-90 step 2 until 90:
    rp(V627+(8u*cosd(i),8u*sind(-i),0))...
    endfor
    cycle;
  pair auxp, auxP, auxm, auxM;
  auxp = point 0 of spatialhalfcircle(V624,blue,8.05u,false);
  auxP = point 0 of spatialhalfcircle(V626,blue,8.05u,false);
  auxm = point 0 of spatialhalfcircle(V625,blue,8.05u,true);
  auxM = point 0 of spatialhalfcircle(V627,blue,8.05u,true);
  pat5 = P624--auxp;
  pat6 = P625--auxm;
  pat7 = pat3 cutbefore pat5 cutafter pat6;
  pat8 = auxp--auxP;
  pat9 = auxm--auxM;
  unfill pat3;

%%  pat10 = buildcycle( pat1, pat2 );
%%  draw pat10;
%%  draw pat1 cutbefore pat2;
%% THIS IS BAD, VERY BAD!!!! BUILDCYCLE DOES NOT WORK!!!!
%  show length pat1;
%  show pat1 intersectiontimes pat2;
%  show reverse pat1 intersectiontimes pat2;
  numeric auxta, auxtb;
  auxta = xpart (reverse pat1 intersectiontimes pat2);
  auxtb = xpart (pat1 intersectiontimes pat2);
  draw subpath (-auxta,auxtb) of pat1;

  draw pat2;
  draw pat7;
  draw pat4;
  draw pat8;
  draw pat9;
  fakehole(V624,blue*4u, 2u );
  fakehole(V625,blue*4u, 2u );


   path pat[];
   pat1 = for i=-90 step 2 until 90:
     rp(V822+(-2u*cosd(i),2u*sind(i),0))...
     endfor
     for i=-90 step 2 until 90:
     rp(V823+(2u*cosd(i),2u*sind(-i),0))...
     endfor
     cycle;
   pat2 = for i=-90 step 2 until 90:
     rp(V828+(-2u*cosd(i),2u*sind(i),0))...
     endfor
     for i=-90 step 2 until 90:
     rp(V829+(2u*cosd(i),2u*sind(-i),0))...
     endfor
     cycle;
   pat3 = for i=-90 step 2 until 90:
     rp(V824+(-8u*cosd(i),8u*sind(i),0))...
     endfor
     for i=-90 step 2 until 90:
     rp(V825+(8u*cosd(i),8u*sind(-i),0))...
     endfor
     cycle;
   pat4 = for i=-90 step 2 until 90:
     rp(V826+(-8u*cosd(i),8u*sind(i),0))...
     endfor
     for i=-90 step 2 until 90:
     rp(V827+(8u*cosd(i),8u*sind(-i),0))...
     endfor
     cycle;
   pair auxp, auxP, auxm, auxM;
   auxp = point 0 of spatialhalfcircle(V824,blue,8.05u,false);
   auxP = point 0 of spatialhalfcircle(V826,blue,8.05u,false);
   auxm = point 0 of spatialhalfcircle(V825,blue,8.05u,true);
   auxM = point 0 of spatialhalfcircle(V827,blue,8.05u,true);
   pat5 = P824--auxp;
   pat6 = P825--auxm;
   pat7 = pat3 cutbefore pat5 cutafter pat6;
   pat8 = auxp--auxP;
   pat9 = auxm--auxM;
   pat10= pat4 cutbefore P827--auxM cutafter P826--auxP;
   unfill pat7--pat9--pat10--(reverse pat8)--cycle;

 %%  pat10 = buildcycle( pat1, pat2 );
 %%  draw pat10;
 %%  draw pat1 cutbefore pat2;
 %% THIS IS BAD, VERY BAD!!!! BUILDCYCLE DOES NOT WORK!!!!
 %  show length pat1;
 %  show pat1 intersectiontimes pat2;
 %  show reverse pat1 intersectiontimes pat2;
   numeric auxta, auxtb;
   auxta = xpart (reverse pat1 intersectiontimes pat2);
   auxtb = xpart (pat1 intersectiontimes pat2);
   draw subpath (-auxta,auxtb) of pat1;

   draw pat2;
   draw pat7;
   draw pat4;
   draw pat8;
   draw pat9;
   fakehole(V824,blue*4u, 2u );
   fakehole(V825,blue*4u, 2u );


   path pat[];
   pat1 = P209--P212--P213--P214--P211;
   pat10 = P211--P210-- for i=90 step -5 until -85:
     rp(V216+(0,8u*cosd(i),8u*sind(i)))..
     endfor P218--P208--P205;
   pat2 = pat10 cutafter pat1;
   pat11 = P214--P201-- for i=90 step -5 until -85:
     rp(V215+(0,8u*cosd(i),8u*sind(i)))..
     endfor P219--P207;
   pair auxp,auxP;
   auxp = point 0 of spatialhalfcircle(V215,red,8.05u,true);
   auxP = point 0 of spatialhalfcircle(V216,red,8.05u,true);
   pat9 = P215--auxp;  
   pat3 = pat11 cutafter pat9;
   pat4 = auxp--auxP;
   pat6 = P211--P212;
   pat5 = P210--P209;
   unfill buildcycle(pat1,pat11,pat10);
   for i=1 upto 6:
     draw pat[i];
   endfor;
   fakehole(V215,red*4u, 2u );
   fakehole(V217,-blue*4u, 2u );

   unfill refillpat1;
   draw retouchpat1;
   draw retouchpat2;
   unfill refillpat2;
   draw retouchpat3;
   draw retouchpat4;
   draw retouchpat5;

   unfill buildcycle(pat1,pat11,pat10);
   for i=1 upto 6:
     draw pat[i];
   endfor;
   fakehole(V215,red*4u, 2u );
   fakehole(V217,-blue*4u, 2u );

 %   label.top( "40", indicl(40) );
 %   label.top( "42", indicl(42) );
 %   label.top( "43", indicl(43) );
 %   label.top( "44", indicl(44) );
 %   label.top( "55", indicl(55) );
 %   label.top( "56", indicl(56) );
 %   label.top( "57", indicl(57) );
 %   label.top( "58", indicl(58) );

 %   draw indicp(40);
 %   draw indicp(42);
 %   draw indicp(43);
 %   draw indicp(44);
 %   draw indicp(55);
 %   draw indicp(56);
 %   draw indicp(57);
 %   draw indicp(58);

    path pat[];
    pat1 = for i=-90-dang step 2 until 90-dang-2:
      rp(V1022+(-2u*cosd(i),2u*sind(i),0))...
      endfor rp(V1022+(-2u*cosd(90-dang),2u*sind(90-dang),0));
    pat11 = for i=-90-dang step 2 until 90-dang-2:
      rp(V1023+(2u*cosd(i),2u*sind(-i),0))...
      endfor rp(V1023+(2u*cosd(90-dang),2u*sind(-90+dang),0));
    pat2 = for i=-90-dang step 2 until 90-dang:
      rp(V1028+(-2u*cosd(i),2u*sind(i),0))...
      endfor
      for i=-90-dang step 2 until 90-dang:
      rp(V1029+(2u*cosd(i),2u*sind(-i),0))...
      endfor
      cycle;
    pat3 = for i=-90-dang step 2 until 90-dang:
      rp(V1024+(-8u*cosd(i),8u*sind(i),0))...
      endfor
      for i=-90-dang step 2 until 90-dang:
      rp(V1025+(8u*cosd(i),8u*sind(-i),0))...
      endfor
      cycle;
    pat4 = for i=-90-dang step 2 until 90-dang:
      rp(V1026+(-8u*cosd(i),8u*sind(i),0))...
      endfor
      for i=-90-dang step 2 until 90-dang:
      rp(V1027+(8u*cosd(i),8u*sind(-i),0))...
      endfor
      cycle;
    pair auxp, auxP, auxm, auxM;
    auxp = point 0 of spatialhalfcircle(V1024,blue,8.05u,true);
    auxP = point 0 of spatialhalfcircle(V1026,blue,8.05u,true);
    auxm = point 0 of spatialhalfcircle(V1025,blue,8.05u,false);
    auxM = point 0 of spatialhalfcircle(V1027,blue,8.05u,false);
    pat5 = P1024--auxp;
    pat6 = P1025--auxm;
    pat7 = pat3 cutbefore pat6 cutafter pat5;
    pat8 = auxp--auxP;
    pat9 = auxm--auxM;
    pat10= pat4 cutbefore P1027--auxM cutafter P1026--auxP;
    unfill pat7--pat8--pat10--(reverse pat9)--cycle;

    draw (point (length pat1) of pat1 -- point 0 of pat11) cutafter pat2;
    draw pat1 cutbefore pat2;
    draw pat2;
    draw pat7;
    draw pat4;
    draw pat8;
    draw pat9;
    fakehole(V1024,blue*4u, 2u );
    fakehole(V1025,blue*4u, 2u );




endfig;

end.















