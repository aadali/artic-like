% torus.mp
% L. Nobre G. 
% 2014

input featpost3Dplus2D;

% Draw a torus with NP dots, a section of 3.14*RS^2
% and a total diameter of 2*(RB+RS). The dots are placed 
% on a single closed line that turns around the $z$ axis
% NB times and around the section NS times. 
% figure(1) exemplifies the use of signalvertex.

def flask( expr TheVal ) =
  begingroup
    numeric RB, RS, RP, NB, NS, NM, phi, theta, first, second, third; 
    RB=1;
    RS=0.3;
    RP=0.02;
    NB=8;
    NS=5;
    NM=4;
    theta=360*NB*TheVal;
    phi=360*NS*TheVal;
    first =(RB+(RS+RP*sind(NM*theta))*cosd(phi))*cosd(theta);
    second=(RB+(RS+RP*sind(NM*theta))*cosd(phi))*sind(theta);
    third =    (RS+RP*sind(NM*theta))*sind(phi);
    ( (first,second,third) )
  endgroup
enddef;

def blask( expr TheVal ) =
  begingroup
    numeric RB, RS, NB, NS, phi, theta, first, second, third; 
    RB=1;
    RS=0.1;
    NB=1;
    NS=24;
    theta=360*NB*TheVal;
    phi=360*NS*TheVal;
    first =(RB+RS*cosd(phi))*cosd(theta);
    second=(RB+RS*cosd(phi))*sind(theta);
    third =    RS*sind(phi);
    ( (first,second,third) )
  endgroup
enddef;

beginfig(5);
  Spread := 200;
  numeric NP, i, varfrac;
  path cl;
  string comm;
  pen bip, sp;
  varfrac = 1.25;
  bip = pencircle scaled 15pt;
  sp = pencircle scaled 9pt;
  NP=600;
  cl = for i=1 upto NP: rp(flask(i/NP)).. endfor cycle;
  draw cl;
  for i=1 upto NP:
    comm:="draw subpath ("
      & decimal(i-0.5)
      & ","
      & decimal(i+0.5)
      & ") of cl withpen bip; undraw subpath ("
      & decimal(i-varfrac)
      & ","
      & decimal(i+varfrac)
      & ") of cl withpen sp;";
    getready( comm, flask( i/NP ) );
  endfor;
  closedline( true, 300, 0.5, 1 )( blask );
  doitnow;
endfig;

beginfig(4);
  f := (16,2,3.8); 
  Spread:=40;
  pickup pencircle scaled 3pt;
  smoothtorus( black-3*blue, blue, 4, 1 );
  smoothtorus( black, blue, 4, 1 );
  smoothtorus( black+3*blue, blue, 4, 1 );
endfig;

beginfig(3);
  Spread := 100;
  f := (7,3,9);
  HoriZon := -2;
  LightSource := (1,-1,15);
  numeric sinthe, costhe, ux, uz, lx, lz, co, i;
  costhe = (4*sqrt(5)-sqrt(3))/11;
  sinthe = sqrt(1-(costhe**2));
  ux = costhe*sqrt(2)/2;
  uz = sinthe/2;
  lx = sqrt(4/3)*sinthe;
  lz = -sqrt(2/3)*costhe;
  co = sqrt(2)*(2*sqrt(5)+5*sqrt(3))/22;
  V1 := (ux,0,-uz); 
  V2 := (lx,0,-lz);
  V3 := (0,ux,uz); V15:= (0,-ux,uz);
  V4 := (0,lx,lz); V16:= (0,-lx,lz);
  V5 := (co,co,0); V14:= (-co,-co,0);
  V6 := (-ux,0,-uz);
  V7 := (-lx,0,-lz);
  V8 := (-co,co,0); V13:= (co,-co,0);
  V9 := improvertex( V5, 1, V4, 1, V3, sqrt(3)-1, (0,1,0) );
  V10:= (-X(V9),-Y(V9),Z(V9));
  V11:= (Y(V9),-X(V9),-Z(V9));
  V12:= (-Y(V9),X(V9),-Z(V9));
  makeface1(1,2,3);
  makeface2(1,2,5);
  makeface3(2,3,5);
  makeface4(9,4,5);
  makeface5(3,4,1);
  makeface6(1,4,5);
  makeface7(3,6,7);
  makeface8(3,7,8);
  makeface9(6,7,8);
  makeface10(9,8,4);
  makeface11(4,6,8);
  makeface12(3,4,6); 
  makeface13(3,5,9);
  makeface14(3,8,9);
  makeface15(7,8,12);
  makeface16(7,14,12);
  makeface17(6,8,12);
  makeface18(6,14,12);
  makeface19(2,5,11);
  makeface20(2,13,11);
  makeface21(1,5,11);
  makeface22(1,13,11);
  makeface23(10,15,13);
  makeface24(10,15,14); 
  makeface25(10,16,13);
  makeface26(10,16,14);
  makeface27(2,1,15);
  makeface28(6,7,15);
  makeface29(13,15,2);
  makeface30(14,15,7);
  makeface31(16,1,15);
  makeface32(16,6,15);
  for i=1 upto 32:
    FC[i] := TableColors;
    FCD[i]:= true;
  endfor;
  ShadowOn := true;
  setthearena(15,6);
  fill_faces();
  draw_invisible( false, false, TableC3, TableC0 );
%   drawoptions( withcolor red );
%   for i=1 upto 16:
%     label.bot( decimal( i ), rp( V[i] ) );
%   endfor;
endfig;

% Draw a torus with faces.
                  
beginfig(2);

    PrintStep := 7;

    f := ( 3,1,1 );

    RB:=1;
    RS:=0.25;
    NB:=3;
    NS:=5;
    jB:=360/NB;
    jS:=360/NS;
    iS:=jS;
    iB:=jB;
    counter:=0;
    for i=1 upto NB:
        for j=1 upto NS:
	    counter := incr(counter);
	    npf[counter]:= 4;
            phi  :=j*jS;
            theta:=i*jB;
            f1   :=(RB+(RS)*cosd(phi))*cosd(theta);
            s1   :=(RB+(RS)*cosd(phi))*sind(theta);
            t1   :=    (RS)*sind(phi);
	    F[counter]p[1]:=(f1,s1,t1);
            f2   :=(RB+(RS)*cosd(phi))*cosd(theta+iB);
            s2   :=(RB+(RS)*cosd(phi))*sind(theta+iB);
            t2   :=    (RS)*sind(phi);
	    F[counter]p[2]:=(f2,s2,t2);
            f3   :=(RB+(RS)*cosd(phi+iS))*cosd(theta+iB);
            s3   :=(RB+(RS)*cosd(phi+iS))*sind(theta+iB);
            t3   :=    (RS)*sind(phi+iS);
	    F[counter]p[3]:=(f3,s3,t3);
            f4   :=(RB+(RS)*cosd(phi+iS))*cosd(theta);
            s4   :=(RB+(RS)*cosd(phi+iS))*sind(theta);
            t4   :=    (RS)*sind(phi+iS);
	    F[counter]p[4]:=(f4,s4,t4);             
        endfor;
    endfor;
    NF:=counter;
%    faceraytrace(0.3, blue);
%    draw_all_test(false);
    sharpraytrace;
endfig;

beginfig(1);
  color tmpcolor;
  string tmpstr;
  Nobjects := 0;
  f := ( 3,1,1 );
  Spread := 110;
  NP:=624;
  RB:=1;
  RS:=0.2;
  R3:=0.1;
  NB:=8;
  NS:=5;
  jB:=360*NB/NP;
  jS:=360*NS/NP;
  for i=1 upto NP:
    phi   :=i*jS;
    theta :=i*jB;
    first :=(RB+(RS+R3*sind(4*theta))*cosd(phi))*cosd(theta);
    second:=(RB+(RS+R3*sind(4*theta))*cosd(phi))*sind(theta);
    third :=    (RS+R3*sind(4*theta))*sind(phi);
    tmpcolor := (first,second,third);
    tmpstr := cstr( tmpcolor );
    getready("signalvertex(" & tmpstr &
      ",1,black);signalvertex(" & tmpstr & ",0.7,white);", tmpcolor );
  endfor;
  doitnow;
endfig;

beginfig(6);
  f := (6,2,4);
%  f := (6,2,7);      % NO CUSP.
  Spread:=70;
  HoriZon := -2;
  ShadowOn := true;
  numeric bigray, smaray, angstep;
  color toruscenter, torusmoment, nearaxe, sideaxe, viewline;
  color circlecenter, circlemoment;
  numeric ang, ind, i, anglim;
  path cpath, apath, opath, ipath, wp, ep;
  pair outerp[], innerp[], refpair;
  pen markpen;
  boolean cuspcond;
  picture holepic;
  markpen = pencircle scaled 5pt;
  pickup markpen;
  bigray = 2;
  smaray = 1;
  angstep= 2;
  toruscenter = black;
  torusmoment = (0,1,2);
  refpair = unitvector( rp(toruscenter+torusmoment)-rp(toruscenter) );
  viewline = f-toruscenter;
  ind = 0;
  sideaxe = bigray*ncrossprod( torusmoment, viewline );
  nearaxe = bigray*ncrossprod( sideaxe, torusmoment );
  anglim = 0.5*angstep-180;
  if ShadowOn:
    for ang=anglim step angstep until -anglim:
      circlecenter:= nearaxe*cosd(ang)+sideaxe*sind(ang);
      circlemoment:=-nearaxe*sind(ang)+sideaxe*cosd(ang);
      filldraw circleshadowpath(circlecenter,circlemoment,smaray);
    endfor;
  fi;
  for ang=anglim step angstep until -anglim:
    ind := incr( ind );
    circlecenter:= nearaxe*cosd(ang)+sideaxe*sind(ang);
    circlemoment:=-nearaxe*sind(ang)+sideaxe*cosd(ang);
    cpath:=spatialhalfcircle(circlecenter,circlemoment,smaray,true);
    if ang>=0:
      outerp[ind]=point 0 of cpath;
      innerp[ind]=point (length cpath) of cpath;
    elseif ang<0:
      innerp[ind]=point 0 of cpath;
      outerp[ind]=point (length cpath) of cpath;
    fi;
  endfor;
  opath = for i=1 upto ind: outerp[i].. endfor cycle;
  ipath = for i=1 upto ind: innerp[i].. endfor cycle;
  holepic = currentpicture;
  clip holepic to ipath;
  unfill opath;
  draw holepic;
  draw opath withcolor green;
  draw ipath withcolor red;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% There is an analytic way of getting the angle of the cusp point! %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  i := ceiling(1+180/angstep);
  forever:
    i := incr( i );
    cuspcond := refpair dotprod innerp[i+1] < refpair dotprod innerp[i];
    exitif cuspcond or ( i > ind-1 );
  endfor;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ep = outerp[ind-i]--innerp[ind-i];
  wp = innerp[i]--outerp[i];
  unfill buildcycle(opath,ep,reverse ipath,wp);
  draw opath withcolor green;
  draw subpath (i,ind-i) of ipath withcolor red;
endfig;

beginfig(7); %%%%%%% copy of figure 6
  f := (6,2,4);
  Spread:=70;
  HoriZon := -2;
  ShadowOn := true;
  numeric bigray, smaray;
  color toruscenter, torusmoment;
  pen markpen;
  markpen = pencircle scaled 5pt;
  bigray = 2;
  smaray = 1;
  toruscenter = black;
  torusmoment = (0,1,2);
  setthestage( 17, 8 );
  pickup markpen;
  smoothtorus( toruscenter, torusmoment, bigray, smaray );
endfig;

end.


