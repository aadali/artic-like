% torusexperiment.mp
% L. Nobre G. 
% 2014

% DO NOT CHANGE PARAMETERS. YOU HAVE BEEN WARNED!

prologues := 1;
input featpost3Dplus2D;
%SphericalDistortion := true;
%ParallelProj := true;
beginfig(1);
  Spread:=70;
  numeric bigray, smaray, angstep;
  color toruscenter, torusmoment, nearaxe, sideaxe, viewline;
  color circlecenter, circlemoment;
  numeric ang, ind, i, anglim;
  path cpath, apath, opath, ipath, wp, ep;
  pair outerp[], innerp[], refpair;
  pen markpen;
  boolean cuspcond;
  picture holepic;
  color hvec, fakef, visualcircenter, hvecf, perpvec, axisvec;
  numeric rratio, hratio, visualray, resultang;
  numeric fakedist, fakeangle, fakeviewlimitangle, tmpangle;
  markpen = pencircle scaled 5pt;
  pickup markpen;
  bigray = 5;
  smaray = 1;
  %message "Villarceau angle: " & decimal(angle(bigray+-+smaray,smaray));
  message "Actual view angle: " & decimal(angle(X(f)++Y(f),Z(f)));
  angstep= 0.5;
  toruscenter = black;
  torusmoment = (0,0,1);
  % f := (3,-6,3);
  %%%%%%%%%%%%%%%%%%%%%%%%%%% f := 3*(0,4.24706,0.67267); fakedist = 7.3;
  f := (7,4,2); fakedist = 4.3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fakeangle = 9;
  fakeviewlimitangle = 2;
%  fakef = (0,4,3);
%  fakef = f;
%  fakef = (-Y(f),X(f),Z(f));
  fakef = (0,fakedist*cosd(fakeangle),fakedist*sind(fakeangle));
  %show fakef;
  %endfig;
  %end.
  if cdotprod(f-viewcentr,torusmoment)<0 :
    torusmoment := -torusmoment;
  fi;
  refpair = unitvector( rp(toruscenter+torusmoment)-rp(toruscenter) );
  viewline = f-toruscenter;
  ind = 0;
  sideaxe = bigray*ncrossprod( torusmoment, viewline );
  nearaxe = bigray*ncrossprod( sideaxe, torusmoment );
  anglim = 0.5*angstep-180;
  hvec = cdotprod(f-toruscenter,N(torusmoment))*N(torusmoment);
  hratio = conorm(hvec)/bigray;
  rratio = conorm(f-toruscenter-hvec)/bigray;
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
  resultang = anglim+angstep*(i-1);
  write decimal(smaray/bigray) & " " & decimal(rratio) & " " & decimal(hratio) & " " & decimal(resultang) to "torusexperimentoneline.dat";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ep = outerp[ind-i]--innerp[ind-i];
  wp = innerp[i]--outerp[i];
  unfill buildcycle(opath,ep,reverse ipath,wp);
  draw opath withcolor green;
  draw subpath (i,ind-i) of ipath withcolor red;
  hvecf = cdotprod(fakef-toruscenter,N(torusmoment))*N(torusmoment);
  visualcircenter = 0.5[toruscenter+hvecf,fakef];
  axisvec = fakef-visualcircenter;
  visualray = conorm( axisvec );
  perpvec = ncrossprod(axisvec,torusmoment)*visualray;
  color tmpf, tmpc, tmpp[], tmpq[], tmpa, tmpv;
  numeric theta, alpha, beta;
  tmpv = N(torusmoment);
  i := 0;
  drawoptions( withcolor blue withpen pencircle scaled 5pt );
  for theta=2*anglim step angstep until -2*anglim:
    tmpf := visualcircenter+planarrotation( axisvec, perpvec, theta );  
    tmpa := N( planarrotation( axisvec, perpvec, theta/2 ) );
    tmpc := toruscenter+tmpa*bigray;
    beta := getangle(tmpa,tmpf-tmpc);
    alpha := angle( smaray, conorm( tmpf-tmpc ) +-+ smaray ); %%%%%%%%%%%%
    tmpp[i] := tmpc + smaray*planarrotation( -tmpa, tmpv, 180-alpha-beta );
    tmpq[i] := tmpc + smaray*planarrotation( -tmpa, tmpv, 180+alpha-beta );
    if i>0:
      tmpangle := getangle(tmpp[i]-tmpp[i-1],tmpp[i]-fakef);
     %write decimal(i) & " " & decimal(tmpangle) to "torusexperimentangles.dat";
      if (tmpangle<fakeviewlimitangle) or (180-tmpangle<fakeviewlimitangle):
	draw rp(tmpp[i])--rp(fakef) withpen pencircle scaled 0pt;
      fi;
    fi;
    i := incr(i);
  endfor;
  draw rp(fakef);
  if smaray>bigray:
    draw rp(toruscenter+tmpv*(smaray+-+bigray));
    draw rp(toruscenter-tmpv*(smaray+-+bigray));
  fi;
  ind := i-1;
  drawoptions();
  pickup pencircle scaled 0pt;  
  draw rp(fakef)--rp((X(fakef),Y(fakef),0))--rp((0,bigray,0))--cycle dashed evenly;
  draw goodcirclepath((0,bigray,0),red,smaray) dashed evenly;
  draw for i=0 upto ind: rp(tmpp[i]).. endfor cycle;
  draw for i=0 upto ind: rp(tmpq[i]).. endfor cycle;
  pickup pencircle scaled 2pt;
  draw spatialhalfcircle( toruscenter, torusmoment, bigray+smaray, true );
  %draw rigorouscircle( toruscenter, torusmoment, bigray-smaray );
  draw rigorouscircle( toruscenter+smaray*tmpv, torusmoment, bigray );
  %draw rigorouscircle( toruscenter-smaray*tmpv, torusmoment, bigray );
  picture tmppicture;
  path zoomwindow;
  zoomwindow = (550,230)--(750,230)--(750,330)--(550,330)--cycle;
  tmppicture = currentpicture;
  clip tmppicture to zoomwindow;
endfig;

beginfig(2);
  draw tmppicture scaled 4;
endfig;

end.


