% splineperspective.mp
% L. Nobre G. and Troy Henderson
% 2007 -- 2012

prologues := 1;
defaultfont := "cmss17";
color f, viewcentr;
boolean	 ParallelProj;

            f := (3,5,4);   % This f is the point of view in 3D
	    
	    viewcentr := black; % This is the aim of the view

	    ParallelProj := false; % Kind of perspective %

    def X(expr A) = 
      if color A: redpart A else: cyanpart A fi
    enddef;
  
    def Y(expr A) =
      if color A: greenpart A else: magentapart A fi
    enddef;
  
    def Z(expr A) =  
      if color A: bluepart A else: yellowpart A fi
    enddef;

    def conorm(expr A) = 
      ( X(A) ++ Y(A) ++ Z(A) )  
    enddef;

    def N(expr A) =
        begingroup
            save M, exitcolor;
            numeric M;
            color exitcolor;
            M = conorm( A );
            if M > 0:
                exitcolor = ( X(A)/M, Y(A)/M, Z(A)/M );
            else:
                exitcolor := black;
            fi;
            ( exitcolor )
        endgroup
    enddef;

    def cdotprod(expr A, B) = 
        ( X(A)*X(B) + Y(A)*Y(B) + Z(A)*Z(B) )
    enddef;

    def ccrossprod(expr A, B) = 
        ( Y(A)*Z(B) - Z(A)*Y(B), 
          Z(A)*X(B) - X(A)*Z(B), 
          X(A)*Y(B) - Y(A)*X(B) )
    enddef;

% The dotproduct of two normalized vectors is the cosine of the angle 
% they form.

    def ndotprod(expr A, B) = 
        begingroup
            save a, b;
            color a, b;
            a = N(A);
            b = N(B);
            ( ( X(a)*X(b) + Y(a)*Y(b) + Z(a)*Z(b) ) )
        endgroup
    enddef;

% The normalized crossproduct of two vectors. 
% Also check getangle below.

    def ncrossprod(expr A, B) = 
        N( ccrossprod( A, B ) )
    enddef;

% Haahaa! Trigonometry. 

    def getangle(expr A, B) =
        begingroup 
            save coss, sine;
            numeric coss, sine;
            coss := cdotprod( A, B );
            sine := conorm( ccrossprod( A, B ) );
            ( angle( coss, sine ) )
        endgroup
    enddef;

    def rp(expr R) =
        begingroup

	  save v, u;
	  save verti, horiz, eta, squarf, radio;
	  color v, u;
	  numeric verti, horiz, eta, squarf, radio;
	  
	  v = N( (-Y(f-viewcentr), X(f-viewcentr), 0) );
	  u = ncrossprod( f-viewcentr, v );
	  
	  horiz = cdotprod( R-viewcentr, v );
	  verti = cdotprod( R-viewcentr, u );

	  if ParallelProj:
	    eta = 1;
	  else:
	    squarf = cdotprod( f-viewcentr, f-viewcentr );
	    radio = cdotprod( R-viewcentr, f-viewcentr );
	    eta = 1 - radio/squarf;
	  fi;
	  ( 150*(horiz,verti)/eta )

        endgroup
    enddef;

    def cartaxes(expr axex, axey, axez) =
        begingroup
            save orig, axxc, ayyc, azzc;
            color orig, axxc, ayyc, azzc;
            orig = (0,0,0);
            axxc = (axex,0,0);
            ayyc = (0,axey,0);
            azzc = (0,0,axez);
            drawarrow rp(orig)..rp(axxc);
            drawarrow rp(orig)..rp(ayyc);
            drawarrow rp(orig)..rp(azzc);
            label.bot( "x" ,rp(axxc));       %%%%%%%%%%%%%%%%%%%%%%%%% 
            label.bot( "y" ,rp(ayyc));       %%   Some  Labels...   %%
            label.lft( "z" ,rp(azzc));       %%%%%%%%%%%%%%%%%%%%%%%%% 
        endgroup
    enddef;

    def line( expr Ang ) =
      begingroup
	numeric a, b, c;
	a = (2-(1 ++ cosd(Ang))*cosd(3*Ang))*cosd(Ang);
	b = (2-(1 ++ cosd(Ang))*cosd(3*Ang))*sind(Ang);
	c =1.5+(1 ++ cosd(Ang))*sind(3*Ang);
	( (a,b,c) )
      endgroup
    enddef;
    
% Evaluate a cubic polynomial of the "standard" Bezier form at t
vardef evalbezier(expr p,t) =
  save _a,_b,_c,_d;
  numeric _a,_b,_c,_d;
  _a:=(1-t)**3;
  _b:=3*((1-t)**2)*t;
  _c:=3*(1-t)*(t**2);
  _d:=t**3;
  (point 0 of p)*_a + (postcontrol 0 of p)*_b + (precontrol 1 of p)*_c +
  (point 1 of p)*_d
enddef;

% Evaluate the derivative of a cubic polynomial of the "standard"
% Bezier form at t
vardef evalbezierderivative(expr p,t) =
  save _a,_b,_c;
  pair _a,_b,_c;
  _a:=3*((point 1 of p) - 3*(precontrol 1 of p) + 3*(postcontrol 0 of p)
    -(point 0 of p));
  _b:=6*((precontrol 1 of p) - 2*(postcontrol 0 of p) + (point 0 of p));
  _c:=3*((postcontrol 0 of p) - (point 0 of p));
  _a*(t**2) + _b*t + _c
enddef;

% Evaluate a rational function of the "standard" cubic NURBS form at t
vardef evalnurbs(expr p,w,t) =
  save _q,_r;
  path _q,_r;
  _q:=((cyanpart w)*(point 0 of p))..
  controls ((magentapart w)*(postcontrol 0 of p))
  and ((yellowpart w)*(precontrol 1 of p)) .. ((blackpart w)*(point 1 of p));
  _r:=(cyanpart w,0) ..
  controls (magentapart w,0) and (yellowpart w,0) .. (blackpart w,0);
  evalbezier(_q,t)/(xpart evalbezier(_r,t))
enddef;

% Evaluate the derivative of a rational function of the "standard"
% cubic NURBS form at t
vardef evalnurbsderivative(expr p,w,t) =
  save _a,_b,_c,_d,_q,_r;
  pair _a,_b;
  numeric _c,_d;
  path _q,_r;
  _q:=((cyanpart w)*(point 0 of p)) ..
  controls ((magentapart w)*(postcontrol 0 of p))
  and ((yellowpart w)*(precontrol 1 of p)) .. ((blackpart w)*(point 1 of p));
  _r:=(cyanpart w,0) ..
  controls (magentapart w,0) and (yellowpart w,0) .. (blackpart w,0);
  _a:=evalbezier(_q,t);
  _b:=evalbezierderivative(_q,t);
  _c:=xpart evalbezier(_r,t);
  _d:=xpart evalbezierderivative(_r,t);
  (_b*_c-_a*_d)/(_c**2)
enddef;

% Fit a cubic polynomial of the "standard" Bezier form to a
% rational function of the 
% "standard" cubic NURBS form with function and derivative agreement
% at tmin and tmax
vardef nurbstobezier(expr p,w,tmin,tmax) =
  save _a,_b,_c,_d,_e;
  pair _a,_b,_c,_d;
  numeric _e;
  _e:=(tmax-tmin)/3;
  _a:=evalnurbs(p,w,tmin);
  _b:=_a + _e*evalnurbsderivative(p,w,tmin);
  _d:=evalnurbs(p,w,tmax);
  _c:=_d - _e*evalnurbsderivative(p,w,tmax);
  _a .. controls _b and _c .. _d
enddef;

% Reparameterize a cubic polynomial of the "standard" Bezier form by mapping
% the interval [tmin,tmax] to [0,1]
vardef beziertobezier(expr p,tmin,tmax) =
  nurbstobezier(p,(1,1,1,1),tmin,tmax)
enddef;

% Evalute the L^2[0,1] norm of a cubic polynomial of the "standard"
% Bezier form
vardef beziernorm(expr p) =
  save _a,_b,_c,_d,_i,_xabs,_yabs,_A,_B,_C,_D,_I;
  numeric _a,_b,_c,_d,_i,_xabs,_yabs,_A,_B,_C,_D,_I;
  _xabs:=max(
    abs(xpart point 0 of p),
    abs(xpart postcontrol 0 of p),
    abs(xpart precontrol 1 of p),
    abs(xpart point 1 of p));
  _yabs:=max(
    abs(ypart point 0 of p),
    abs(ypart postcontrol 0 of p),
    abs(ypart precontrol 1 of p),
    abs(ypart point 1 of p));
  if (_xabs > 0):
    _a:=xpart((point 1 of p) - 3*(precontrol 1 of p)
      + 3*(postcontrol 0 of p) - (point 0 of p))/_xabs;
    _b:=3*xpart((precontrol 1 of p) - 2*(postcontrol 0 of p)
      + (point 0 of p))/_xabs;
    _c:=3*xpart((postcontrol 0 of p) - (point 0 of p))/_xabs;
    _d:=xpart(point 0 of p)/_xabs;
    _i:=(_a**2)/7 + ((_b)**2 + 2*_a*_c)/5 + (_a*_b + 2*_b*_d + (_c**2))/3 + (_a*_d + _b*_c)/2 + (_c*_d + (_d**2));
  else:
    _i:=0;
  fi;
  if (_yabs > 0):
    _A:=ypart((point 1 of p) - 3*(precontrol 1 of p)
      + 3*(postcontrol 0 of p) - (point 0 of p))/_yabs;
    _B:=3*ypart((precontrol 1 of p) - 2*(postcontrol 0 of p)
      + (point 0 of p))/_yabs;
    _C:=3*ypart((postcontrol 0 of p) - (point 0 of p))/_yabs;
    _D:=ypart(point 0 of p)/_yabs;
    _I:=(_A**2)/7 + ((_B)**2 + 2*_A*_C)/5
    + (_A*_B + 2*_B*_D + (_C**2))/3 + (_A*_D + _B*_C)/2 + (_C*_D + (_D**2));
  else:
    _I:=0;
  fi;
  (_xabs*sqrt(_i)) ++ (_yabs*sqrt(_I))
enddef;

% Fit a cubic Bezier spline to a rational function of the "standard"
% cubic NURBS form by iteratively refining the Bezier curve.
% p is a 4 point path containing the 4 cubic NURBS (2D) control points
% w is a cmykcolor containing the 4 cubic NURBS weights
% EPS is the tolerance to stop refining each branch of the Bezier spline
vardef fitnurbswithbezier(expr p,w,EPS) =
  save _a,_b,_c,_e,_error,_k,_q;
  numeric _a,_b,_c,_error,_k;
  path _q,_q[],_e;
  _a:=0;
  _b:=1;
  _k:=1/sqrt(2);
  _q:=(point 0 of p);
  _q[4]:=nurbstobezier(p,w,_a,_b);
  forever:
    exitunless(_a<1);
    _q[1]:=_q[4];
    _c:=_b-_k*((_b-_a)**2);
    _q[2]:=beziertobezier(_q[1],_a,_c);
    _q[3]:=nurbstobezier(p,w,_a,_c);
    _q[4]:=_q[3];
    _e:=((point 0 of _q[2])-(point 0 of _q[3])) ..
    controls ((postcontrol 0 of _q[2])-(postcontrol 0 of _q[3]))
    and ((precontrol 1 of _q[2])-(precontrol 1 of _q[3])) ..
    ((point 1 of _q[2])-(point 1 of _q[3]));
    _error:=beziernorm(_e)/beziernorm(_q[3]);
%    show _error;
    if (_error > EPS):
      _b:=_c;
    else:
      _q[2]:=beziertobezier(_q[1],_c,_b);
      _q[3]:=nurbstobezier(p,w,_c,_b);
      _e:=((point 0 of _q[2])-(point 0 of _q[3])) ..
      controls ((postcontrol 0 of _q[2])-(postcontrol 0 of _q[3]))
      and ((precontrol 1 of _q[2])-(precontrol 1 of _q[3])) ..
      ((point 1 of _q[2])-(point 1 of _q[3]));
      _error:=beziernorm(_e)/beziernorm(_q[3]);
      if (_error > EPS):
	_q:=_q .. controls (postcontrol 0 of _q[4])
	and (precontrol 1 of _q[4]) .. (point 1 of _q[4]);
	_a:=_c;
	_q[4]:=_q[3];
      else:
	_q:=_q .. controls (postcontrol 0 of _q[1])
	and (precontrol 1 of _q[1]) .. (point 1 of _q[1]);
	_a:=_b;
	_q[4]:=nurbstobezier(p,w,_a,1);
      fi;
      _b:=1;
    fi;
  endfor;
  _q
enddef;

% This macro is used to provide a path to draw the NURBS
% It returns a path of length N passing through N+1 equally spaced
% (in time) points along the NURBS connected by line segments
vardef samplednurbs(expr p,w,N) =
  save _a,_b,_c,_d,_n,_t,_q;
  numeric _a,_b,_c,_d,_n,_t;
  path _q;
  _q:=(point 0 of p);
  for _n=1 upto N:
    _t:=_n/N;
    _a:=(cyanpart w)*((1-_t)**3);
    _b:=3*(magentapart w)*((1-_t)**2)*_t;
    _c:=3*(yellowpart w)*(1-_t)*(_t**2);
    _d:=(blackpart w)*(_t**3);
    _q:=_q .. ((_a*(point 0 of p)+_b*(postcontrol 0 of p)
	+_c*(precontrol 1 of p)+_d*(point 1 of p))/(_a+_b+_c+_d));
  endfor;
  ( _q )
enddef;

% The code below is a development by Troy from an original by Przemek Koprowski

vardef rationalbezier(expr A,B,C,D) =
    begingroup
        save P,Q,E,a,b,c,d,r,EPS;
        color P[],Q[],E;
        pair a,b,c,d;
        EPS:=1/10;
        a:=(redpart A,greenpart A)/(bluepart A);
        b:=(redpart B,greenpart B)/(bluepart B);
        c:=(redpart C,greenpart C)/(bluepart C);
        d:=(redpart D,greenpart D)/(bluepart D);
        r:=max(abs(a-b),abs(a-c),abs(a-d),abs(b-c),abs(b-d),abs(c-d));
        if (max(abs(b-1/2[a,c]),abs(c-1/2[b,d])) > EPS*r):
            P[0]:=A;
            P[1]:=1/2[A,B];
               E:=1/2[B,C];
            Q[2]:=1/2[C,D];
            Q[3]:=D;
            P[2]:=1/2[P[1],E];
            Q[1]:=1/2[E,Q[2]];
            P[3]:=1/2[P[2],Q[1]];
            Q[0]:=P[3];
            rationalbezier(P[0],P[1],P[2],P[3]) & rationalbezier(Q[0],Q[1],Q[2],Q[3])
        else:
            a .. controls b and c .. d
        fi
    endgroup
enddef;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here's where the fun begins %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def casteljau( expr Za, Zb, Zc, Zd, Pt ) = %%%%%%%%%%%%%%%%%%% 2D or 3D
  begingroup
    save A, B, C, D;
    numeric A, B, C, D;
    A = (1-Pt)**3;
    B = 3*((1-Pt)**2)*Pt;
    C = 3*(1-Pt)*(Pt**2);
    D = Pt**3;
    ( (A*Za+B*Zb+C*Zc+D*Zd) )
  endgroup
enddef;

def twothr( expr Z ) = ( xpart Z, ypart Z, 0 ) enddef;

def twotwo( expr Z ) = rp( twothr( Z ) ) enddef;

def xoy( expr Z ) = rp( ( X(Z), Y(Z), 0 ) ) enddef;

def yoz( expr W ) = rp( ( 0, Y(W), Z(W) ) ) enddef;

def xoz( expr W ) = rp( ( X(W), 0, Z(W) ) ) enddef;

def nextthirty( expr Za, Zb, Zc, Zd, Pt ) = %%% input 3D and return 2D
  begingroup
    save A, B, C, D, Tot, P;
    numeric A, B, C, D, Tot;
    color P;
    P = N( f - viewcentr );
    A = ((1-Pt)**3)*cdotprod( P, f-Za );
    B = 3*((1-Pt)**2)*Pt*cdotprod( P, f-Zb );
    C = 3*(1-Pt)*(Pt**2)*cdotprod( P, f-Zc );
    D = (Pt**3)*cdotprod( P, f-Zd );
    Tot = A+B+C+D;
    ( (A*rp(Za)+B*rp(Zb)+C*rp(Zc)+D*rp(Zd))/Tot )
  endgroup
enddef;

vardef nurbstobezierold (expr p,w) =
  save _a,_b,_c,_d,_j,_n,_r,_s,_t,_A,_B,_Aold,_Bold,_C,_D,_EPS,_J,_N;
  _EPS:=0.00001;
  _J:=10;
  
  _Aold:=0;
  _Bold:=0;
  _A:=1;
  _B:=1;
  _s:=((_A-_Aold)++(_B-_Bold))/(_A++_B);
  _j:=1; _r:=0;
  forever:
    exitunless((_s>_EPS) and (_j<_J));
    _j:=_j+1;
    _N:=2**_j;
    _Aold:=_A;
    _Bold:=_B;
    _D:=_N+1/_N-21/_N/_N/_N-1/_N/_N/_N/_N/_N+20/_N/_N/_N/_N/_N/_N/_N;
    _C:=120*(2+2/_N/_N-5/_N/_N/_N/_N)/_D;
    _D:=60*(3+3/_N/_N+10/_N/_N/_N/_N)/_D;
    _c:=5/_N/_N/_N/_N;
    _a:=2+2/_N/_N-_c;
    _b:=2-3/_N/_N+_c;
    _c:=1+6/_N/_N+_c;
    _A:=(-2*(cyanpart p)*_a+(blackpart p)*_b)/_c;
    _B:=((cyanpart p)*_b-2*(blackpart p)*_a)/_c;
    for _n=0 upto _N:
      _t:=_n/_N;
      _a:=(1-_t)**3;
      _b:=((1-_t)**2)*_t;
      _c:=(1-_t)*(_t**2);
      _d:=_t**3;
      _r:=((cyanpart w)*(cyanpart p)*_a + 3*(magentapart
	  w)*(magentapart p)*_b + 3*(yellowpart w)*(yellowpart p)*_c +
	(blackpart w)*(blackpart p)*_d)/((cyanpart w)*_a + 3*(magentapart
	  w)*_b + 3*(yellowpart w)*_c + (blackpart w)*_d);
      _A:=_A+(_C*_b-_D*_c)*_r;
      _B:=_B+(_C*_c-_D*_b)*_r;
    endfor;
    _s:=((_A-_Aold)++(_B-_Bold))/(_A++_B);
  endfor;
  (_A,_B)/3
enddef;

def nurbsapprox( expr Pa, Pb, Pc, Pd ) =
  begingroup
    color Pn;
    numeric wa, wb, wc, wd;
    path returnpath;
    pair xpair, ypair, ba, bb, bc, bd;
    cmykcolor xcontrols, ycontrols;
    Pn = N( f - viewcentr );
    wa = cdotprod( Pn, f-Pa );
    wb = cdotprod( Pn, f-Pb );
    wc = cdotprod( Pn, f-Pc );
    wd = cdotprod( Pn, f-Pd );
    xcontrols = (xpart rp(Pa),xpart rp(Pb),xpart rp(Pc),xpart rp(Pd));
    ycontrols = (ypart rp(Pa),ypart rp(Pb),ypart rp(Pc),ypart rp(Pd));
    xpair = nurbstobezierold( xcontrols, (wa,wb,wc,wd) );
    ypair = nurbstobezierold( ycontrols, (wa,wb,wc,wd) );
    ba = rp( Pa );
    bb = (xpart xpair, xpart ypair);
    bc = (ypart xpair, ypart ypair);
    bd = rp( Pd );
    %show ba;
    %show bb;
    %show bc;
    %show bd;
    %show wa;
    %show wb;
    %show wc;
    %show wd;
    returnpath = ba..controls bb and bc..bd;
    (returnpath)
  endgroup
enddef;

def fitthreednurbswithtwodbezier( expr pa, pb, pc, pd, EPS ) =
  begingroup
    save _a,_b,_c,_e,_error,_k,_q,w,wa,wb,wc,wd,pn,za, zb, zc, zd;
    numeric _a,_b,_c,_error,_k,wa,wb,wc,wd;
    path p,_q,_q[],_e;
    color pn;
    pair za, zb, zc, zd;
    cmykcolor w;
    za = rp(pa); show za;
    zb = rp(pb); show zb;
    zc = rp(pc); show zc;
    zd = rp(pd); show zd;
    p = za .. controls zb and zc .. zd;
    pn = N( f - viewcentr );
    wa = cdotprod( pn, f-pa ); show wa;
    wb = cdotprod( pn, f-pb ); show wb;
    wc = cdotprod( pn, f-pc ); show wc;
    wd = cdotprod( pn, f-pd ); show wd;
    w = ( wa, wb, wc, wd );
    _a:=0;
    _b:=1;
    _k:=1/sqrt(2);
    _q:=(point 0 of p);
    _q[4]:=nurbstobezier(p,w,_a,_b);
    forever:
      exitunless(_a<1);
      _q[1]:=_q[4];
      _c:=_b-_k*((_b-_a)**2);
      _q[2]:=beziertobezier(_q[1],_a,_c);
      _q[3]:=nurbstobezier(p,w,_a,_c);
      _q[4]:=_q[3];
      _e:=((point 0 of _q[2])-(point 0 of _q[3])) ..
      controls ((postcontrol 0 of _q[2])-(postcontrol 0 of _q[3]))
      and ((precontrol 1 of _q[2])-(precontrol 1 of _q[3])) ..
      ((point 1 of _q[2])-(point 1 of _q[3]));
      _error:=beziernorm(_e)/beziernorm(_q[3]);
      if (_error > EPS):
	_b:=_c;
      else:
	_q[2]:=beziertobezier(_q[1],_c,_b);
	_q[3]:=nurbstobezier(p,w,_c,_b);
	_e:=((point 0 of _q[2])-(point 0 of _q[3])) ..
	controls ((postcontrol 0 of _q[2])-(postcontrol 0 of _q[3]))
	and ((precontrol 1 of _q[2])-(precontrol 1 of _q[3])) ..
	((point 1 of _q[2])-(point 1 of _q[3]));
	_error:=beziernorm(_e)/beziernorm(_q[3]);
	if (_error > EPS):
	  _q:=_q .. controls (postcontrol 0 of _q[4])
	  and (precontrol 1 of _q[4]) .. (point 1 of _q[4]);
	  _a:=_c;
	  _q[4]:=_q[3];
	else:
	  _q:=_q .. controls (postcontrol 0 of _q[1])
	  and (precontrol 1 of _q[1]) .. (point 1 of _q[1]);
	  _a:=_b;
	  _q[4]:=nurbstobezier(p,w,_a,1);
	fi;
	_b:=1;
      fi;
    endfor;
    ( _q )
  endgroup
enddef;

beginfig(0);
  draw for i=1 step 6 until 360: xoy(line(i)).. endfor cycle;
  draw for i=1 step 6 until 360: rp(line(i)).. endfor cycle;
  for i=1 step 15 until 360:
    draw rp(line(i)) withpen pencircle scaled 2mm;
  endfor;  
endfig;

beginfig(4);
% p contains the 4 control points of the rational function of the
% "standard" cubic NURBS form
  path p;
  p:=(297.63725,297.63725) .. controls (132.98871,286.67885) and (180.62535,152.16249) .. (429.54399,226.31157);

% w contains the 4 weights for the rational function of the
% "standard" cubic NURBS form
  cmykcolor w;
  w:=(2.15756,1.6709,0.8598,1.34647);
  
% EPS represents the minimum "acceptable error" to stop refining any
% given branch of the Bezier
  Err:=0.067;

% q represents the Bezier spline fit to the rational function of the
% "standard" cubic NURBS form
  path q;
  q:=fitnurbswithbezier(p,w,Err);
%  q:=fitnurbswithbezier(reverse p,(blackpart w,yellowpart w,magentapart w,cyanpart w),Err);

% draw the NURBS by sampling it at many points and connecting the
% samples via line segments
  draw samplednurbs(p,w,20) withcolor red withpen pencircle scaled 2bp;

% draw the Bezier spline and its knots
  draw q;
  for n=0 upto length q:
    draw fullcircle scaled 2 shifted point n of q withcolor blue;
  endfor;
endfig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f := 0.35*(3,5,2);           
beginfig(1);
  numeric tu, num, i, fac;
  pen pencontrol, penspline, penalytic;
  color colcontrol, colspline, colorytic, colormark;
  color w[];
  pencontrol = pencircle scaled 4pt;
  penspline = pencircle scaled 2pt;
  penalytic = pencircle scaled 1pt;
  colcontrol = black;
  colspline = red;
  colorytic = blue+green;
  colormark = (0.8,0.8,0.1);
  tu = 6cm;
  num = 50;
  fac = 1.2;
  transform T;
  T = identity scaled tu;
  z21 = origin; 
  z22 = (1,0);
  z23 = (1,1);
  z24 = (0,1);
  z1 = z21 transformed T; 
  z2 = z22 transformed T;
  z3 = z23 transformed T;
  z4 = z24 transformed T;
  z6 = (fac,0) transformed T;
  z8 = (0,fac) transformed T;
  drawarrow z1--z6;
  drawarrow z1--z8;
  label.lrt( "x", z6 );
  label.ulft( "y", z8 );
  dotlabels.urt(1,2,3,4);
  z11 = twotwo( z21 );
  z12 = twotwo( z22 );
  z13 = twotwo( z23 );
  z14 = twotwo( z24 );
  w1 = twothr( z21 );
  w2 = twothr( z22 );
  w3 = twothr( z23 );
  w4 = twothr( z24 );  
  cartaxes( fac, fac, 0.3*fac );
  draw z12--z13--z14 dashed evenly;
  draw z2--z3--z4 dashed evenly;
  % 1) Next line: MetaPost intrinsic path.
  draw z1..controls z2 and z3..z4 withpen penspline withcolor colspline;
  % 2) Next line: my implementation of the MetaPost intrinsic path.
  draw z1 for i=1 upto num: ..casteljau(z1,z2,z3,z4,i/num) endfor
                                  withpen penalytic withcolor colorytic;
  % 5) Next line: hopefully, how it should be done. Yeah! Way to go!
  draw z11 for i=1 upto num: ..nextthirty(w1,w2,w3,w4,i/num) endfor
                                  withpen pencontrol withcolor colcontrol;
  % 4) Next line: MetaPost intrinsic path of perspectived control points.
  draw z11..controls z12 and z13..z14
				  withpen penspline withcolor colspline;
  % 3) Next line: what should be drawn in perspective.
  draw z11 for i=1 upto num: ..rp(casteljau(w1,w2,w3,w4,i/num)) endfor
                                  withpen penalytic withcolor colorytic;
  % 6) Next line: Troy's approximation
  draw nurbsapprox(w1,w2,w3,w4) withcolor colormark;
  
				      
endfig;

f := 1.05*(3,5,2);           

beginfig(2);
  color w[];
  w1 = (1,0,0);
  w2 = (0,0,1);
  w3 = (0,1,0);
  w4 = (1,1,1);
  w5 = (1,1,0);
  w6 = (1,0,1);
  w7 = (0,1,1);
  cartaxes( fac, fac, fac );
  draw rp(w1)--rp(w2)--rp(w3)--rp(w4)--rp(w5)--rp(w1)--rp(w6)--
       rp(w2)--rp(w7)--rp(w3)--rp(w5) dashed withdots;
  draw rp(w6)--rp(w4)--rp(w7) dashed withdots;
  draw xoy(w1) for i=1 upto num:
    ..xoy(casteljau(w1,w3,black,w5,i/num))
  endfor withcolor colormark;
  draw xoy(w1) for i=1 upto num: ..xoy(casteljau(w1,w2,w3,w4,i/num)) endfor;
  draw xoz(w1) for i=1 upto num: ..xoz(casteljau(w1,w2,w3,w4,i/num)) endfor;
  draw rp(w1) for i=1 upto num: ..nextthirty(w1,w2,w3,w4,i/num) endfor
                                  withpen pencontrol withcolor colcontrol;
  draw rp(w1) for i=1 upto num: ..rp(casteljau(w1,w2,w3,w4,i/num)) endfor
                                  withpen penalytic withcolor colorytic;
  draw nextthirty(w1,w2,w3,w4,0.5) withpen pencontrol withcolor red;	      
endfig;

beginfig(5);
  color w[];
  for i=1 upto 4:
    w[i]=(uniformdeviate(1),uniformdeviate(1),uniformdeviate(1));
    draw rp(w[i]) withpen pencontrol;
    draw xoy(w[i]) withpen penspline;
    draw xoz(w[i]) withpen penspline;
    draw yoz(w[i]) withpen penspline;
  endfor;
  cartaxes( fac, fac, fac );
  draw rp(w1)--rp(w2)--rp(w3)--rp(w4) withpen penspline dashed evenly;
  draw xoy(w1)--xoy(w2)--xoy(w3)--xoy(w4) withpen penalytic dashed evenly;
  draw xoz(w1)--xoz(w2)--xoz(w3)--xoz(w4) withpen penalytic dashed evenly;
  draw yoz(w1)--yoz(w2)--yoz(w3)--yoz(w4) withpen penalytic dashed evenly;
  draw xoy(w1) for i=1 upto num: ..xoy(casteljau(w1,w2,w3,w4,i/num)) endfor;
  draw xoz(w1) for i=1 upto num: ..xoz(casteljau(w1,w2,w3,w4,i/num)) endfor;
  draw yoz(w1) for i=1 upto num: ..yoz(casteljau(w1,w2,w3,w4,i/num)) endfor;
  draw rp(w1) for i=1 upto num: ..nextthirty(w1,w2,w3,w4,i/num) endfor
                                  withpen pencontrol withcolor colcontrol;
%%%%%%%%%%%%% the following line may crash;
%%  draw fitthreednurbswithtwodbezier(w1,w2,w3,w4,0.2)
%%				  withpen penalytic withcolor red;	      
endfig;

beginfig(3);
  path xyp, xzp;
  xyp = origin..tension 2 and 0.75..right...(right+up+right+up);
  xzp = origin..tension 2 and 0.75..(right+up)...(right+up+right);
  z1 = postcontrol 0 of xyp;
  z2 = postcontrol 0 of xzp;
  z3 = precontrol 1 of xyp;
  z4 = precontrol 1 of xzp;
%  show (xpart z1);
%  show (xpart z2); 
%  show (xpart z3);
%  show (xpart z4);
  draw xyp;
  draw xzp;
endfig;

end.
