% nurbstobezier.mp
% Troy Henderson
% 2007

prologues := 1;

% Evaluate a cubic polynomial of the "standard" Bezier form at t
vardef evalbezier(expr p,t) =
  save _a,_b,_c,_d;
  numeric _a,_b,_c,_d;
  _a:=(1-t)**3;
  _b:=3*((1-t)**2)*t;
  _c:=3*(1-t)*(t**2);
  _d:=t**3;
  (point 0 of p)*_a + (postcontrol 0 of p)*_b + (precontrol 1 of p)*_c + (point 1 of p)*_d
enddef;

% Evaluate the derivative of a cubic polynomial of the "standard"
% Bezier form at t
vardef evalbezierderivative(expr p,t) =
  save _a,_b,_c;
  pair _a,_b,_c;
  _a:=3*((point 1 of p) - 3*(precontrol 1 of p) + 3*(postcontrol 0 of p) -(point 0 of p));
  _b:=6*((precontrol 1 of p) - 2*(postcontrol 0 of p) + (point 0 of p));
  _c:=3*((postcontrol 0 of p) - (point 0 of p));
  _a*(t**2) + _b*t + _c
enddef;

% Evaluate a rational function of the "standard" cubic NURBS form at t
vardef evalnurbs(expr p,w,t) =
  save _q,_r;
  path _q,_r;
  _q:=((cyanpart w)*(point 0 of p)).. controls ((magentapart w)*(postcontrol 0 of p)) and ((yellowpart w)*(precontrol 1 of p)) .. ((blackpart w)*(point 1 of p));
  _r:=(cyanpart w,0) .. controls (magentapart w,0) and (yellowpart w,0) .. (blackpart w,0);
  evalbezier(_q,t)/(xpart evalbezier(_r,t))
enddef;

% Evaluate the derivative of a rational function of the "standard"
% cubic NURBS form at t
vardef evalnurbsderivative(expr p,w,t) =
  save _a,_b,_c,_d,_q,_r;
  pair _a,_b;
  numeric _c,_d;
  path _q,_r;
  _q:=((cyanpart w)*(point 0 of p)) .. controls ((magentapart w)*(postcontrol 0 of p)) and ((yellowpart w)*(precontrol 1 of p)) .. ((blackpart w)*(point 1 of p));
  _r:=(cyanpart w,0) .. controls (magentapart w,0) and (yellowpart w,0) .. (blackpart w,0);
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
  _xabs:=max(abs(xpart point 0 of p),abs(xpart postcontrol 0 of p),abs(xpart precontrol 1 of p),abs(xpart point 1 of p));
  _yabs:=max(abs(ypart point 0 of p),abs(ypart postcontrol 0 of p),abs(ypart precontrol 1 of p),abs(ypart point 1 of p));
  if (_xabs > 0):
    _a:=xpart((point 1 of p) - 3*(precontrol 1 of p) + 3*(postcontrol 0 of p) - (point 0 of p))/_xabs;
    _b:=3*xpart((precontrol 1 of p) - 2*(postcontrol 0 of p) + (point 0 of p))/_xabs;
    _c:=3*xpart((postcontrol 0 of p) - (point 0 of p))/_xabs;
    _d:=xpart(point 0 of p)/_xabs;
    _i:=(_a**2)/7 + ((_b)**2 + 2*_a*_c)/5 + (_a*_b + 2*_b*_d + (_c**2))/3 + (_a*_d + _b*_c)/2 + (_c*_d + (_d**2));
  else:
    _i:=0;
  fi;
  if (_yabs > 0):
    _A:=ypart((point 1 of p) - 3*(precontrol 1 of p) + 3*(postcontrol 0 of p) - (point 0 of p))/_yabs;
    _B:=3*ypart((precontrol 1 of p) - 2*(postcontrol 0 of p) + (point 0 of p))/_yabs;
    _C:=3*ypart((postcontrol 0 of p) - (point 0 of p))/_yabs;
    _D:=ypart(point 0 of p)/_yabs;
    _I:=(_A**2)/7 + ((_B)**2 + 2*_A*_C)/5 + (_A*_B + 2*_B*_D + (_C**2))/3 + (_A*_D + _B*_C)/2 + (_C*_D + (_D**2));
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
    _e:=((point 0 of _q[2])-(point 0 of _q[3])) .. controls ((postcontrol 0 of _q[2])-(postcontrol 0 of _q[3])) and ((precontrol 1 of _q[2])-(precontrol 1 of _q[3])) .. ((point 1 of _q[2])-(point 1 of _q[3]));
    _error:=beziernorm(_e)/beziernorm(_q[3]);
    show _error;
    if (_error > EPS):
      _b:=_c;
    else:
      _q[2]:=beziertobezier(_q[1],_c,_b);
      _q[3]:=nurbstobezier(p,w,_c,_b);
      _e:=((point 0 of _q[2])-(point 0 of _q[3])) .. controls ((postcontrol 0 of _q[2])-(postcontrol 0 of _q[3])) and ((precontrol 1 of _q[2])-(precontrol 1 of _q[3])) .. ((point 1 of _q[2])-(point 1 of _q[3]));
      _error:=beziernorm(_e)/beziernorm(_q[3]);
      if (_error > EPS):
	_q:=_q .. controls (postcontrol 0 of _q[4]) and (precontrol 1 of _q[4]) .. (point 1 of _q[4]);
	_a:=_c;
	_q[4]:=_q[3];
      else:
	_q:=_q .. controls (postcontrol 0 of _q[1]) and (precontrol 1 of _q[1]) .. (point 1 of _q[1]);
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
    _q:=_q .. ((_a*(point 0 of p)+_b*(postcontrol 0 of p)+_c*(precontrol 1 of p)+_d*(point 1 of p))/(_a+_b+_c+_d));
  endfor;
  ( _q )
enddef;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here's where the fun begins %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beginfig(0);
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
  Err:=0.040;

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

end
