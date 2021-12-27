% iamsorrykarl.mp
% 2006
% L. Nobre G.

%input featpost3Dplus2D;

% bash: export TEX=latex
% tcsh: setenv TEX latex

prologues := 1;

beginfig(1);
  numeric u, nstripes, aspecratio, minskip, maxskip, shiftlim, shiftmin;
  numeric nbreaks, i, j, ymax, shiftx, lefh, righ, actualx, newx;
  numeric prevlefh, prevrigh, xmax, auxan, auxr, cdiam;
  path lefp[], rigp[], strp[];
  pair shif, reflef, refrig, auxa, auxb;
  boolean exitcondition;
  u = 1cm;
  nstripes = 7;
  aspecratio = 1.7;
  minskip = 0.5u; % must not be too small
  maxskip = 4u;
  shiftlim = 0.7u;
  shiftmin = 0.25u;
  cdiam = 0.65u;
  xmax = (2*nstripes-1)*u;
  ymax = aspecratio*xmax;
  prevlefh = ymax;
  prevrigh = ymax;
  actualx = u;
  for i=1 upto nstripes:
    lefp[i] = ((2*i-1)*u,ymax);
    rigp[i] = ((2*i-1)*u,ymax)--(2*i*u,ymax);
  endfor;
  exitcondition = false;
  forever:
    auxr := uniformdeviate(2)-1;
    if auxr > 0:
      newx := actualx + shiftmin + (shiftlim-shiftmin)*auxr;
    else:
      newx := actualx - shiftmin + (shiftlim-shiftmin)*auxr;
    fi;
    lefh := prevlefh - minskip - uniformdeviate(maxskip-minskip);
    righ := prevrigh - minskip - uniformdeviate(maxskip-minskip);
    if lefh < minskip:
      lefh := minskip;
      exitcondition := true;
    fi;
    if righ < minskip:
      righ := minskip;
      exitcondition := true;
    fi;    
    reflef := (newx,lefh);
    refrig := (newx+xmax,righ);
    auxan := angle(refrig-reflef);
    shif := (newx-actualx)*(1,sind(auxan)/cosd(auxan));
    for i=1 upto nstripes:
      auxa := point (length lefp[i]) of lefp[i];
      auxb := point (length rigp[i]) of rigp[i];
      pair poil, poir;
      poil = auxa+whatever*down;
      poil = whatever[reflef,refrig];
      poir = auxb+whatever*down;
      poir = whatever[reflef,refrig];
      lefp[i] := lefp[i]--poil--(poil+shif);
      rigp[i] := rigp[i]--poir--(poir+shif);
    endfor;
    actualx := newx;
    prevlefh := lefh;
    prevrigh := righ;
    exitif exitcondition;
  endfor;
  for i=1 upto nstripes:
    auxa := point (length lefp[i]) of lefp[i];
    auxb := point (length rigp[i]) of rigp[i];
    pair poil, poir;
    poil = auxa+whatever*down;
    poil = whatever[origin,right];
    poir = auxb+whatever*down;
    poir = whatever[origin,right];
    lefp[i] := lefp[i]--poil--poir;
    rigp[i] := rigp[i]--poir;
    strp[i] = buildcycle( lefp[i], reverse rigp[i] );
    fill strp[i] withcolor red;
  endfor;
%  for i=nstripes+1 upto nstripes+50:
%    strp[i] = fullcircle scaled cdiam shifted
%              (uniformdeviate(xmax),uniformdeviate(ymax));
%    fill strp[i];
%  endfor;		   
endfig;

end.