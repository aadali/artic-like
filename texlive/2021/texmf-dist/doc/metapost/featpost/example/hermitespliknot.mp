% hermitespliknot.mp
% L. Nobre G. 
% 2006

input featpost3Dplus2D;

% bash: export TEX=latex
% tcsh: setenv TEX latex

verbatimtex
\documentclass{article} 
\usepackage{palatino,mathpazo} 
\begin{document} 
etex

numeric param;
param = 0.7;

Spread := 39;
f := 1.1*(3,5,4);

vardef Pne( expr pr ) = ( (2*pr-3)*(pr**2)+1 ) enddef;
vardef Pwo( expr pr ) = ( (3-2*pr)*(pr**2)   ) enddef;
vardef Tne( expr pr ) = ( pr*(pr*(pr-2)+1)   ) enddef;
vardef Two( expr pr ) = ( (pr-1)*(pr**2)     ) enddef;

vardef Wei( expr pre, thi, nex ) =
  ( 0.5*( nex - pre ) )
enddef;

vardef Hermite( expr pr, zer, one, two, thr ) =
  ( Pne(pr)*one+Pwo(pr)*two+Tne(pr)*Wei(zer,one,two)+Two(pr)*Wei(one,two,thr) )
enddef;

beginfig(1);
  color ki[], actpoi, newpoi;
  numeric inum, bnum, del, halz, totz, pstep, p, bray;
  numeric icen, ap, istart, ifinis, thfac, vhfac;
  inum = 12;
  bnum = 500;  
  del = 8;
  halz = 1.5;
  totz = 3.0;
  pstep = 0.002;
  bray = 0.05;
  icen = 1;
  thfac = 3.5;
  vhfac = 0.8;
  ki0 = (-6,0,0);
  ki1 = (-5,0,0);
  ki2 = (-1.5-1.5*cosd((180-del)*param)-cosd((180+del)*param),
    -sind((180+del)*param)+1.5*sind((180-del)*param),halz);
  ki3 = (-1.5-1.5*cosd((180-del)*param),
    1.5*sind((180-del)*param),totz);
  ki4 = (-1.5-1.5*cosd((180-del)*param)+cosd((180+del)*param),
    sind((180+del)*param)+1.5*sind((180-del)*param),halz);
  ki5 = (cosd((180-del)*param),sind((180-del)*param),0);
  ki6 = (0,0,0);
  ki7 = (-cosd((180-del)*param),-sind((180-del)*param),0);
  ki8 = (1.5+1.5*cosd((180-del)*param)-cosd((180+del)*param),
    -sind((180+del)*param)+1.5*sind((180-del)*param),halz);
  ki9 = (1.5+1.5*cosd((180-del)*param),
    1.5*sind((180-del)*param),totz);
  ki10= (1.5+1.5*cosd((180-del)*param)+cosd((180+del)*param),
    sind((180+del)*param)+1.5*sind((180-del)*param),halz);
  ki11= (5,0,0);
  ki12= (6,0,0);

  istart = icen;
  ifinis = icen+1;
  ap := 0;
  forever:
    ap := ap+pstep;
    actpoi := Hermite(ap,ki[istart-1],ki[istart],ki[ifinis],ki[ifinis+1]);
    exitif conorm(actpoi-ki[istart])>bray;
  endfor;
  getready("signalvertex("
      & cstr(actpoi)
      & "," & decimal(thfac)
      & ",(0,0,0));signalvertex("
      & cstr(actpoi)
      & "," & decimal(thfac*vhfac)
      & ",(1,1,1));", actpoi);
  p = 1;
  forever:
    ap := ap+pstep;
    if ap>1:
      ap := ap-1;
      istart := ifinis;
      ifinis := istart+1;
    fi;
    exitif (ifinis>inum-1) or (p>bnum);
    newpoi := Hermite(ap,ki[istart-1],ki[istart],ki[ifinis],ki[ifinis+1]);
    if conorm(newpoi-actpoi)>2*bray:
      actpoi := newpoi;
      getready("signalvertex("
	  & cstr(actpoi)
	  & "," & decimal(thfac)
	  & ",(0,0,0));signalvertex("
	  & cstr(actpoi)
	  & "," & decimal(thfac*vhfac)
	  & ",(1,1,1));", actpoi);
      p := incr( p );
    fi;
  endfor;
  doitnow;
%  for p=0 upto inum:
%    dotlabel.top( decimal(p), rp(ki[p]) );
%  endfor;
%  produce_auto_scale;
endfig;

verbatimtex 
\end{document} 
etex

end;
