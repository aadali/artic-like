% Copyright 2016, Palle Jorgensen, hamselv(at)pallej.dk
%
% This file is part of the mparrows MetaPost module.
% 
% mparrows is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% Foobar is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with Foobar. If not, see <http://www.gnu.org/licenses/>.
verbatimtex
%&latex
\documentclass[11pt]{scrbook}
\usepackage{tgpagella}
\usepackage[T1]{fontenc}

\begin{document}
etex;

input mparrows;

filenametemplate "%j%c.mps";

for s within btex text etex:
  defaultfont:=fontpart.s;
  defaultscale:=xxpart.s;
  exitif true;
endfor

for i = 0 upto 7:
  beginfig(i);
    setarrows(i);
    drawarrow (origin{right}..{up}(0,2mm));
    drawdblarrow ((0,-2mm){up}..origin{right}..{up}(0,2mm)) shifted (5mm,0);
    currentpicture:=currentpicture scaled 5;
  endfig;
endfor

setarrows(barbed);
for i = 0 upto 10:
  beginfig(i+10);
    barbedarrowindent:=2*i/10;
    drawarrow origin{dir 20}..{dir 20}(1cm,0);
    currentpicture:=currentpicture scaled 5;
    label.rt("barbedarrowindent=" & decimal.barbedarrowindent, (55mm,0mm));
  endfig;
endfor;

beginfig(30);
  barbedarrowindent:=-.6;
  drawarrow origin{dir 20}..{dir 20}(1cm,0);
  currentpicture:=currentpicture scaled 5;
  label.rt("barbedarrowindent=" & decimal.barbedarrowindent, (55mm,0mm));
endfig;

setarrows(orig);

beginfig(31);
  drawarrow (origin{right}..{up}(0,2mm));
  draw _apth withcolor red;
  currentpicture:=currentpicture scaled 5;
endfig;

setarrows(default);

beginfig(32);
  drawarrow (origin{right}..{up}(0,2mm));
  draw _apth withcolor red;
  currentpicture:=currentpicture scaled 5;
endfig;

setarrows(barbedunfilled);

beginfig(33);
  barbedarrowindent:=.6;
  drawarrow (origin{right}..{up}(0,2mm));
  draw (makepath pencircle scaled 2ahlength) shifted (0,2mm) withcolor red withpen currentpen scaled .2;
  draw (makepath pencircle scaled (2*barbedarrowindent*ahlength)) shifted (0,2mm) withcolor green withpen currentpen scaled .2;
  currentpicture:=currentpicture scaled 5;

endfig;


end.
