%%% file fig-clash-svg-dvips.mp
%%% Copyright 2009, 2011 Stephan Hennig
%
% This work may be distributed and/or modified under the conditions of
% the LaTeX Project Public License, either version 1.3 of this license
% or (at your option) any later version.  The latest version of this
% license is in http://www.latex-project.org/lppl.txt
%
outputtemplate := "%j-%c.mps";
input mpcolornames
dvipsnames;
picture disc; disc := image(fill fullcircle scaled 50);
beginfig(1);
  draw disc withcolor Lavender;
  begingroup
    svgnames;
    draw disc shifted (75,0) withcolor Lavender;
  endgroup;
  draw disc shifted (150,0) withcolor Lavender;
endfig;
end
