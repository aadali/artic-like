%%% file proof-mpcolornames.mp
%%% Copyright 2009,2011 Stephan Hennig
%
% This work may be distributed and/or modified under the conditions of
% the LaTeX Project Public License, either version 1.3 of this license
% or (at your option) any later version.  The latest version of this
% license is in http://www.latex-project.org/lppl.txt
%
input mpcolornames

boxh := 7.5bp;
boxw := 12.5bp;
penw := .3bp;
string specname;

def colorbox=
  unitsquare xscaled boxw yscaled boxh
enddef;

def proof(suffix col)=
  outputtemplate := "proof-spec-" & specname & "-" & str col & ".mps";
  beginfig(1);
    fill colorbox withcolor col;
    draw colorbox withpen pencircle scaled penw;
  endfig;
enddef;

endinput
