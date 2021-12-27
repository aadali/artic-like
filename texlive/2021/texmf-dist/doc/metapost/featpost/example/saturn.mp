% saturn.mp
% L. Nobre G.
% 2007

input featpost3Dplus2D;

verbatimtex
\documentclass{article}
\usepackage{newcent}
\pagestyle{empty}
\begin{document}
etex

f := 2.5*(5,0,-1);
Spread := 100;

beginfig(1);
  path planet, twilight, darkness, ringins, ringout, shadein, shadout;
  path fullplanet, fullring, backdark, sunplanet, backring;
  path nearringout, nearringins, hdtv, darkplanet;
  numeric planetray, ringinsray, ringoutray, cutmargin, athmosray;
  color sunposition, axedirection, labelpos, darknesstart, darknessend;
  numeric inclinangle, labelray, begintime, endtime, colorfactor;
  color planetcolor, ringcolor, toptextcolor, deepskycolor, bottextcolor;
  picture storepicture;
  planetray = 1;
  ringinsray = 1.6;
  ringoutray = 2.3;
  labelray = 2.1;
  cutmargin = 0.01;
  inclinangle = 15;
  sunposition = (-7,15,4);
  colorfactor = 1.3;
  bboxmargin := 10mm;
  planetcolor = colorfactor*(255,152,55)/255;
  ringcolor = colorfactor*(221,57,199)/255;
  toptextcolor = colorfactor*(187,164,109)/255;
  deepskycolor = colorfactor*(11,18,84)/255;
  bottextcolor = colorfactor*(15,103,23)/255;
  viewcentr := (0,0.5[planetray,ringinsray],0);
  labelpos = (0,labelray,0);
  athmosray = planetray+cutmargin;
  axedirection = (0,-sind(inclinangle),cosd(inclinangle));
  planet = rigorousfearpath( black, planetray );
  twilight = spatialhalfcircle( black, sunposition, athmosray, true );
  darknesstart = athmosray*ncrossprod( sunposition, axedirection );
  darknessend = darknesstart - labelray*N( sunposition );
  darkness = pathofstraightline( darknesstart, darknessend );
  darknesstart := -darknesstart;
  darknessend := darknesstart - labelray*N( sunposition );
  backdark = pathofstraightline( darknesstart, darknessend );
  ringins = goodcirclepath( black, axedirection, ringinsray );
  ringout = goodcirclepath( black, axedirection, ringoutray );
  endtime = xpart ( planet intersectiontimes twilight );
  begintime = xpart ( reverse planet intersectiontimes twilight );
  sunplanet = subpath (-begintime-cutmargin,endtime+cutmargin) of planet;
  darkplanet = planet cutbefore twilight cutafter twilight;
  fullplanet = buildcycle( sunplanet, reverse twilight );
  fill fullplanet withcolor planetcolor;
  endtime := xpart ( ringins intersectiontimes darkness );
  begintime := xpart ( reverse ringins intersectiontimes sunplanet );
  nearringins = subpath (-begintime-cutmargin,endtime+cutmargin) of ringins;
  endtime := xpart ( ringout intersectiontimes darkness );
  begintime := xpart ( reverse ringout intersectiontimes sunplanet );
  nearringout = subpath (-begintime-cutmargin,endtime+cutmargin) of ringout;
  backring = buildcycle( backdark, ringout, reverse darkplanet, reverse ringins );
  fullring = buildcycle( nearringout, reverse darkness, reverse nearringins, sunplanet );
  fill backring withcolor ringcolor;
  fill fullring withcolor ringcolor;
  picture thestatement;
  thestatement = thelabel(
    btex
      \begin{tabular}{c}
        \textbf{\textsf{{\Large T}$\!$ecnumbria}} \\
        \textsf{\scriptsize a sombra da tecnologia}
      \end{tabular}
    etex,
    origin);
  pickup pencircle scaled 1pt;
  draw thestatement scaled 2 shifted rp(labelpos) withcolor toptextcolor;
  storepicture = currentpicture;
  z1 = lrcorner currentpicture;
  z2 = urcorner currentpicture;
  z3 = ulcorner currentpicture;
  z4 = llcorner currentpicture;
  z5 = z1+down*bboxmargin;
  z6 = z1+right*bboxmargin;
  z7 = z2+right*bboxmargin;
  z8 = z2+up*bboxmargin;
  z9 = z3+up*bboxmargin;
  z10= z3+left*bboxmargin;
  z11= z4+left*bboxmargin;
  z12= z4+down*bboxmargin;
  z13= 0.5[z5,z12]+down*3*bboxmargin;
  hdtv = z5..z6---z7..z8---z9..z10---z11..z12---cycle;
  message "Aspect ratio: " & decimal( (xpart (z7-z10))/(ypart (z8-z5)) );
  numeric wid;
  wid = (xpart (z7-z10))/(1mm);
  message "Width: " & decimal( wid ) & "mm";
  message "Resolution for MooMiniCard: " & decimal(1.9083*222/wid) & "*100dpi";
  %% I don't know why there is need for the 1.9083 factor...
  fill hdtv withcolor deepskycolor;
  draw storepicture;
%  label(btex \texttt{\LARGE http://tecnumbria.com} etex, z13 );
endfig;

% beginfig(2);
%   prologues := 1;
%   pickup pencircle scaled 9pt;
%   fill hdtv;
%   unfill fullplanet;
%   draw fullplanet;
%   unfill fullring;
%   draw fullring;
%   unfill backring;
%   draw backring;
% endfig;

% beginfig(3);
%   prologues := 1;
%   pickup pencircle scaled 0pt;
%   draw fullplanet withcolor red;
%   draw fullring withcolor green;
%   draw backring withcolor blue;
% endfig;

verbatimtex \end{document} etex

end;

