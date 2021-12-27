% steamcamera.mp
% L. Nobre G.
% 2006

input featpost3Dplus2D;

verbatimtex
\documentclass{article} 
\usepackage{palatino,mathpazo} 
\begin{document} 
etex

Spread := 50;
f := (3,5,2);

beginfig(0);
  numeric slabthick, camethick, blocthick, blochange;
  numeric slabside, cameside, blocside;
  color slabpos, camepos, blocpos;
  numeric markmargin;
  color lowonbloc, lowonmark, lowonextr;
  color higonbloc, higonmark, higonextr;
  color newlevelcenter;
  slabthick = 0.07;
  camethick = 0.7;
  blocthick = 0.6;
  slabside = 1.2;
  cameside = 0.9;
  blocside = 0.7;
  blochange = 0.25;
  markmargin = 0.04;
  slabpos = black;
  camepos = -0.5*(slabthick+camethick)*blue;
  blocpos = 0.5*slabthick*blue;
  lowonbloc = blocpos+(0,0.5*blocside,blocthick-blochange);
  higonbloc = blocpos+(0,0.5*blocside,blocthick);
  lowonmark = blocpos+(-0.5*slabside,0.5*blocside,blocthick-blochange);
  higonmark = blocpos+(-0.5*slabside,0.5*blocside,blocthick);
  lowonextr = lowonmark - markmargin*red;
  higonextr = higonmark - markmargin*red;
  kindofcube( false, false, camepos, 0, 0, 0, cameside, cameside, camethick );
  kindofcube( false, false, slabpos, 0, 0, 0, slabside, slabside, slabthick );
  rigorousdisc( 0, true, blocpos, 0.5*blocside, blocthick*blue );
  newlevelcenter = (blocthick-blochange+0.5*slabthick)*blue;
  draw
    spatialhalfcircle( newlevelcenter, blue, 0.5*blocside, true )
    dashed evenly;
  draw rp(lowonbloc)--rp(lowonextr);
  draw rp(higonbloc)--rp(higonextr);
  drawdblarrow rp(lowonmark)--rp(higonmark);
  label.rt( btex $\Delta y$ etex, rp(0.5[lowonmark,higonmark]) );
endfig;

verbatimtex 
\end{document} 
etex

end.