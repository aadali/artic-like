% Admittedly silly but instructive example;
input fix_page;

def silly_display expr p =
 draw p shifted (.5,.25) withcolor red; draw p;
 pickup pensquare scaled .5pt;
 bboxmargin:=0;
 draw bbox p withcolor blue;
 pickup pencircle scaled .5pt;
 draw lrcorner p -- ulcorner p withcolor blue;
 draw urcorner p -- llcorner p withcolor blue;
enddef;

beginfig(101); silly_display fix_page(1, "storyx.mpy"); endfig;
beginfig(102); silly_display fix_page(2, "storyx.mpy"); endfig;
beginfig(103); silly_display fix_page(3, "storyx.mpy"); endfig;

bye.
