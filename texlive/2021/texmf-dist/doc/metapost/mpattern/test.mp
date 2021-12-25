input mpattern;
prologues:=1;

beginpattern(szachy_i);
   patterncolor(.8);
   fill unitsquare scaled 4mm;
   fill unitsquare scaled 4mm shifted (4mm,4mm);
endpattern;

beginpattern(szachy_ii);
   fill unitsquare scaled 3mm withcolor .7white;
   fill unitsquare scaled 3mm shifted (3mm,3mm) withcolor .7white;
endpattern;

beginfig(1);
   fill fullcircle scaled 30mm withpattern szachy_i;

   picture q;
   q:=currentpicture; clearit;

   draw q shifted (10mm*right) rotated 30 scaled .5;

   fill unitsquare scaled 20mm withpattern szachy_ii;
endfig;

beginpattern(checker);
   u:=1mm;
   fill unitsquare scaled u;
   fill unitsquare scaled u shifted (u,u);
   patternbbox(0,0,2u,2u);
   patterntransform(identity slanted .17);
endpattern;

beginpattern(rotated_checker);
   fill unitsquare scaled 4mm rotated 45 withcolor .7white;
   patterntransform(identity rotated 30);
endpattern;

beginfig(11);
   path p;
   z1=(1cm,1cm);
   p=fullcircle scaled 3cm;
   fill p withpattern checker;
   unfill p shifted z1;
   fill p shifted z1 withpattern rotated_checker;
endfig;

beginpattern(tt);
   fill unitsquare scaled 2mm rotated 45 withcolor .5white;
   z1=llcorner currentpicture;
   z2=lrcorner currentpicture;
   z3=urcorner currentpicture;
   z4=ulcorner currentpicture;
%   draw z1--z2--z3--z4--cycle;
   patterntransform(identity rotated 30 slanted .2 yscaled 3);
   patternstep((1.1*xpart(urcorner currentpicture-ulcorner currentpicture),
   2/3*xpart(urcorner currentpicture-ulcorner currentpicture)));
endpattern;

beginpattern(ttt);
   fill unitsquare scaled 2mm rotated 45 withcolor .5white;
   z1=llcorner currentpicture;
   z2=lrcorner currentpicture;
   z3=urcorner currentpicture;
   z4=ulcorner currentpicture;
%   draw z1--z2--z3--z4--cycle;
   patterntransform(identity rotated 30 slanted .2);
   patternstep((1.1*xpart(urcorner currentpicture-ulcorner currentpicture),
   2/3*xpart(urcorner currentpicture-ulcorner currentpicture)));
endpattern;

beginpattern(tr);
   u:=4mm;
   z1=(0,0);
   z2=(u,0);
   z3=(u,2u);
   z4=(0,2u);
   z12=.5[z1,z2];
   z23=.5[z2,z3];
   z34=.5[z3,z4];
   z41=.5[z4,z1];
   draw z12--z23--z34--z41--cycle;
   draw z23--z41;
   draw z1--z2;
   draw z3--z4;
   patternbbox(z1,z3);
   patterntransform(identity slanted .1 rotated 11);
   patternstep(1.2(x2-x1),1.2(y3-y1));
endpattern;

beginfig(12);
   path p;
   z1=(1cm,1cm);
   p=fullcircle scaled 3cm;
   fill p withpattern tt;
   draw p;
   unfill p shifted z1;
   fill p shifted z1 withpattern tr;
   draw p shifted z1;
endfig;

beginpattern(ll);
   pickup pencircle scaled .3;
   draw origin--4up;
   patternbbox(left,4up+right);
   patterntransform(identity rotated -45);
endpattern;

beginpattern(lll);
   pickup pencircle scaled .3;
   draw origin--4up;
   patternbbox(left,4up+right);
   patterntransform(identity rotated 45);
endpattern;

beginpattern(literki);
   label(btex \font\q=ptmr8r at 5pt\q abcd etex,origin);
   z1=ulcorner currentpicture;
   z2=urcorner currentpicture;
   patterntransform(identity slanted .17);
endpattern;

beginpattern(cm_literki);
   label(btex \font\q=cmr5\q xyz XYZ ZYX zyx etex,origin);
   z1=ulcorner currentpicture;
   z2=urcorner currentpicture;
   patterntransform(identity slanted .17);
endpattern;

beginfig(13);
   label(btex \font\q=ptmr8r at 90pt\q ABCD etex, origin) withpattern ll;
endfig;

beginfig(14);
   picture a; a=btex \font\q=ptmbi8r at 90pt\q ABCD etex;
   draw a withpattern literki;
   undraw a shifted (20up+10right);
   draw a shifted (20up+10right) withpattern checker;
   undraw a shifted (40up+20right);
   draw a shifted (40up+20right) withpattern lll;
   draw a shifted (40up+20right) withpattern ll;
   undraw a shifted (60up+30right);
   draw a shifted (60up+30right) withpattern cm_literki;
endfig;

end.
