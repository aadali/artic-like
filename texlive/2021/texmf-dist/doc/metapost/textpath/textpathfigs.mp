input latexmp
setupLaTeXMP(preamblefile="textpathfigs");
input textpath

prologues := 0;

beginfig(1);
  draw textpath("Greetings from MetaPost!", fullcircle scaled 50bp, 0)
    withcolor (0.2,0.2,0.6);
endfig;


beginfig(2);
path p;
  p := reverse fullcircle rotated -90 scaled 50bp;
  draw textpath("Greetings from MetaPost!", p, 0.5) withcolor (0.2,0.2,0.6);
endfig;


beginfig(3);
path p;
  p := reverse fullcircle rotated -90;
  draw textpathFont("\usefont{T1}{pzc}{m}{n}\huge", "Happy Birthday to",
    p scaled 400bp, 0.5) withcolor (1, 0.6, 0.2);
  draw textpathFont("\usefont{T1}{bch}{m}{n}\large", "Daisy Duck!",
    p scaled 350bp, 0.5) withcolor (0.9, 0.3, 0.1);
endfig;


beginfig(4);
path p;
  p := reverse fullcircle rotated -90 scaled 50bp;
  draw p withcolor .8white;
  draw textpathRaw("Greetings from MetaPost!", p, 0.5) withcolor (0.2,0.2,0.6);
endfig;


beginfig(5);
path p;
  p := reverse halfcircle scaled 80bp;
  draw p withcolor .8white;
  draw textpathRaw("\LARGE$x^2+y^2=r^2$", p, 0.5) withcolor (0.2,0.2,0.6);
endfig;


beginfig(6);
path p;
  p := reverse halfcircle xscaled 150bp yscaled 80;
  draw p withcolor .8white;
  interim textpathShift := 7pt;
  draw textpathRaw("\Large$\frac{(x-u)^2}{a^2}+\frac{(y-v)^2}{b^2}=1$", p, 0.5) withcolor (0.2,0.2,0.6);
endfig;



beginfig(7);
path p;
string f;
  f := "\usefont{T1}{pzc}{m}{n}\Large";
  p := subpath (5.7,6.3) of fullcircle scaled 1400bp;
  draw p withcolor .8white;
  interim textpathRepeat := 3;
  interim textpathStretch := 0;
  interim textpathHSpace := 10pt;
  draw textpathFont(f, "Happy Birthday!", p, 0.5) withcolor (1, 0.6, 0.2);
endfig;


beginfig(8);
path p[];
string f;
  f := "\usefont{T1}{pzc}{m}{n}\Large";
  p1 := subpath (5.7,6.3) of fullcircle scaled 1400bp;
  p2 := p1 shifted (0bp,-30bp);
  draw p1 withcolor .8white;
  draw p2 withcolor .8white;
  interim textpathRepeat := 3;
  interim textpathStretch := 1;
  draw textpathFont(f, "Happy Birthday!", p1, 0.5) withcolor (1, 0.6, 0.2);
  draw textpathFont(f, "Happy Birthday!", p2, 0) withcolor (1, 0.6, 0.2);
endfig;


beginfig(9);
string f;
path p;
numeric r,h,w;
picture pic;
% Font Brush Script Italic is available on CTAN.
  f := "\usefont{T1}{pbsi}{xl}{n}\fontsize{2.1pt}{2.1pt}\selectfont";
  w := 210bp;
  h := .276w;
  r := .19h;
  p := (-.5w,0)--(-.5w,.5h-r)--quartercircle rotated -90 scaled (2r)
    shifted (-.5w,.5h)--(0,.5h);
  p := p--reverse p reflectedabout ((0,-1),(0,1));
  p := p--reverse p reflectedabout ((-1,0),(1,0))--cycle;
  draw p withpen pensquare scaled 3.5pt;
  interim textpathRepeat := 30;
  interim textpathShift := -.6pt;
  draw textpathFont(f, "Fa\ss' Dich kurz!", p, 0.5) withcolor white;
  label(textext
    ("\usefont{T1}{bfu}{mb}{n}\fontsize{22pt}{22pt}\selectfont Telephonzelle"),
    origin);% Bitstream Futura
endfig;


beginfig(10);
string f;
path p;
picture pic;
numeric lenp, lenpic, tb, k;
color daisy;
  daisy := (0.9, 0.3, 0.1);
  f := "\usefont{T1}{bch}{m}{n}\fontsize{48pt}{48pt}\selectfont";
  p := subpath (2.7,5.3) of reverse fullcircle rotated -90 scaled 120bp;
  draw p withcolor .8white;
  interim textpathAbsRotation := 1;
  draw textpathFont(f, "Daisy", p, 0.5) withcolor daisy + 0.85(white-daisy);
  textpathAbsRotation := 0;
  interim textpathDebug := 2;
  pickup pencircle scaled .2bp;
  draw textpathFont(f, "Daisy", p, 0.5) withcolor .6white;
  textpathDebug := 1;
  draw textpathFont(f, "Daisy", p, 0.5) withcolor daisy;
  textpathDebug := 4;
  draw textpathFont(f, "Daisy", p, 0.5) withcolor .6white;
endfig;


beginfig(11);
string f, t;
path heart;
  f := "\usefont{T1}{pzc}{m}{it}\tiny";
  t := "Love";
  heart := halfcircle shifted (-0.5bp,0bp)..{dir-50}(0bp,-1.5bp);
  heart := heart--reverse heart reflectedabout ((0,0),(0,1))--cycle;
  heart := heart scaled 60bp;
  interim textpathRotation := 90;
  interim textpathLetterSpace := 1pt;
  interim textpathRepeat := 30;
  draw textpathFont(f, t, heart, 0.5) withcolor red+0.1green;
endfig;


beginfig(12);
path p;
string f, t;
picture pic;
  f := "\usefont{T1}{fwb}{m}{n}\Large";% From the emerald package
  t := "Don't panic!  Don't panic!  Don't panic!";
  p := origin
  for i:=1 upto 20: ..(i, sind(i*45)) endfor;
  p := p xscaled 20 yscaled 35;
  interim textpathRotation := 0;
  interim textpathAbsRotation := 1;
  interim textpathLetterSpace := 6pt;
  pic := textpathFont(f, t, p, 0);
  background := black;
  interim bboxmargin := 30bp;
  unfill bbox pic;
  draw pic withcolor red+green;
endfig;


beginfig(13);
path p;
  p := (origin--(1bp,0bp)) scaled 75bp;
  draw p withcolor .8white;
  interim textpathClip := 1;
  draw textpath("Greetings from MetaPost!", p, 0.5) withcolor (0.2,0.2,0.6);
  p := p shifted (0,-20bp);
  draw p withcolor .8white;
  textpathClip := 0;
  draw textpath("Greetings from MetaPost!", p, 0.5) withcolor (0.2,0.2,0.6);
  p := p shifted (0,-20bp);
  draw p withcolor .8white;
  interim textpathAutoScale := 1;
  draw textpath("Greetings from MetaPost!", p, 0.5) withcolor (0.2,0.2,0.6);
endfig;


beginfig(14);
path p;
  interim textpathDraft := 1;
  p := (origin--(1bp,0bp)) scaled 75bp;
  draw p withcolor .8white;
  interim textpathClip := 1;
  draw textpath("Greetings from MetaPost!", p, 0.5) withcolor (0.2,0.2,0.6);
  p := p shifted (0,-20bp);
  draw p withcolor .8white;
  textpathClip := 0;
  draw textpath("Greetings from MetaPost!", p, 0.5) withcolor (0.2,0.2,0.6);
  p := p shifted (0,-20bp);
  draw p withcolor .8white;
  interim textpathAutoScale := 1;
  draw textpath("Greetings from MetaPost!", p, 0.5) withcolor (0.2,0.2,0.6);
endfig;


beginfig(15);
path p;
  p := reverse halfcircle xscaled 150bp yscaled 80;
  draw p withcolor .8white;
  interim textpathShift := 9pt;
  interim textpathFancyStrokes := 0;
  draw textpathRaw("\Large$y = v+b\sqrt{1-\frac{(x-u)^2}{a^2}}$", p, 0.5) withcolor (0.2,0.2,0.6);
endfig;


beginfig(16);
path p;
  p := reverse halfcircle xscaled 150bp yscaled 80;
  draw p withcolor .8white;
  interim textpathShift := 9pt;
  interim textpathFancyStrokes := 1;
  interim textpathCureSqrt := 0;
  draw textpathRaw("\Large$y = v+b\sqrt{1-\frac{(x-u)^2}{a^2}}$", p, 0.5) withcolor (0.2,0.2,0.6);
endfig;


beginfig(17);
path p;
  p := reverse halfcircle xscaled 150bp yscaled 80;
  draw p withcolor .8white;
  interim textpathShift := 9pt;
  interim textpathFancyStrokes := 1;
  interim textpathCureSqrt := 1;
  draw textpathRaw("\Large$y = v+b\sqrt{1-\frac{(x-u)^2}{a^2}}$", p, 0.5) withcolor (0.2,0.2,0.6);
endfig;


beginfig(18);
string f;
path p;
  f := "\usefont{T1}{pzc}{m}{it}\fontsize{56pt}{56pt}\selectfont";
  p := reverse fullcircle rotated 33 scaled 150bp;
  interim textpathLetterSpace := -3pt;
  draw textpathFont(f, "Pro \&", p, 1) withcolor (1, 0.6, 0.2);
  interim textpathLetterSpace := 0pt;
  f := "\usefont{T1}{pzc}{m}{it}\fontsize{48pt}{48pt}\selectfont";
  p := fullcircle rotated 33 scaled 145bp;
  draw textpathFont(f, "Con\<{\kern-2pt}tra", p, 0.0) withcolor (1, 0.6, 0.2);
endfig;


beginfig(19);
path p;
  p := reverse fullcircle rotated -90 scaled 50bp;
  interim textpathDebug := 1 + 2;
  draw textpath("Greetings from MetaPost!", p, 0.5) withcolor (0.2,0.2,0.6);
endfig;


beginfig(20);
path p;
  p := reverse fullcircle rotated -90 scaled 50bp;
  interim textpathDebug := 2 + 4;
  draw textpath("Greetings from MetaPost!", p, 0.5) withcolor (0.2,0.2,0.6);
endfig;


end
