%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Structure of a MetaPost file using drv
%   1. Preamble
%        input drv;
%        verbatimtex%&latex
%        <LaTeX preamble>
%        \begin{document}
%        etex;
%  2. Figures
%       <optional drv tunings>
%       beginfig(<index>)
%       <judgment & inference declarations>
%       draw drv_tree;
%       <optional extra MetaPost code>
%       endfig;
%  3. Postamble
%       end
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1. Preamble
%

input drv;

verbatimtex%&latex

% Any piece of code related to font selection in the preamble of
% `drv-guide.tex' (document class, font package, font selection commands ...)
% should be reported here.

\documentclass[twoside,11pt]{article}
\usepackage[T1]{fontenc}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{stmaryrd}

% Choose your favourite fonts combination (if available).
% See "http://ctan.tug.org/tex-archive/info/Free_Math_Font_Survey".
%
%\usepackage{cmbright}			% CM text & CM Bright math
%\usepackage{ccfonts,eulervm}		% Concrete text & Euler math
%\usepackage{ccfonts}			% Concrete text & math
%\usepackage[math]{iwona}		% Iwona text & math
%\usepackage[math]{kurier}		% Kurier text & math
%\usepackage[math]{anttor}		% Antykwa Torunska text & math
%\usepackage{kmath,kerkis}		% Kerkis text & math
%\usepackage{fouriernc}			% New Century Schoolbook & Fourier math
%\usepackage{pxfonts}			% Palatino & pxfonts math
%\usepackage{mathpazo}			% Palatino & Pazo math
%\usepackage{mathpple}			% Palatino & Euler math
\usepackage[varg]{txfonts}		% Times & txfonts math
%\usepackage{mathtime}			% Times & Belleek math
%\usepackage{mathptmx}			% Times & Symbol math
%\usepackage{mbtimes}			% Omega Serif text & Omega math
%\usepackage{arev}			% Arev Sans text with Arev math
%\usepackage[charter]{mathdesign}	% Bitstream Charter & Math Design math
%\usepackage{comicsans}			% Comic Sans text & math
%\usepackage[garamond]{mathdesign}	% Garamond & Math Design math
%\usepackage{fourier}			% Utopia & Fourier math
%\usepackage[utopia]{mathdesign}	% Utopia & Math Design math

\newcommand{\limp}{\multimapinv}
\newcommand{\rimp}{\multimap}

\begin{document}
etex;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 2. Figures
%

%
% Uncomment the following line to experiment the radial mode.
% drv_radial_mode on;
%

%
% title page
%

beginfig(100)
jgm 0 "\Gamma, \Delta, \Theta, \Pi, \Upsilon\vdash F";
  jgm 1 "\Pi\vdash (A\to D)\to E";
  jgm 2 "\Gamma, \Delta, \Theta, (A\to D)\to E, \Upsilon\vdash F";
    jgm 3 "\Gamma, \Delta, \Theta\vdash A\to D";
      jgm 4 "A, \Gamma, \Delta, \Theta\vdash D";
        jgm 5 "A, \Gamma, \Delta\vdash B\wedge C";
          jgm 6 "A, \Gamma\vdash B";
          jgm 7 "\Delta\vdash C";
        jgm 8 "B\wedge C, \Theta\vdash D";
    jgm 9 "E, \Upsilon\vdash F";
nfr 0 (1, 2) ("\text{cut}", 1);
  nfr 1 () ("\pi", 4);
  nfr 2 (3, 9) ("\to_L", 1);
    nfr 3 (4) ("\to_R", 1);
      nfr 4 (5, 8) ("\text{cut}", 1);
        nfr 5 (6, 7) ("\wedge_R", 1);
          nfr 6 () ("\gamma", 2);
          nfr 7 () ("\delta", 1);
        nfr 8 () ("\theta", 3);
    nfr 9 () ("\upsilon", 2);
draw drv_tree;
endfig;

%
% first example
%

beginfig(110)
jgm 0 "A\vdash B";
jgm 1 "B\vdash C";
jgm 2 "A\vdash C";
nfr 0 () ("f", 1);
nfr 1 () ("g", 1);
nfr 2 (0, 1) ("\circ", 1);
draw drv_tree;
endfig;

%
% sub-judgments
%

beginfig(111)
jgm 0 "A\vdash B";
jgm 1 "A", "\vdash", "B";
nfr 0 () ("f", 1);
nfr 1 (0) ("f", 1);
draw drv_tree;
endfig;


%
% inference line styles
%

beginfig(120)
jgm 0 "\text{none}";
jgm 1 "\text{simple}";
jgm 2 "\text{double}";
jgm 3 "\text{dotted}";
jgm 4 "\text{dashed}";
jgm 5 "\text{waved}";
jgm 6 "\text{\TeX-dotted}";
jgm 7 "\text{default}";
nfr 0 ( ) ("\leftarrow 0", 0);
nfr 1 (0) ("\leftarrow 1", 1);
nfr 2 (1) ("\leftarrow 2", 2);
nfr 3 (2) ("\leftarrow 3", 3);
nfr 4 (3) ("\leftarrow 4", 4);
nfr 5 (4) ("\leftarrow 5", 5);
nfr 6 (5) ("\leftarrow 6", 6);
nfr 7 (6) ("\leftarrow\_", _);
draw drv_tree;
endfig;

%
% declarations order
%

beginfig(130)
% preorder declarations
jgm 0 "0";
  jgm 1 "00";
    jgm 2 "000";
    jgm 3 "001";
    jgm 4 "002";
  jgm 5 "01";
    jgm 6 "010";
    jgm 7 "011";
    jgm 8 "012";
  jgm 9 "02";
    jgm 10 "020";
    jgm 11 "021";
    jgm 12 "022";
nfr 0 (1, 5, 9) ("a", _);
  nfr 1 (2, 3, 4) ("b", _);
    nfr 2 () ("c", _);
    nfr 3 () ("d", _);
    nfr 4 () ("e", _);
  nfr 5 (6, 7, 8) ("f", _);
    nfr 6 () ("g", _);
    nfr 7 () ("h", _);
    nfr 8 () ("i", _);
  nfr 9 (10, 11, 12) ("j", _);
    nfr 10 () ("k", _);
    nfr 11 () ("l", _);
    nfr 12 () ("m", _);
draw drv_tree;
endfig;

beginfig(131)
% postorder declarations
    jgm 0 "000";
    jgm 1 "001";
    jgm 2 "002";
  jgm 3 "00";
    jgm 4 "010";
    jgm 5 "011";
    jgm 6 "012";
  jgm 7 "01";
    jgm 8 "020";
    jgm 9 "021";
    jgm 10 "022";
  jgm 11 "02";
jgm 12 "0";
    nfr 0 () ("a", _);
    nfr 1 () ("b", _);
    nfr 2 () ("c", _);
  nfr 3 (0, 1, 2) ("d", _);
    nfr 4 () ("e", _);
    nfr 5 () ("f", _);
    nfr 6 () ("g", _);
  nfr 7 (4, 5, 6) ("h", _);
    nfr 8 () ("i", _);
    nfr 9 () ("j", _);
    nfr 10 () ("k", _);
  nfr 11 (8, 9, 10) ("l", _);
nfr 12 (3, 7, 11) ("m", _);
draw drv_tree;
endfig;

%
% dcl
%

beginfig(140)
dcl 0 () ("f", 1) "A\vdash B";
dcl 1 () ("g", 1) "B\vdash C";
dcl 2 (0, 1) ("\circ", 1) "A\vdash C";
draw drv_tree;
endfig;

beginfig(141)
dcl 0 () ("f", 1) "A", "\vdash", "B";
draw drv_tree;
endfig;

beginfig(141)
dcl 0 () ("f", 1) "A\vdash B";
dcl 1 (0) ("f", 1) "A", "\vdash", "B";
draw drv_tree;
endfig;

beginfig(150)
dcl 0 (1, 5, 9) ("a", _)     "0";
  dcl 1 (2, 3, 4) ("b", _)     "00";
    dcl 2 () ("c", _)            "000";
    dcl 3 () ("d", _)            "001";
    dcl 4 () ("e", _)            "002";
  dcl 5 (6, 7, 8) ("f", _)     "01";
    dcl 6 () ("g", _)            "010";
    dcl 7 () ("h", _)            "011";
    dcl 8 () ("i", _)            "012";
  dcl 9 (10, 11, 12) ("j", _)  "02";
    dcl 10 () ("k", _)            "020";
    dcl 11 () ("l", _)            "021";
    dcl 12 () ("m", _)            "022";
draw drv_tree;
endfig;

%
% bxd
%

beginfig(160)
dcl 0 (1, 4) ("", _)     "a";
  dcl 1 (2) ("", _)      "a";
    dcl 2 (3) ("", _)    "a";
      dcl 3 () ("", _)   "aaaaaaa";
  dcl 4 () ("", _)       "aaaaa";
draw drv_tree;
endfig;

beginfig(161)
dcl 0 (bxd 1, 4) ("", _) "a";
  dcl 1 (2) ("", _)      "a";
    dcl 2 (3) ("", _)    "a";
      dcl 3 () ("", _)   "aaaaaaa";
  dcl 4 () ("", _)       "aaaaa";
draw drv_tree;
endfig;

vardef auxbox[]=
  if drv_radialmode_on:
    save l_ngl, r_ngl, min_rad, max_rad, tnum;
    numeric l_ngl, r_ngl, min_rad, max_rad, tnum;
    l_ngl:=ngl_part sub_lmost_ppt[@];
    r_ngl:=ngl_part sub_rmost_ppt[@];
    tnum:=if bot_jdg[@]: - fi 1;
    min_rad:=iln[ccl[@]].rad+iln[ccl[@]].hg+min_clr_scl*min_clr[drv_sty]/2;
    max_rad:=sub_max_rad[@]+min_clr_scl*min_clr[drv_sty]/2;
    ((dir l_ngl) scaled min_rad--(dir l_ngl) scaled max_rad
       {dir (l_ngl+tnum*90)} .. {dir (r_ngl+tnum*90)}
     (dir r_ngl) scaled max_rad--(dir r_ngl) scaled min_rad
       {dir (r_ngl-tnum*90)} .. {dir (l_ngl-tnum*90)}
     cycle) shifted jdg[@].org
  else:
    ((sub_min_x[@], ypart iln[ccl[@]].c+iln[ccl[@]].hg)--
     (sub_min_x[@], sub_max_y[@])--
     (sub_max_x[@], sub_max_y[@])--
     (sub_max_x[@], ypart iln[ccl[@]].c+iln[ccl[@]].hg)--
     cycle) shifted (0, min_clr_scl*min_clr[drv_sty]/2)
  fi
enddef;

beginfig(162)
dcl 0 (bxd 1, 4) ("", _) "a";
  dcl 1 (2) ("", _)	 "a";
    dcl 2 (3) ("", _)	 "a";
      dcl 3 () ("", _)	 "aaaaaaa";
  dcl 4 () ("", _)	 "aaaaa";
draw drv_tree;
draw auxbox[1] withcolor (1, 0, 0);
endfig;

%
% mvd
%

beginfig(170)
jgm 1 "aaa";
jgm 2 "bbb";
jgm 3 "ccc";
jgm 4 "d";
nfr 1 () ("", _);
nfr 2 () ("", _);
nfr 3 () ("", _);
nfr 4 (mvd 1 (2, 3), 2, 3) ("", _);
draw drv_tree;
endfig;

beginfig(171)
jgm 1 "aaa";
jgm 2 "bbb";
jgm 3 "ccc";
jgm 4 "d";
nfr 1 () ("", _);
nfr 2 () ("", _);
nfr 3 () ("", _);
nfr 4 (1, mvd 2 (2, 4), 3) ("", _);
draw drv_tree;
endfig;

beginfig(172)
jgm 1 "aaa";
jgm 2 "bbb";
jgm 3 "ccc";
jgm 4 "d";
nfr 1 () ("", _);
nfr 2 () ("", _);
nfr 3 () ("", _);
nfr 4 (1, 2, mvd 3 (2, 6)) ("", _);
draw drv_tree;
endfig;

beginfig(180)
jgm 0 "\textsc{Size matters -- Part 1}";
jgm 1 "\text{Here is a rather long judgment"& % string concatenation
      " that I don't want to shorten.}";
jgm 2 "\text{Will the derivation tree fit on the page?}";
jgm 3 "\text{It does.}";
nfr 0 () ("", 0);
nfr 1 (0) ("", 1);
nfr 2 () ("", 1);
nfr 3 (mvd 1 (2, 3), 2) ("", 1);
draw drv_tree;
endfig;

%
% Nfr
%

beginfig(190)
jgm 0 "a";
jgm 1 "b";
jgm 2 "c";
jgm 3 "d";
Nfr 0 () ("0", "", "", _);
Nfr 1 () ("1", "", "", _);
Nfr 2 (0, 1) ("2", "E", "", _);
Nfr 3 (2) ("3", "", "F", _);
draw drv_tree;
endfig;

%
% Dcl
%

beginfig(200)
Dcl 0 () ("", "", "", _) "a";
Dcl 1 (0) ("", "", "B", _) "c";
Dcl 2 () ("", "", "", _) "d";
Dcl 3 (1, 2) ("", "E", "", _) "f";
draw drv_tree;
endfig;

%
% Mvd
%

beginfig(210)
jgm 1 "aaa";
jgm 2 "bbb";
jgm 3 "ccc";
nfr 1 () ("", _);
nfr 2 () ("", _);
nfr 3 (Mvd 1 (2, "d", "", 3), 2) ("", _);
draw drv_tree;
endfig;

%
% drv_font_size
%

drv_font_size "\tiny";

beginfig(220)
jgm 1 "a";
jgm 2 "b";
jgm 3 "c";
nfr 1 () ("1", 1);
nfr 2 () ("2", 1);
Nfr 3 (1, 2) ("3", "4", "", 1);
draw drv_tree;
endfig;

drv_font_size "\small";

beginfig(221)
jgm 1 "a";
jgm 2 "b";
jgm 3 "c";
nfr 1 () ("1", 1);
nfr 2 () ("2", 1);
Nfr 3 (1, 2) ("3", "4", "", 1);
draw drv_tree;
endfig;

drv_font_size "\Large";

beginfig(222)
jgm 1 "a";
jgm 2 "b";
jgm 3 "c";
nfr 1 () ("1", 1);
nfr 2 () ("2", 1);
Nfr 3 (1, 2) ("3", "4", "", 1);
draw drv_tree;
endfig;

drv_font_size "\normalsize";

%
% drv_math_style (drv, );
%

beginfig(230)
jgm 1 "a";
jgm 2 "b";
jgm 3 "c";
nfr 1 () ("1", 1);
nfr 2 () ("2", 1);
Nfr 3 (1, 2) ("3", "4", "", 1);
draw drv_tree;
endfig;

drv_math_style (drv, "\textstyle");

beginfig(231)
jgm 1 "a";
jgm 2 "b";
jgm 3 "c";
nfr 1 () ("1", 1);
nfr 2 () ("2", 1);
Nfr 3 (1, 2) ("3", "4", "", 1);
draw drv_tree;
endfig;

drv_math_style (drv, "\scriptstyle");

beginfig(232)
jgm 1 "a";
jgm 2 "b";
jgm 3 "c";
nfr 1 () ("1", 1);
nfr 2 () ("2", 1);
Nfr 3 (1, 2) ("3", "4", "", 1);
draw drv_tree;
endfig;

drv_math_style (drv, "\displaystyle");

%
% drv_math_style (jdg, );
%

drv_math_style (jdg, "\displaystyle");

beginfig(240)
jgm 1 "\bigwedge_{i\in I} A_i";
jgm 2 "b";
jgm 3 "c";
nfr 1 () ("1", 1);
nfr 2 () ("2", 1);
Nfr 3 (1, 2) ("3", "4", "", 1);
draw drv_tree;
endfig;

drv_math_style (jdg, "\textstyle");

beginfig(241)
jgm 1 "\bigwedge_{i\in I} A_i";
jgm 2 "b";
jgm 3 "c";
nfr 1 () ("1", 1);
nfr 2 () ("2", 1);
Nfr 3 (1, 2) ("3", "4", "", 1);
draw drv_tree;
endfig;

drv_math_style (jdg, "\scriptstyle");

beginfig(242)
jgm 1 "\bigwedge_{i\in I} A_i";
jgm 2 "b";
jgm 3 "c";
nfr 1 () ("1", 1);
nfr 2 () ("2", 1);
Nfr 3 (1, 2) ("3", "4", "", 1);
draw drv_tree;
endfig;

drv_math_style (jdg, "\textstyle");

%
% drv_math_style (ilb, );
%

drv_math_style (ilb, "\textstyle");

beginfig(250)
jgm 1 "a";
jgm 2 "b";
jgm 3 "c";
nfr 1 () ("1", 1);
nfr 2 () ("2", 1);
nfr 3 (1, 2) ("3", 1);
draw drv_tree;
endfig;

drv_math_style (ilb, "\scriptstyle");

beginfig(251)
jgm 1 "a";
jgm 2 "b";
jgm 3 "c";
nfr 1 () ("1", 1);
nfr 2 () ("2", 1);
nfr 3 (1, 2) ("3", 1);
draw drv_tree;
endfig;

drv_math_style (ilb, "\scriptscriptstyle");

beginfig(252)
jgm 1 "a";
jgm 2 "b";
jgm 3 "c";
nfr 1 () ("1", 1);
nfr 2 () ("2", 1);
nfr 3 (1, 2) ("3", 1);
draw drv_tree;
endfig;

drv_math_style (ilb, "\scriptstyle");

%
% drv_scale (clr, );
%

drv_scale (clr, 0);

beginfig(260)
dcl 1 () ("", 1) "\Bigl(a\Bigr)";
dcl 2 (1) ("", 1) "a";
draw drv_tree;
endfig;

drv_scale (clr, 1);

beginfig(261)
dcl 1 () ("", 1) "\Bigl(a\Bigr)";
dcl 2 (1) ("", 1) "a";
draw drv_tree;
endfig;

drv_scale (clr, 2.5);

beginfig(262)
dcl 1 () ("", 1) "\Bigl(a\Bigr)";
dcl 2 (1) ("", 1) "a";
draw drv_tree;
endfig;

drv_scale (clr, 4);

beginfig(263)
dcl 1 () ("", 1) "\Bigl(a\Bigr)";
dcl 2 (1) ("", 1) "a";
draw drv_tree;
endfig;

drv_scale (clr, 1);

%
% drv_scale (prm, );
%

drv_scale (prm, 0);

beginfig(270)
dcl 1 () ("", 1) "a";
dcl 2 () ("", 1) "a";
dcl 3 (1, 2) ("", 1) "a";
draw drv_tree;
endfig;

drv_scale (prm, 1);

beginfig(271)
dcl 1 () ("", 1) "a";
dcl 2 () ("", 1) "a";
dcl 3 (1, 2) ("", 1) "a";
draw drv_tree;
endfig;

drv_scale (prm, 2.5);

beginfig(272)
dcl 1 () ("", 1) "a";
dcl 2 () ("", 1) "a";
dcl 3 (1, 2) ("", 1) "a";
draw drv_tree;
endfig;

drv_scale (prm, 4);

beginfig(273)
dcl 1 () ("", 1) "a";
dcl 2 () ("", 1) "a";
dcl 3 (1, 2) ("", 1) "a";
draw drv_tree;
endfig;

drv_scale (prm, 1);

%
% drv_scale (jgm, );
%

drv_scale (jgm, 0);

beginfig(280)
dcl 1 () ("", 1) "a";
dcl 2 () ("", 1) "a";
dcl 3 (1, 2) ("", 1) "a";
draw drv_tree;
endfig;

drv_scale (jgm, 1);

beginfig(281)
dcl 1 () ("", 1) "a";
dcl 2 () ("", 1) "a";
dcl 3 (1, 2) ("", 1) "a";
draw drv_tree;
endfig;

drv_scale (jgm, 2.5);

beginfig(282)
dcl 1 () ("", 1) "a";
dcl 2 () ("", 1) "a";
dcl 3 (1, 2) ("", 1) "a";
draw drv_tree;
endfig;

drv_scale (jgm, 4);

beginfig(283)
dcl 1 () ("", 1) "a";
dcl 2 () ("", 1) "a";
dcl 3 (1, 2) ("", 1) "a";
draw drv_tree;
endfig;

drv_scale (jgm, 1);

%
% drv_scale (ilb, );
%

drv_scale (ilb, 0);

beginfig(290)
dcl 1 () ("b", 1) "a";
dcl 2 () ("b", 1) "a";
dcl 3 (1, 2) ("b", 1) "a";
draw drv_tree;
endfig;

drv_scale (ilb, 1);

beginfig(291)
dcl 1 () ("b", 1) "a";
dcl 2 () ("b", 1) "a";
dcl 3 (1, 2) ("b", 1) "a";
draw drv_tree;
endfig;

drv_scale (ilb, 2.5);

beginfig(292)
dcl 1 () ("b", 1) "a";
dcl 2 () ("b", 1) "a";
dcl 3 (1, 2) ("b", 1) "a";
draw drv_tree;
endfig;

drv_scale (ilb, 4);

beginfig(293)
dcl 1 () ("b", 1) "a";
dcl 2 () ("b", 1) "a";
dcl 3 (1, 2) ("b", 1) "a";
draw drv_tree;
endfig;

drv_scale (ilb, 1);

%
% drv_junction_style
%

drv_junction_style 0;

beginfig(300)
dcl 1 ()	("", 1) "aaaaaaaaaaaa";
dcl 2 (1)	("", 1) "a";
dcl 3 (2)	("", 1) "a";
dcl 4 (3)	("", 1) "a";
dcl 5 (4)	("", 1) "aaaaaa";
dcl 6 ()	("", 1) "aaaaaa";
dcl 7 (6)	("", 1) "a";
dcl 8 (7)	("", 1) "a";
dcl 9 (5, 8)	("", 1) "a";
draw drv_tree;
endfig;

drv_junction_style 1;

beginfig(301)
dcl 1 ()	("", 1) "aaaaaaaaaaaa";
dcl 2 (1)	("", 1) "a";
dcl 3 (2)	("", 1) "a";
dcl 4 (3)	("", 1) "a";
dcl 5 (4)	("", 1) "aaaaaa";
dcl 6 ()	("", 1) "aaaaaa";
dcl 7 (6)	("", 1) "a";
dcl 8 (7)	("", 1) "a";
dcl 9 (5, 8)	("", 1) "a";
draw drv_tree;
draw auxbox[8] withcolor .9(1, 1, 1);
endfig;

drv_junction_style 2;

beginfig(302)
dcl 1 ()	("", 1) "aaaaaaaaaaaa";
dcl 2 (1)	("", 1) "a";
dcl 3 (2)	("", 1) "a";
dcl 4 (3)	("", 1) "a";
dcl 5 (4)	("", 1) "aaaaaa";
dcl 6 ()	("", 1) "aaaaaa";
dcl 7 (6)	("", 1) "a";
dcl 8 (7)	("", 1) "a";
dcl 9 (5, 8)	("", 1) "a";
draw drv_tree;
draw auxbox[5] withcolor .9(1, 1, 1);
draw auxbox[8] withcolor .9(1, 1, 1);
endfig;

drv_junction_style 1;

%
% drv_alignment_style
%

drv_alignment_style l;

beginfig(310)
dcl 1 () ("", 1) "a";
dcl 2 () ("", 1) "a";
dcl 3 (1, 2) ("", 1) "a";
dcl 4 () ("", 1) "a";
dcl 5 () ("", 1) "a";
dcl 6 (4, 5) ("", 1) "a";
dcl 7 (3, 6) ("", 1) "a";
draw drv_tree;
endfig;

drv_alignment_style c;

beginfig(311)
dcl 1 () ("", 1) "a";
dcl 2 () ("", 1) "a";
dcl 3 (1, 2) ("", 1) "a";
dcl 4 () ("", 1) "a";
dcl 5 () ("", 1) "a";
dcl 6 (4, 5) ("", 1) "a";
dcl 7 (3, 6) ("", 1) "a";
draw drv_tree;
endfig;

drv_alignment_style r;

beginfig(312)
dcl 1 () ("", 1) "a";
dcl 2 () ("", 1) "a";
dcl 3 (1, 2) ("", 1) "a";
dcl 4 () ("", 1) "a";
dcl 5 () ("", 1) "a";
dcl 6 (4, 5) ("", 1) "a";
dcl 7 (3, 6) ("", 1) "a";
draw drv_tree;
endfig;

drv_alignment_style c;

%
% drv_labels_position
%

drv_labels_position l;

beginfig(320)
dcl 1 () ("b", 1) "a";
dcl 2 () ("b", 1) "a";
dcl 3 (1, 2) ("b", 1) "a";
draw drv_tree;
endfig;

drv_labels_position r;

beginfig(321)
dcl 1 () ("b", 1) "a";
dcl 2 () ("b", 1) "a";
dcl 3 (1, 2) ("b", 1) "a";
draw drv_tree;
endfig;

%
% drv_roots_position
%

drv_roots_position t;

beginfig(330)
dcl 1 () ("", 1) "a";
dcl 2 () ("", 1) "a";
dcl 3 (1, 2) ("", 1) "a";
dcl 4 () ("", 1) "a";
dcl 5 () ("", 1) "a";
dcl 6 (4, 5) ("", 1) "a";
dcl 7 (3, 6) ("", 1) "a";
draw drv_tree;
drv_axis (iln, 3);
endfig;

drv_roots_position b;

beginfig(331)
dcl 1 () ("", 1) "a";
dcl 2 () ("", 1) "a";
dcl 3 (1, 2) ("", 1) "a";
dcl 4 () ("", 1) "a";
dcl 5 () ("", 1) "a";
dcl 6 (4, 5) ("", 1) "a";
dcl 7 (3, 6) ("", 1) "a";
draw drv_tree;
endfig;

%
% drv_axis_reference
%

beginfig(340)
dcl 1 () ("", 1) "a";
dcl 2 () ("", 1) "a";
dcl 3 (1, 2) ("", 1) "a";
draw drv_tree;
endfig;

drv_axis_reference jdg;

beginfig(341)
dcl 1 () ("", 1) "a";
dcl 2 () ("", 1) "a";
dcl 3 (1, 2) ("", 1) "a";
draw drv_tree;
endfig;

drv_axis_reference iln;

%
% drv_left_delimiter
%

drv_left_delimiter "(";

beginfig(350)
dcl 1 () ("", 1) "a";
Dcl 2 (1) ("", "", "B", 1) "c";
dcl 3 () ("", 1) "d";
Dcl 4 (2, 3) ("", "E", "", 1) "f";
draw drv_tree;
endfig;

drv_left_delimiter "\lfloor";

beginfig(351)
dcl 1 () ("", 1) "a";
Dcl 2 (1) ("", "", "B", 1) "c";
dcl 3 () ("", 1) "d";
Dcl 4 (2, 3) ("", "E", "", 1) "f";
draw drv_tree;
endfig;

drv_left_delimiter ".";

beginfig(352)
dcl 1 () ("", 1) "a";
Dcl 2 (1) ("", "", "B", 1) "c";
dcl 3 () ("", 1) "d";
Dcl 4 (2, 3) ("", "E", "", 1) "f";
draw drv_tree;
endfig;

drv_left_delimiter "\lbrace";

%
% drv_right_delimiter
%

drv_right_delimiter "\rangle";

beginfig(360)
dcl 1 () ("", 1) "a";
Dcl 2 (1) ("", "", "B", 1) "c";
dcl 3 () ("", 1) "d";
Dcl 4 (2, 3) ("", "E", "", 1) "f";
draw drv_tree;
endfig;

drv_right_delimiter "\uparrow";

beginfig(361)
dcl 1 () ("", 1) "a";
Dcl 2 (1) ("", "", "B", 1) "c";
dcl 3 () ("", 1) "d";
Dcl 4 (2, 3) ("", "E", "", 1) "f";
draw drv_tree;
endfig;

drv_right_delimiter ".";

beginfig(362)
dcl 1 () ("", 1) "a";
Dcl 2 (1) ("", "", "B", 1) "c";
dcl 3 () ("", 1) "d";
Dcl 4 (2, 3) ("", "E", "", 1) "f";
draw drv_tree;
endfig;

drv_right_delimiter "\rbrace";

%
% drv_box_mode
%

drv_box_mode on;

beginfig(370)
dcl 0 ()	("", 1) "aaaaa";
dcl 1 ()	("", 1) "aaaaaaa";
dcl 2 (1)	("", 1) "a";
dcl 3 (2)	("", 1) "a";
dcl 4 (0, 3)	("", 1) "a";
draw drv_tree;
endfig;

beginfig(371)
dcl 0 ()	("", 1) "aaaaa";
dcl 1 ()	("", 1) "aaaaaaa";
dcl 2 (1)	("", 1) "a";
dcl 3 (2)	("", 1) "a";
dcl 4 (0, 3)	("", 1) "a";
draw drv_tree;
for idx=0, 1, 2, 3:
  draw auxbox[idx] withcolor (1, 0, 0);
endfor
endfig;

drv_box_mode off;

beginfig(372)
dcl 0 ()	("", 1) "aaaaa";
dcl 1 ()	("", 1) "aaaaaaa";
dcl 2 (1)	("", 1) "a";
dcl 3 (2)	("", 1) "a";
dcl 4 (0, 3)	("", 1) "a";
draw drv_tree;
endfig;

%
% drv_fraction_mode
%

beginfig(380)
jgm 0 "\overbrace{A, \Gamma\vdash B}";
jgm 1 "B, \Delta\vdash C";
jgm 2 "A, \Gamma, \Delta\vdash C";
jgm 3 "\Gamma, \Delta\vdash A\to C";
nfr 0 () ("", 0);
nfr 1 () ("", 1);
nfr 2 (0, 1) ("\text{cut}", 1);
nfr 3 (2) ("\to_R", 1);
draw drv_tree;
endfig;

drv_fraction_mode off;

beginfig(381)
jgm 0 "\overbrace{A, \Gamma\vdash B}";
jgm 1 "B, \Delta\vdash C";
jgm 2 "A, \Gamma, \Delta\vdash C";
jgm 3 "\Gamma, \Delta\vdash A\to C";
nfr 0 () ("", 0);
nfr 1 () ("", 1);
nfr 2 (0, 1) ("\text{cut}", 1);
nfr 3 (2) ("\to_R", 1);
draw drv_tree;
endfig;

drv_fraction_mode on;

beginfig(382)
jgm 0 "\overbrace{A, \Gamma\vdash B}";
jgm 1 "B, \Delta\vdash C";
jgm 2 "A, \Gamma, \Delta\vdash C";
jgm 3 "\Gamma, \Delta\vdash A\to C";
nfr 0 () ("", 1);
nfr 1 () ("", 1);
nfr 2 (0, 1) ("\text{cut}", 1);
nfr 3 (2) ("\to_R", 1);
draw drv_tree;
draw iln[0] withcolor background withpen drv_pen scaled 2;
draw iln[0] withcolor (1, 0, 0);
endfig;

%
% drv_proof_mode
%

drv_proof_mode on;

beginfig(390)
jgm 0 "A", "\vdash", "A";
jgm 1 "B", "\vdash", "B";
jgm 2 "A", ",", "A", "\multimap", "B", "\vdash", "B";
jgm 3 "A", "\multimap", "B", "\vdash", "A", "\multimap", "B";
nfr 0 () ("1", 1);
nfr 1 () ("1", 1);
nfr 2 (0, 1) ("\multimap_{L}", 1);
nfr 3 (2) ("\multimap_{R}", 1);
draw drv_tree;
endfig;

drv_proof_mode off;

beginfig(391)
jgm 0 "A", "\vdash", "A";
jgm 1 "B", "\vdash", "B";
jgm 2 "A", ",", "A", "\multimap", "B", "\vdash", "B";
jgm 3 "A", "\multimap", "B", "\vdash", "A", "\multimap", "B";
nfr 0 () ("1", 1);
nfr 1 () ("1", 1);
nfr 2 (0, 1) ("\multimap_{L}", 1);
nfr 3 (2) ("\multimap_{R}", 1);
draw drv_tree;
endfig;

%
% a fraction
%

drv_font_size "\huge";
drv_scale (jgm, 0);

beginfig(400)
dcl 0 ( ) ("", 0) "\gamma";
dcl 1 (0) ("", 1) "\delta";
draw drv_tree;
endfig;

drv_font_size "\normalsize";
drv_scale (jgm, 1);

%
% drv_bbox
%

beginfig(410)
dcl 0 (1, 5) ("", _)  "a";
  dcl 1 (2, 3, 4) ("", _)   "b";
    dcl 2 () ("", _)          "c";
    dcl 3 () ("", _)          "d";
    dcl 4 () ("", _)          "e";
  dcl 5 () ("", _)          "f";
fill drv_bbox 1 withcolor (255, 230, 205)/255; % rgb color
draw drv_tree;
endfig;

beginfig(411)
dcl 0 (bxd 1, 5) ("", _)  "a";
  dcl 1 (2, 3, 4) ("", _)   "b";
    dcl 2 () ("", _)          "c";
    dcl 3 () ("", _)          "d";
    dcl 4 () ("", _)          "e";
  dcl 5 () ("", _)          "f";
fill drv_bbox 1 withcolor (255, 230, 205)/255; % rgb color
draw drv_tree;
endfig;

%
% drv_axis
%

beginfig(420)
Dcl 0 ( ) ("",   "", "", _) "a";
Dcl 1 (0) ("", "{}", "", _) "b";
draw drv_tree;
drv_axis (iln, 0);
endfig;

beginfig(421)
Dcl 0 ( ) ("",   "", "", _) "a";
Dcl 1 (0) ("", "{}", "", _) "b";
draw drv_tree;
drv_axis (jdg, 1);
endfig;

beginfig(422)
Dcl 0 ( ) ("",   "", "", _) "a";
Dcl 1 (0) ("", "{}", "", _) "b";
draw drv_tree;
drv_axis (dlm, 1);
endfig;

%
% NFR
%

beginfig(430)
jgm 0 "a";
NFR 0 () ("1", "2", "3", "4", _, _, _);
draw drv_tree;
endfig;

%
% DCL
%

beginfig(440)
DCL 0 ()     ("1", "", "", "", _, _, _) "a";
DCL 1 ()     ("", "2", "", "", _, _, _) "b";
DCL 2 (0, 1) ("", "", "3", "", _, _, _) "c";
DCL 3 (2)    ("", "", "", "4", _, _, _) "d";
draw drv_tree;
endfig;

%
% MVD
%

beginfig(450)
jgm 0 "a";
  jgm 1 "b";
    jgm 2 "cccccc";
    jgm 3 "ddd";
  jgm 4 "eeeeeeeee";
nfr 0 (1, MVD 4 (5, "", "F", r, 4)) ("", _);
  nfr 1 (MVD 2 (2, "G", "", l, 3), 3) ("", _);
    nfr 2 () ("", _);
    nfr 3 () ("", _);
  nfr 4 () ("", _);
draw drv_tree;
endfig;

beginfig(460)
jgm 0 "\textsc{Size matters -- Part 2}";
jgm 1 "\text{Here is an even longer judgment"&
      " that I don't want to shorten either.}";
jgm 2 "\text{This time I'm pretty sure that the"&
      " derivation tree won't fit on the page.}";
jgm 3 "\text{It does! Amazing.}";
nfr 0 () ("", 0);
nfr 1 (0) ("", 1);
nfr 2 () ("", 1);
nfr 3 (MVD 1 (2, "", "", l, 3), 2) ("", 1);
draw drv_tree;
endfig;

%
% Components and central points
%

drv_font_size "\large";

beginfig(470) % components
DCL 5 () ("", "", "", "", _, _, 0) "A";
DCL 6 (5) ("(1)", "(2)", "(3)", "(4)", _, _, 1) "B", "C";
DCL 7 (MVD 6 (2, "(5)", "(6)", _, 6)) ("", "", "", "", _, _, 1) "D";
drv_freeze; % usually called by drv_tree
draw jdg[5];                          % judgment                   A
draw sbj[6][0] withcolor (3, 0, 0)/3; % sub-judgment               B
draw sbj[6][1] withcolor (2, 0, 0)/3; % sub-judgment               C
draw iln[6]    withcolor (0, 4, 0)/4; % inference line
draw l_ilb[6]  withcolor (0, 3, 0)/4; % left inference label      (1)
draw r_ilb[6]  withcolor (0, 2, 0)/4; % right inference label     (2)
draw l_dlm[6]  withcolor (0, 0, 3)/3; % left delimiter
draw l_dlb[6]  withcolor (0, 0, 2)/3; % left delimiter label      (3)
draw r_dlm[6]  withcolor (3, 0, 3)/3; % right delimiter
draw r_dlb[6]  withcolor (2, 0, 2)/3; % right delimiter label     (4)
draw phm[6]    withcolor (0, 3, 3)/4; % phantom steps path
draw l_plb[6]  withcolor (0, 2, 2)/4; % left phantom steps label  (5)
draw r_plb[6]  withcolor (0, 1, 1)/4; % right phantom steps label (6)
draw jdg[7];                          % judgment                   D
draw iln[7];                          % inference line
endfig;

beginfig(471) % central points
DCL 5 () ("", "", "", "", _, _, 0) "A";
DCL 6 (5) ("(1)", "(2)", "(3)", "(4)", _, _, 1) "B", "C";
DCL 7 (MVD 6 (2, "(5)", "(6)", _, 6)) ("", "", "", "", _, _, 1) "D";
pickup pencircle scaled 2.5;
draw drv_tree    withcolor (3, 3, 3)/4;
draw jdg[5].c;                          % judgment                   A
draw sbj[6][0].c withcolor (3, 0, 0)/3; % sub-judgment               B
draw sbj[6][1].c withcolor (2, 0, 0)/3; % sub-judgment               C
draw iln[6].c    withcolor (0, 4, 0)/4; % inference line
draw l_ilb[6].c  withcolor (0, 3, 0)/4; % left inference label      (1)
draw r_ilb[6].c  withcolor (0, 2, 0)/4; % right inference label     (2)
draw l_dlm[6].c  withcolor (0, 0, 3)/3; % left delimiter
draw l_dlb[6].c  withcolor (0, 0, 2)/3; % left delimiter label      (3)
draw r_dlm[6].c  withcolor (3, 0, 3)/3; % right delimiter
draw r_dlb[6].c  withcolor (2, 0, 2)/3; % right delimiter label     (4)
draw phm[6].c    withcolor (0, 3, 3)/4; % phantom steps path
draw l_plb[6].c  withcolor (0, 2, 2)/4; % left phantom steps label  (5)
draw r_plb[6].c  withcolor (0, 1, 1)/4; % right phantom steps label (6)
draw jdg[7].c;                          % judgment                   D
draw iln[7].c;                          % inference line
endfig;

drv_font_size "\normalsize";

%
% Dimensions
%

vardef aux_drawdblarrow expr pth=
save ah;
path ah;
ah:=ahlength*(dir (180-ahangle/2))--(0, 0)--
    ahlength*(dir (180+ahangle/2));
if direction infinity of pth<>(0, 0):
  draw pth;
  draw ah rotated (angle (direction infinity of pth))
    shifted (point infinity of pth);
  draw ah rotated (180+angle (direction -infinity of pth))
    shifted (point -infinity of pth);
fi
enddef;

drv_box
  "jdg_ma[]",
  "jdg_cp[]",
  "jdg_hg[]",
  "jdg_dp[]",
  "iln_cp[]",
  "n_hg",
  "d_dp";
for idx=0, 1:
  box_init jdg_ma[idx](mth_tex[1]"\textit{math axis}");
  box_init jdg_cp[idx](mth_tex[2]"\texttt{jdg["&decimal idx&"].c}");
  box_init jdg_hg[idx](mth_tex[1]"\texttt{jdg["&decimal idx&"].hg}");
  box_init jdg_dp[idx](mth_tex[1]"\texttt{jdg["&decimal idx&"].dp}");
  box_init iln_cp[idx](mth_tex[2]"\texttt{iln["&decimal idx&"].c}");
endfor
box_init n_hg(mth_tex[1]"\texttt{num\_hg}");
box_init d_dp(mth_tex[1]"\texttt{den\_dp}");

drv_font_size "\fontsize{100}{0}";

beginfig(480)
dcl 0 ()  ("", 0) "\gamma";
dcl 1 (0) ("", 1) "\delta";
draw drv_tree withcolor .75*background;
interim linejoin:=beveled;
interim ahangle:=60;
lref[0]:=xpart iln[1].l-20pt;
lref[1]:=xpart iln[1].l-10pt;
lref[2]:=xpart iln[1].l-25pt;
rref[0]:=xpart iln[1].r+20pt;
rref[1]:=xpart iln[1].r+10pt;
rref[2]:=xpart iln[1].r+25pt;
for idx=0, 1:
  forsuffixes sfx=l, c ,r:
    draw jdg[idx].sfx withpen pencircle scaled 2.5;
  endfor
  lx[idx]:=xpart jdg[idx].l;
  rx[idx]:=xpart jdg[idx].r;
  ty[idx]:=ypart jdg[idx].c+jdg[idx].hg;
  cy[idx]:=ypart jdg[idx].c;
  by[idx]:=ypart jdg[idx].c+jdg[idx].dp;
  for pth=
    (lx[idx], ty[idx])--(rref[0], ty[idx]),
    (lref[0], cy[idx])--(rref[0], cy[idx]),
    (lx[idx], by[idx])--(rref[0], by[idx]),
    (lx[idx], ty[idx])--(lx[idx], by[idx]),
    (rx[idx], ty[idx])--(rx[idx], by[idx]):
      draw pth dashed evenly scaled .5;
  endfor
  for pth=
    (rref[1], ty[idx])--(rref[1], cy[idx]),
    (rref[1], cy[idx])--(rref[1], by[idx]):
      aux_drawdblarrow pth;
  endfor
  jdg_ma[idx].r=(lref[2], cy[idx]);
  jdg_cp[idx].c=jdg[idx].c shifted (0, -7.5pt);
  jdg_hg[idx].l=(rref[2], (ty[idx]+cy[idx])/2);
  jdg_dp[idx].l=(rref[2], (cy[idx]+by[idx])/2);
  forsuffixes sfx=jdg_ma, jdg_cp, jdg_hg, jdg_dp:
    box_freeze sfx[idx];
    draw sfx[idx];
  endfor
endfor
draw iln[1].c withpen pencircle scaled 2.5;
iln_cp[1].c=iln[1].c shifted (0, -7.5pt);
box_freeze iln_cp[1];
draw iln_cp[1];
iln_y:=ypart iln[1].c;
draw (lref[0], iln_y)--(rref[0], iln_y) dashed evenly scaled .5;
aux_drawdblarrow (lref[1], cy[0])--(lref[1], iln_y);
aux_drawdblarrow (lref[1], iln_y)--(lref[1], cy[1]);
n_hg.r=(lref[2], (cy[0]+iln_y)/2);
d_dp.r=(lref[2], (iln_y+cy[1])/2);
forsuffixes sfx=n_hg, d_dp:
  box_freeze sfx;
  draw sfx;
endfor
endfig;

drv_font_size "\normalsize";

%
% drv_styled
%

beginfig(490)
jgm 4 "A", "\vdash", "A";
jgm 5 "B", "_2", "\vdash", "B", "_3";
jgm 6 "A", ",", "A", "\multimap", "B", "_1", "\vdash", "B", "_4";
jgm 7 "A", "\multimap", "B", "_0", "\vdash", "A", "\multimap", "B", "_5";
nfr 4 () ("1", _);
nfr 5 () ("1", _);
nfr 6 (4, 5) ("\multimap_{L}", _);
nfr 7 (6) ("\multimap_{R}", _);
drv_freeze;                               % B
draw (sbj[7][2].c shifted (0, -num_hg) .. % 0
      sbj[7][2].c {up} ..                 % 0
      sbj[6][4].c ..                      % 1
      sbj[5][0].c .. tension 1.2 ..       % 2
      sbj[5][3].c ..                      % 3
      sbj[6][7].c ..                      % 4
      sbj[7][7].c {down} ..               % 5
      sbj[7][7].c shifted (0, -num_hg))   % 5
drv_styled 2 withcolor (1, 0, 0);
draw drv_tree;
endfig;

%
% ``User specified'' junction style
%

beginfig(500)
jgm 0 "{\cdot}";
jgm 1 "{\cdot}";
jgm 2 "\text{You may check that"&
  if drv_radialmode_on:
    " the angle between the two dots above is 45$^\circ$.}"
  else:
    " the distance between the two dots above is 5 cm.}"
  fi;
NFR 0 () ("", "", "", "", _, _, _);
NFR 1 () ("", "", "", "", _, _, _);
NFR 2 (0, 1) ("", "", "", "", 3, _, _); % caution: 3
if drv_radialmode_on:
  jdg[1].lng=jdg[0].rng-45;
else:
  xpart jdg[1].c=xpart jdg[0].c+5cm;
fi
draw drv_tree;
endfig;

%
% ``User specified'' alignment style
%

beginfig(510)                                          % "\vdash":
jgm 0 "B, A, \Gamma", "\vdash", "C";                   % sbj[0][1]
jgm 1 "A, \Gamma", "\vdash", "B\multimap C";           % sbj[1][1]
jgm 2 "\Gamma", "\vdash", "A\multimap(B\multimap C)";  % sbj[2][1]
NFR_opt 0 () () () (_, _, 0);
NFR_opt 1 (0) ("\multimap_R") () (_, u, 1);            % caution: u
NFR_opt 2 (1) ("\multimap_R") () (_, u, 1);            % caution: u
if drv_radialmode_on:
  sbj[1][1].cng=sbj[0][1].cng;
  sbj[2][1].cng=sbj[1][1].cng;
else:
  xpart sbj[1][1].c=xpart sbj[0][1].c;
  xpart sbj[2][1].c=xpart sbj[1][1].c;
fi
draw drv_tree;
endfig;

beginfig(511)                                          % "\vdash":
jgm 0 "B, A, \Gamma", "\vdash", "C";                   % sbj[0][1]
jgm 1 "A, \Gamma", "\vdash", "B\multimap C";           % sbj[1][1]
jgm 2 "\Gamma", "\vdash", "A\multimap(B\multimap C)";  % sbj[2][1]
NFR_opt 0 () () () (_, _, 0);
NFR_opt 1 (0) ("\multimap_R") () (_, u, 1);            % caution: u
NFR_opt 2 (1) ("\multimap_R") () (_, u, 1);            % caution: u
if drv_radialmode_on:
  sbj[1][1].lng=sbj[0][1].rng;
  sbj[2][1].lng=sbj[1][1].rng;
else:
  xpart sbj[1][1].l=xpart sbj[0][1].r;
  xpart sbj[2][1].l=xpart sbj[1][1].r;
fi
draw drv_tree;
endfig;

%
% debugging
%

beginfig(520)
jgm 0 "A\vdash B";
jgm 1 "B\vdash C";
jgm 2 "A\vdash C";
jgm 3 "C\vdash D";
jgm 4 "A\vdash D";
nfr 0 () ("f", _);
nfr 1 () ("g", _);
nfr 2 (0, 1) ("\circ", _);
nfr 3 () ("h", _);
nfr 4 (2, 3) ("\circ", _);
draw drv_tree;
endfig;

%
% derivation forests
%

beginfig(530)
% first tree
dcl 10 () ("", _) "a";
% second tree
dcl 20 (21, 22) ("", _) "d";
  dcl 21 () ("", _) "b";
  dcl 22 () ("", _) "c";
draw drv_tree;
endfig;

beginfig(531)
% first tree
dcl 10 () ("", _) "a";
% second tree
dcl 20 (21, 22) ("", _) "d";
  dcl 21 () ("", _) "b";
  dcl 22 () ("", _) "c";
% relative positionning
if drv_radialmode_on:
  jdg[10].rad=jdg[22].rad;
else:
  ypart jdg[10].c=ypart jdg[22].c;
fi
draw drv_tree;
endfig;

beginfig(532)
% first tree
dcl 10 () ("", _) "a";
% second tree
dcl 20 (21, 22) ("", _) "d";
  dcl 21 () ("", _) "b";
  dcl 22 () ("", _) "c";
% relative positionning
if drv_radialmode_on:
  iln[10].rng=iln[20].lng;
else:
  xpart iln[10].r=xpart iln[20].l;
fi
draw drv_tree;
endfig;

%
% drv_root
%

beginfig(533)
% first tree
dcl 10 () ("", _) "a";
% second tree
dcl 20 (21, 22) ("", _) "d";
  dcl 21 () ("", _) "b";
  dcl 22 () ("", _) "c";
drv_root (20, t); % root at top!
draw drv_tree;
endfig;

drv_left_delimiter "\downarrow";
drv_right_delimiter "\uparrow";

beginfig(540)
% first tree
jgm 10 "A\vdash D";
Nfr 10 (11, 14) ("\circ", "h\circ (g\circ f)", "", 1);
  dcl 11 (12, 13) ("\circ", 1) "A\vdash C";
    dcl 12 () ("f", 2) "A\vdash B";
    dcl 13 () ("g", 3) "B\vdash C";
  dcl 14 () ("h", 4) "C\vdash D";
% second tree
jgm 20 "\phantom{A\vdash D}";
Nfr 20 (21, 22) ("\circ", "", "(h\circ g)\circ f", 1);
  dcl 21 () ("f", 2) "A\vdash B";
  dcl 22 (23, 24) ("\circ", 1) "B\vdash D";
    dcl 23 () ("g", 3) "B\vdash C";
    dcl 24 () ("h", 4) "C\vdash D";
drv_root (20, t); % root at top!
% overlapping
jdg[10].c=jdg[20].c;
draw drv_tree;
endfig;

drv_left_delimiter "\lbrace";
drv_right_delimiter "\rbrace";

%
% drv_radial_mode
%

drv_radial_mode on;

beginfig(550)
dcl 0 () ("f", 2) "A\vdash B";
dcl 1 () ("g", 3) "B\vdash C";
dcl 2 () ("h", 4) "C\vdash D";
dcl 3 (1, 2) ("\circ", 1) "B\vdash D";
dcl 4 (0, 3) ("\circ", 1) "A\vdash D";
draw drv_tree;
endfig;

drv_roots_position t;

beginfig(551)
dcl 0 () ("f", 2) "A\vdash B";
dcl 1 () ("g", 3) "B\vdash C";
dcl 2 () ("h", 4) "C\vdash D";
dcl 3 (1, 2) ("\circ", 1) "B\vdash D";
dcl 4 (0, 3) ("\circ", 1) "A\vdash D";
draw drv_tree;
endfig;

drv_roots_position b;
drv_radial_mode off;

beginfig(552)
dcl 0 () ("f", 2) "A\vdash B";
dcl 1 () ("g", 3) "B\vdash C";
dcl 2 () ("h", 4) "C\vdash D";
dcl 3 (1, 2) ("\circ", 1) "B\vdash D";
dcl 4 (0, 3) ("\circ", 1) "A\vdash D";
draw drv_tree;
endfig;

drv_roots_position t;

beginfig(553)
dcl 0 () ("f", 2) "A\vdash B";
dcl 1 () ("g", 3) "B\vdash C";
dcl 2 () ("h", 4) "C\vdash D";
dcl 3 (1, 2) ("\circ", 1) "B\vdash D";
dcl 4 (0, 3) ("\circ", 1) "A\vdash D";
draw drv_tree;
endfig;

drv_roots_position b;

%
% drv_scale (crv, );
%

drv_radial_mode on;
drv_scale (crv, 0.5);

beginfig(560)
dcl 0 () ("f", 2) "A\vdash B";
dcl 1 () ("g", 3) "B\vdash C";
dcl 2 (0, 1) ("\circ", 1) "A\vdash C";
draw drv_tree;
endfig;

drv_scale (crv, 1);

beginfig(561)
dcl 0 () ("f", 2) "A\vdash B";
dcl 1 () ("g", 3) "B\vdash C";
dcl 2 (0, 1) ("\circ", 1) "A\vdash C";
draw drv_tree;
endfig;

drv_scale (crv, 2);

beginfig(562)
dcl 0 () ("f", 2) "A\vdash B";
dcl 1 () ("g", 3) "B\vdash C";
dcl 2 (0, 1) ("\circ", 1) "A\vdash C";
draw drv_tree;
endfig;

drv_scale (crv, 4);

beginfig(563)
dcl 0 () ("f", 2) "A\vdash B";
dcl 1 () ("g", 3) "B\vdash C";
dcl 2 (0, 1) ("\circ", 1) "A\vdash C";
draw drv_tree;
endfig;

drv_scale (crv, 1);
drv_radial_mode off;

%
% drv_azimuth
%

drv_radial_mode on;
drv_azimuth 120;

beginfig(570)
dcl 0 () ("f", 2) "A\vdash B";
dcl 1 () ("g", 3) "B\vdash C";
dcl 2 (0, 1) ("\circ", 1) "A\vdash C";
draw drv_tree;
endfig;

drv_azimuth 90;

beginfig(571)
dcl 0 () ("f", 2) "A\vdash B";
dcl 1 () ("g", 3) "B\vdash C";
dcl 2 (0, 1) ("\circ", 1) "A\vdash C";
draw drv_tree;
endfig;

drv_azimuth 60;

beginfig(572)
dcl 0 () ("f", 2) "A\vdash B";
dcl 1 () ("g", 3) "B\vdash C";
dcl 2 (0, 1) ("\circ", 1) "A\vdash C";
draw drv_tree;
endfig;

drv_azimuth 90;
drv_radial_mode off;

%
% Size matters - Part 3
%

drv_radial_mode on;

beginfig(580)
dcl 0 () ("", 0) "\textsc{Size matters -- Part 3}";
dcl 1 (0) ("", 1) "\text{How could this very long judgment fit on the page"&
                  " without being desperately bent into a circle arc ?}";
draw drv_tree;
endfig;

drv_radial_mode off;

%
% User specified junction and alignment styles
%

drv_radial_mode on;

beginfig(590)                                          % "\vdash":
jgm 0 "B, A, \Gamma", "\vdash", "C";                   % sbj[0][1]
jgm 1 "A, \Gamma", "\vdash", "B\multimap C";           % sbj[1][1]
jgm 2 "\Gamma", "\vdash", "A\multimap(B\multimap C)";  % sbj[2][1]
NFR_opt 0 () () () (_, _, 0);
NFR_opt 1 (0) ("\multimap_R") () (_, u, 1);            % caution: u
NFR_opt 2 (1) ("\multimap_R") () (_, u, 1);            % caution: u
sbj[1][1].cng=sbj[0][1].cng;
sbj[2][1].cng=sbj[1][1].cng;
draw drv_tree;
endfig;

beginfig(591)                                          % "\vdash":
jgm 0 "B, A, \Gamma", "\vdash", "C";                   % sbj[0][1]
jgm 1 "A, \Gamma", "\vdash", "B\multimap C";           % sbj[1][1]
jgm 2 "\Gamma", "\vdash", "A\multimap(B\multimap C)";  % sbj[2][1]
NFR_opt 0 () () () (_, _, 0);
NFR_opt 1 (0) ("\multimap_R") () (_, u, 1);            % caution: u
NFR_opt 2 (1) ("\multimap_R") () (_, u, 1);            % caution: u
sbj[1][1].lng=sbj[0][1].rng;
sbj[2][1].lng=sbj[1][1].rng;
draw drv_tree;
endfig;

drv_radial_mode off;

%
% example - triple unit
%

drv_verbatimtex
  "\renewcommand{\limp}{\thinspace{\multimapinv}\thinspace}"&
  "\renewcommand{\rimp}{\thinspace{\multimap}\thinspace}";

string id, Llimp, Rlimp, Lrimp, Rrimp;
id:="\text{id}";
Llimp:="\overset{\multimapinv}{\eta}";
Rlimp:="\overset{\multimapinv}{\varepsilon}";
Lrimp:="\overset{\multimap}{\eta}";
Rrimp:="\overset{\multimap}{\varepsilon}";

beginfig(600)
jgm 10 "1\limp(a\rimp 1)\otimes\bigl(1\limp(a\rimp 1)\bigr)\rimp 1\vdash 1";
  jgm 11 "1\limp(a\rimp 1)\vdash 1\limp(a\rimp 1)";
    jgm 12 "1\limp(a\rimp 1)\otimes a\rimp 1\vdash 1";
      jgm 13 "a\rimp 1\vdash a\rimp 1";
        jgm 14 "a\otimes a\rimp 1\vdash 1";
          jgm 15 "a\vdash a";
          jgm 16 "1\vdash 1";
      jgm 17 "1\vdash 1";
  jgm 18 "1\vdash 1";
nfr 10 (11, 18) (Lrimp, 1);
  nfr 11 (12) (Rlimp, 1);
    nfr 12 (13, 17) (Llimp, 1);
      nfr 13 (14) (Rrimp, 1);
        nfr 14 (15, 16) (Lrimp, 1);
          nfr 15 () (id, 1);
          nfr 16 () (id, 1);
      nfr 17 () (id, 1);
  nfr 18 () (id, 1);
draw drv_tree;
endfig;

beginfig(601)
jgm 20 "1\limp(a\rimp 1)\otimes\bigl(1\limp(a\rimp 1)\bigr)\rimp 1\vdash 1";
  jgm 21 "\bigl(1\limp(a\rimp 1)\bigr)\rimp 1\vdash a\rimp 1";
    jgm 22 "a\otimes\bigl(1\limp(a\rimp 1)\bigr)\rimp 1\vdash 1";
      jgm 23 "a\vdash 1\limp(a\rimp 1)";
        jgm 24 "a\otimes a\rimp 1\vdash 1";
          jgm 25 "a\vdash a";
          jgm 26 "1\vdash 1";
      jgm 27 "1\vdash 1";
  jgm 28 "1\vdash 1";
nfr 20 (21, 28) (Llimp, 1);
  nfr 21 (22) (Rrimp, 1);
    nfr 22 (23, 27) (Lrimp, 1);
      nfr 23 (24) (Rlimp, 1);
        nfr 24 (25, 26) (Lrimp, 1);
          nfr 25 () (id, 1);
          nfr 26 () (id, 1);
      nfr 27 () (id, 1);
  nfr 28 () (id, 1);
draw drv_tree;
endfig;

%
% example - Restall
%

string Lor, Ror[], Land[], Rand;
Lor:="\vee L";
Ror1:="\vee R_1";
Ror2:="\vee R_2";
Land1:="\wedge L_1";
Land2:="\wedge L_2";
Rand:="\wedge R";

beginfig(610)
jgm 0 "p\wedge(q\vee(r_1\wedge r_2))\vdash p\wedge(q\vee(r_1\wedge r_2))";
  jgm 1 "p\wedge(q\vee(r_1\wedge r_2))\vdash p";
    jgm 2 "p\vdash p";
  jgm 3 "p\wedge(q\vee(r_1\wedge r_2))\vdash q\vee(r_1\wedge r_2)";
    jgm 4 "q\vee(r_1\wedge r_2)\vdash q\vee(r_1\wedge r_2)";
      jgm 5 "q\vdash q\vee(r_1\wedge r_2)";
        jgm 6 "q\vdash q";
      jgm 7 "r_1\wedge r_2\vdash q\vee(r_1\wedge r_2)";
        jgm 8 "r_1\wedge r_2\vdash r_1\wedge r_2";
          jgm 9 "r_1\wedge r_2\vdash r_1";
            jgm 10 "r_1\vdash r_1";
          jgm 11 "r_1\wedge r_2\vdash r_2";
            jgm 12 "r_2\vdash r_2";
nfr 0 (1, 3) (Rand, 1);
  nfr 1 (2) (Land1, 1);
    nfr 2 () ("", 0);
  nfr 3 (bxd 4) (Land2, 1);		% fig. 611 -> nfr 3 (mvd 4 (1, 3))
    nfr 4 (5, 7) (Lor, 1);		% fig. 611 -> Nfr
      nfr 5 (6) (Ror1, 1);
        nfr 6 () ("", 0);
      nfr 7 (bxd 8) (Ror2, 1);		% fig. 611 -> nfr 7 (mvd 8 (1, 3))
        nfr 8 (9, 11) (Rand, 1);	% fig. 611 -> Nfr
          nfr 9 (10) (Land1, 1);
            nfr 10 () ("", 0);
          nfr 11 (12) (Land2, 1);
            nfr 12 () ("", 0);
fill drv_bbox 4 withcolor (255, 230, 205)/255;
fill drv_bbox 8 withcolor (236, 211, 185)/255;
draw drv_tree;
endfig;

beginfig(611)
jgm 0 "p\wedge(q\vee(r_1\wedge r_2))\vdash p\wedge(q\vee(r_1\wedge r_2))";
  jgm 1 "p\wedge(q\vee(r_1\wedge r_2))\vdash p";
    jgm 2 "p\vdash p";
  jgm 3 "p\wedge(q\vee(r_1\wedge r_2))\vdash q\vee(r_1\wedge r_2)";
    jgm 4 "q\vee(r_1\wedge r_2)\vdash q\vee(r_1\wedge r_2)";
      jgm 5 "q\vdash q\vee(r_1\wedge r_2)";
        jgm 6 "q\vdash q";
      jgm 7 "r_1\wedge r_2\vdash q\vee(r_1\wedge r_2)";
        jgm 8 "r_1\wedge r_2\vdash r_1\wedge r_2";
          jgm 9 "r_1\wedge r_2\vdash r_1";
            jgm 10 "r_1\vdash r_1";
          jgm 11 "r_1\wedge r_2\vdash r_2";
            jgm 12 "r_2\vdash r_2";
nfr 0 (1, 3) (Rand, 1);
  nfr 1 (2) (Land1, 1);
    nfr 2 () ("", 0);
  nfr 3 (mvd 4 (1, 3)) (Land2, 1);
    Nfr 4 (5, 7) (Lor, "\text{Id\,}_{q\vee(r_1\wedge r_2)}", "", 1);
      nfr 5 (6) (Ror1, 1);
        nfr 6 () ("", 0);
      nfr 7 (mvd 8 (1, 3)) (Ror2, 1);
        Nfr 8 (9, 11) (Rand, "\text{Id\,}_{r_1\wedge r_2}", "", 1);
          nfr 9 (10) (Land1, 1);
            nfr 10 () ("", 0);
          nfr 11 (12) (Land2, 1);
            nfr 12 () ("", 0);
fill drv_bbox 4 withcolor (255, 230, 205)/255;
fill drv_bbox 8 withcolor (236, 211, 185)/255;
draw drv_tree;
endfig;

%
% example - Roegel
%

drv_junction_style 0;
drv_scale (jdg, 0.9);
drv_scale (prm, 0.9);

vardef mix[]=
  "\textit{mix}("&(decimal @)&")"
enddef;

beginfig(620)
jgm 0 "\Gamma_A, \Gamma'\vdash\Delta, \Delta'_A";
  jgm 1 "\Gamma_A, \Gamma', \Gamma_A, \Gamma', \Gamma_A, \Gamma'\vdash"&
    "\Delta, \Delta'_A, \Delta, \Delta'_A, \Delta, \Delta'_A";
    jgm 2 "\Gamma_A, \Gamma', \Gamma_A\vdash"&    % a space character is
      " C, \Delta, \Delta'_A, \Delta, \Delta'_A"; % needed before 'C'!
      jgm 3 "\Gamma_A, \Gamma'\vdash B, C, \Delta, \Delta'_A";
        jgm 4 "\Gamma, B\vee C\vdash \Delta";
          jgm 5 "\Gamma, B\vdash \Delta";
            jgm 6 "\Pi_1";
          jgm 7 "\Gamma, C\vdash \Delta";
            jgm 8 "\Pi_2";
        jgm 9 "\Gamma'\vdash B, C, \Delta'";
          jgm 10 "\Pi_3";
      jgm 11 "\Gamma_A, \Gamma', B\vdash\Delta, \Delta'_A";
        jgm 12 "\Gamma, B\vdash \Delta";
          jgm 13 "\Pi_1";
        jgm 14 "\Gamma'\vdash B\vee C, \Delta'";
          jgm 15 "\Gamma'\vdash B, C, \Delta'";
            jgm 16 "\Pi_3";
    jgm 17 "\Gamma_A, \Gamma', C\vdash\Delta, \Delta'_A";
      jgm 18 "\Gamma, C\vdash \Delta";
        jgm 19 "\Pi_2";
      jgm 20 "\Gamma'\vdash B\vee C, \Delta'";
        jgm 21 "\Gamma'\vdash B, C, \Delta'";
          jgm 22 "\Pi_3";
nfr 0 (1) ("\text{contr$_g$, contr$_d$}", 1);
  nfr 1 (2, 17) (mix5, 1);
    nfr 2 (mvd 3 (5, 6), 11) (mix3 , 1);
      nfr 3 (4, 9) (mix1, 1);
        nfr 4 (5, 7) ("\vee_g", 1);
          nfr 5 (mvd 6 (1, 6)) ("", 0);
            nfr 6 () ("", 0);
          nfr 7 (mvd 8 (1, 6)) ("", 0);
            nfr 8 () ("", 0);
        nfr 9 (mvd 10 (1, 6)) ("", 0);
          nfr 10 () ("", 0);
      nfr 11 (12, 14) (mix2, 1);
        nfr 12 (mvd 13 (1, 6)) ("", 0);
          nfr 13 () ("", 0);
        nfr 14 (15) ("\vee_d", 1);
          nfr 15 (mvd 16 (1, 6)) ("", 0);
            nfr 16 () ("", 0);
    nfr 17 (18, 20) (mix4, 1);
      nfr 18 (mvd 19 (1, 6)) ("", 0);
        nfr 19 () ("", 0);
      nfr 20 (21) ("\vee_d", 1);
        nfr 21 (mvd 22 (1, 6)) ("", 0);
          nfr 22 () ("", 0);
draw drv_tree;
endfig;

drv_junction_style 1;
drv_scale (jdg, 1);
drv_scale (prm, 1);

%
% example - Lutz
%

drv_math_style (ilb, "\textstyle");
drv_labels_position (ilb, l);
drv_scale (clr, 1.5);
drv_axis_reference jdg;

string neg, par, otm;
neg:="^\bot";
par:="\mathbin{\bindnasrepma}"; % from stmaryrd.sty
otm:="\otimes";

beginfig(630)
    %  0    1    2    3    4    5    6    7    8    9
jgm 0 "a", neg, par, "a";
jgm 1 "a", neg, par, "(", "a", otm, "(", "a", par, "a", neg, ")", ")";
jgm 2 "a", neg, par, "(", "a", otm, "a", ")", par, "a", neg;
jgm 3 "a", neg, par, "(", "a", otm, "a", ")", par, "(",
      "(", "a", par, "a", neg, ")", otm, "a", neg, ")";
jgm 4 "a", neg, par, "(", "a", otm, "a", ")", par, "a",
      par, "(", "a", neg, otm, "a", neg, ")";
nfr 0 ( ) ("\textsf{id}", 1);
nfr 1 (0) ("\textsf{id}", 1);
nfr 2 (1) ("\textsf{s}",  1);
nfr 3 (2) ("\textsf{id}", 1);
nfr 4 (3) ("\textsf{s}",  1);
drv_freeze;
drawoptions (withpen drv_pen scaled 2 withcolor (255, 94, 94)/255);
draw
  sbj[4][12].c shifted (0, -num_hg) .. {up}
  sbj[4][12].c ..
  sbj[3][13].c .. tension .9 ..
  sbj[3][11].c ..
  sbj[4][9].c {down} ..
  sbj[4][9].c shifted (0, -num_hg);
draw
  sbj[4][15].c shifted (0, -num_hg) .. {up}
  sbj[4][15].c ..
  sbj[3][17].c .. tension 1.2 ..
  sbj[2][9].c ..
  sbj[1][9].c {up} .. tension 1.5.. {down}
  sbj[1][7].c ..
  sbj[2][6].c .. tension 1.2 ..
  sbj[3][6].c ..
  sbj[4][6].c {down} ..
  sbj[4][6].c shifted (0, -num_hg);
draw
  sbj[4][0].c shifted (0, -num_hg) .. {up}
  sbj[4][0].c ..
  sbj[3][0].c .. tension 1.2 ..
  sbj[2][0].c ..
  sbj[1][0].c .. tension 1.2 ..
  sbj[0][0].c {up} .. tension 1.8.. {down}
  sbj[0][3].c .. tension 1.2 ..
  sbj[1][4].c ..
  sbj[2][4].c .. tension 1.2 ..
  sbj[3][4].c ..
  sbj[4][4].c {down} ..
  sbj[4][4].c shifted (0, -num_hg);
drawoptions ();
draw drv_tree; % fig. 631 -> draw jdg[4];
endfig;

beginfig(631)
    %  0    1    2    3    4    5    6    7    8    9
jgm 0 "a", neg, par, "a";
jgm 1 "a", neg, par, "(", "a", otm, "(", "a", par, "a", neg, ")", ")";
jgm 2 "a", neg, par, "(", "a", otm, "a", ")", par, "a", neg;
jgm 3 "a", neg, par, "(", "a", otm, "a", ")", par, "(",
      "(", "a", par, "a", neg, ")", otm, "a", neg, ")";
jgm 4 "a", neg, par, "(", "a", otm, "a", ")", par, "a",
      par, "(", "a", neg, otm, "a", neg, ")";
nfr 0 ( ) ("\textsf{id}", 1);
nfr 1 (0) ("\textsf{id}", 1);
nfr 2 (1) ("\textsf{s}",  1);
nfr 3 (2) ("\textsf{id}", 1);
nfr 4 (3) ("\textsf{s}",  1);
drv_freeze;
drawoptions (withpen drv_pen scaled 2 withcolor (255, 94, 94)/255);
draw
  sbj[4][12].c shifted (0, -num_hg) .. {up}
  sbj[4][12].c ..
  sbj[3][13].c .. tension .9 ..
  sbj[3][11].c ..
  sbj[4][9].c {down} ..
  sbj[4][9].c shifted (0, -num_hg);
draw
  sbj[4][15].c shifted (0, -num_hg) .. {up}
  sbj[4][15].c ..
  sbj[3][17].c .. tension 1.2 ..
  sbj[2][9].c ..
  sbj[1][9].c {up} .. tension 1.5.. {down}
  sbj[1][7].c ..
  sbj[2][6].c .. tension 1.2 ..
  sbj[3][6].c ..
  sbj[4][6].c {down} ..
  sbj[4][6].c shifted (0, -num_hg);
draw
  sbj[4][0].c shifted (0, -num_hg) .. {up}
  sbj[4][0].c ..
  sbj[3][0].c .. tension 1.2 ..
  sbj[2][0].c ..
  sbj[1][0].c .. tension 1.2 ..
  sbj[0][0].c {up} .. tension 1.8.. {down}
  sbj[0][3].c .. tension 1.2 ..
  sbj[1][4].c ..
  sbj[2][4].c .. tension 1.2 ..
  sbj[3][4].c ..
  sbj[4][4].c {down} ..
  sbj[4][4].c shifted (0, -num_hg);
drawoptions ();
draw jdg[4]; % fig. 630 -> draw drv_tree;
endfig;

drv_math_style (ilb, "\scriptstyle");
drv_labels_position (ilb, r);
drv_scale (clr, 1);
drv_axis_reference iln;

%
% example - a continued fraction
%

drv_math_style (ilb, "\textstyle");
drv_scale (jdg, 0);
drv_scale (ilb, 0);
drv_labels_position (ilb, l);
drv_roots_position t;
drv_box_mode on;

beginfig(640)
dcl 0 (1) ("1+", 1) "a";
dcl 1 (2) ("2+", 1) "b";
dcl 2 (3) ("3+", 1) "c";
dcl 3 (4) ("4+", 1) "d";
dcl 4 ( ) ("", 0) "\cdots";
draw drv_tree;
endfig;

drv_math_style (ilb, "\scriptstyle");
drv_scale (jdg, 1);
drv_scale (ilb, 1);
drv_labels_position (ilb, r);
drv_roots_position b;
drv_box_mode off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 3. Postamble
%

end
