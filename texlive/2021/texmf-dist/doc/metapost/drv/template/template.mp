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
% `template.tex' (document class, font package, font selection commands ...)
% should be reported here.

\documentclass[twoside,11pt]{article}
\usepackage[T1]{fontenc}
\usepackage{amsmath}
\usepackage{amssymb}

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

\begin{document}
etex;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 2. Figures
%

% Default settings
%
% drv_font_size "\normalsize";
% drv_math_style (drv, "\displaystyle");
% drv_math_style (jdg, "\textstyle");
% drv_math_style (ilb, "\scriptstyle");
% drv_math_style (dlb, "\textstyle");
% drv_math_style (plb, "\textstyle");
% drv_scale (clr, 1);
% drv_scale (prm, 1);
% drv_scale (jdg, 1);
% drv_scale (ilb, 1);
% drv_junction_style 1;
% drv_alignment_style c;
% drv_path_style (iln, 1);
% drv_path_style (phm, 3);
% drv_labels_position (ilb, r);
% drv_labels_position (dlb, l);
% drv_labels_position (plb, l);
% drv_roots_position b;
% drv_axis_reference iln;
% drv_left_delimiter "\lbrace";
% drv_right_delimiter "\rbrace";
% drv_box_mode off;
% drv_fraction_mode on;
% drv_proof_mode off;
% drv_radial_mode off;
% drv_labels_mode (ilb, on);
% drv_labels_mode (dlb, on);
% drv_labels_mode (plb, on);
% drv_scale (crv, 1);
% drv_azimuth 90;

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

string Lor, Ror[], Land[], Rand;
Lor:="\vee L";
Ror1:="\vee R_1";
Ror2:="\vee R_2";
Land1:="\wedge L_1";
Land2:="\wedge L_2";
Rand:="\wedge R";

beginfig(110)
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 3. Postamble
%

end
