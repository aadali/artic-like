input drv;

verbatimtex%&latex

% Any piece of code related to font selection in the preambule of
% `coq-sample.tex' (document class, font package, font selection commands ...)
% should be reported here.

\documentclass[8pt]{beamer}
\usepackage[T1]{fontenc}
\usepackage{amsmath}
\usepackage[varg]{txfonts}
\newcommand{\var}{\text{Var}}
\newcommand{\app}{\text{App}}
\newcommand{\lam}{\text{Lam}}
\begin{document}
etex;

drv_font_size "\small";
drv_labels_position l;

for index=0 upto 6:
  beginfig(index)
  jgm 0	"\vdash", "(P\to(Q\to R))\to(Q\to(P\to R))";
    jgm 1 "P\to(Q\to R)\vdash", "Q\to(P\to R)";
      jgm 2 "P\to(Q\to R), Q\vdash", "P\to R";
        jgm 3 "P\to(Q\to R), Q, P\vdash", "R";
          jgm 4 "P\to(Q\to R), Q, P\vdash", "Q\to R";
            jgm 5 "P\to(Q\to R), Q, P\vdash", "P\to(Q\to R)";
            jgm 6 "P\to(Q\to R), Q, P\vdash", "P";
          jgm 7 "P\to(Q\to R), Q, P\vdash", "Q";
  nfr 0 (1) ("\lam", 1);
    nfr 1 (2) ("\lam", 1);
      nfr 2 (3) ("\lam", 1);
        nfr 3 (4, 7) ("\app", 1);
          nfr 4 (5, 6) ("\app", 1);
            nfr 5 () ("\var", 1);
            nfr 6 () ("\var", 1);
          nfr 7 () ("\var", 1);
  setbounds currentpicture to drv_bbox 0;
    % All the pictures will have the same bounding box
    % as the full derivation tree.
  draw sbj[0][0];
  if index=0:
    draw sbj[0][1] withcolor blue;
  else:
    draw sbj[0][1];
    draw iln[0];
    draw l_ilb[0];
    draw sbj[1][0];
    if index=1:
      draw sbj[1][1] withcolor blue;
    else:
      draw sbj[1][1];
      draw iln[1];
      draw l_ilb[1];
      draw sbj[2][0];
      if index=2:
        draw sbj[2][1] withcolor blue;
      else:
        draw sbj[2][1];
        draw iln[2];
        draw l_ilb[2];
        draw sbj[3][0];
        if index=3:
          draw sbj[3][1] withcolor blue;
        else:
          draw sbj[3][1];
          draw iln[3];
          draw l_ilb[3];
          draw jdg[4];
          draw iln[4];
          draw l_ilb[4];
            draw jdg[5];
            draw iln[5];
            draw l_ilb[5];
            draw sbj[6][0];
          draw sbj[7][0];
          if index=4:
            draw sbj[6][1] withcolor blue;
            draw sbj[7][1] withcolor blue;
          else:
            draw sbj[6][1];
            draw iln[6];
            draw l_ilb[6];
            if index=5:
              draw sbj[7][1] withcolor blue;
            else:
              draw sbj[7][1];
              draw iln[7];
              draw l_ilb[7];
            fi
          fi
        fi
      fi
    fi
  fi
  endfig;
endfor

end
