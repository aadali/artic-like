% micromu.mp
% L. Nobre G., Ryan Van Wagoner, Jens-Uwe Morawski
% 2006

verbatimtex
%&latex
\documentclass{article}
\usepackage[T1]{fontenc}
\usepackage[latin1]{inputenc}
\usepackage{textcomp}
\usepackage{pifont}
\newcommand{\upmu}{\Pisymbol{psy}{109}}
\begin{document}
etex

beginfig(1);
  transform T;
  numeric referenceangle, tangentd, shiftlen;
  picture micromu;
  referenceangle = 14;
  shiftlen = 5mm;
  tangentd = sind(referenceangle)/cosd(referenceangle);
  T = identity slanted (-tangentd);
  micromu = thelabel( btex \framebox{\upmu\ \textmu\ $\mu$.} etex, origin );
  draw micromu shifted (0,shiftlen);
  draw micromu transformed T;
endfig;

verbatimtex 
\end{document} 
etex

end.

