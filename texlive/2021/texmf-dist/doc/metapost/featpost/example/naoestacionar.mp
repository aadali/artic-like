% naoestacionar.mp
% L. Nobre G.
% 2008

prologues := 3;

verbatimtex
\documentclass{article}
\usepackage{palatino}
\begin{document}
etex

beginfig(1);
        numeric stampnum, initang, i, stepang, stampoffset, amplfactor;
        picture stamp;
        bboxmargin := 9mm;
        stampnum = 36;
        stampoffset = 0.92mm;
	amplfactor = 20;
        label.urt( btex \begin{tabular}{c}
	  N\~{A}O \\
	  estacionar
	  \end{tabular} etex scaled amplfactor, origin );
        stamp = currentpicture;
        unfill bbox stamp;
        initang = 180/stampnum;
        stepang = 2*initang;
        drawoptions();
        for i=0 upto (stampnum-1):
            draw stamp shifted (stampoffset*dir(initang+i*stepang));
        endfor;
        undraw stamp;
endfig;

verbatimtex
\end{document}
etex

end.
