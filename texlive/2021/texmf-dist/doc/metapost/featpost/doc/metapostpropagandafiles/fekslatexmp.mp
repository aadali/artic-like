% fekslatexmp.mp
verbatimtex \documentclass{beamer}\begin{document} etex
input latexmp;
setupLaTeXMP(
  class="beamer",
  mode=rerun,
  packages="calc,ifthen,mflogo" );
beginfig(1);
  numeric u, angstep, i, m;
  u = 45mm;
  m = 0.1u;
  angstep = 7;
  draw fullcircle scaled 2u; 
  for i=angstep step angstep until 9angstep:
    drawarrow origin--(u*dir(i));
    dotlabel.urt( textext("$\alpha_{" & decimal(i/angstep) & "}$" ), u*dir(i) );
  endfor;
  clip currentpicture to (-m,-m)--(u+m,-m)--(u+m,u+m)--(-m,u+m)--cycle;
endfig;
verbatimtex \end{document} etex
end
