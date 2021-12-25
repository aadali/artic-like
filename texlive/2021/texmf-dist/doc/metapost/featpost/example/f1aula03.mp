% f1aula03.mp
% L. Nobre G.
% 2010 removed package icomma

verbatimtex
\documentclass{beamer}
\usepackage{amsmath,amssymb}
\usepackage[squaren]{SIunits}
\usepackage[latin1]{inputenc}
\usepackage[portuges]{babel}
\newcommand{\veri}{\ensuremath{\hat{\text{i}}}}
\newcommand{\verj}{\ensuremath{\hat{\text{j}}}}
\newcommand{\verk}{\ensuremath{\hat{\text{k}}}}
\newcommand{\veru}{\ensuremath{\hat{\text{u}}}}
\newcommand{\emf}{\ensuremath{\mathcal{E}}}
\newcommand{\vecnorm}[1]{\left|\left|\vec{#1}\right|\right|}
\renewcommand{\arraystretch}{1.45}
\begin{document}
\frame{
etex

ahangle := 36;
ahlength := 3.5mm;
labeloffset := 3.5mm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% {} []
beginfig(1);
  pen basepen, vecpen;
  basepen = pencircle scaled 0.5mm;
  vecpen = pencircle scaled 1.15mm;
  path mainaxe;
  numeric malen, mastep, matic, lscale;
  malen = 10cm;
  mastep = 5mm;
  matic = 2mm;
  lscale = 2.5;
  mainaxe = origin--(malen,0);
  drawarrow mainaxe withpen basepen;
  label.bot( btex $x$ etex scaled lscale, (malen,0) );
  pair refpa, refpb, refpc;
  refpa = (6mastep,0);
  refpb = (12mastep,0);
  refpc = (14mastep,0);
  drawarrow refpb--refpa withpen vecpen withcolor red;
  label.top( btex $\vec{v}$ etex scaled lscale, 0.5[refpa,refpb] );
  drawarrow refpb--refpc withpen vecpen withcolor blue;
  label.top( btex \veri etex scaled lscale, 0.5[refpc,refpb] );
  numeric i;
  for i=mastep step 2mastep until malen-2mastep:
    draw (i,-0.5matic)--(i,0.5matic) withpen basepen;
    draw (i+mastep,-matic)--(i+mastep,matic) withpen basepen;
  endfor;
  label.bot( btex -3 etex scaled lscale, refpa );
  label.bot( btex 0 etex scaled lscale, refpb );
  label.bot( btex 1 etex scaled lscale, refpc );
endfig;

beginfig(2);
  drawarrow mainaxe withpen basepen;
  label.bot( btex $x$ etex scaled lscale, (malen,0) );
  drawarrow refpb--refpa withpen vecpen withcolor red;
  label.top( btex $\Delta x$ etex scaled lscale, 0.5[refpa,refpb] );
  draw (6mastep,-matic)--(6mastep,matic) withpen basepen;
  draw (12mastep,-matic)--(12mastep,matic) withpen basepen;
  label.bot( btex fim etex scaled lscale, refpa );
  label.bot( btex início etex scaled lscale, refpb );
endfig;

beginfig(3);
  pair refpd;
  refpd = (malen,-6mastep);
  drawarrow refpb--refpa withpen vecpen withcolor red;
  drawarrow refpd--refpa withpen basepen;
  drawarrow refpd--refpb withpen basepen;
  label.top( btex $\Delta \vec{x}$ etex scaled lscale, 0.5[refpa,refpb] );
  label.llft( btex $\vec{x}_2$ etex scaled lscale, 0.5[refpd,refpa] );
  label.urt( btex $\vec{x}_1$ etex scaled lscale, 0.5[refpd,refpb] );
endfig;

beginfig(4)
  numeric an;
  an = angle(refpd);
  path angpa, cutp, calfa, cbeta;
  angpa = refpd--origin--(malen,0);
  cutp = origin--(28*mastep*dir(-an));
  calfa = halfcircle scaled (24mastep);
  calfa := calfa cutafter cutp;
  cbeta = halfcircle scaled (28mastep);
  cbeta := cbeta cutafter cutp;
  drawoptions( withpen basepen );
  drawarrow angpa;
  drawarrow fullcircle scaled (12mastep) cutafter angpa;
  drawarrow calfa rotated an;
  drawarrow reverse (cbeta rotated an);
  label.bot( btex $\theta$ etex scaled lscale, (0,-6mastep) );
  label.llft( btex $\alpha$ etex scaled lscale, refpb );
  label.lrt( btex $\beta$ etex scaled lscale, refpc );  
endfig;

verbatimtex
}
\end{document}
etex

end.
