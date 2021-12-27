input automata;

verbatimtex
% generic layout
\documentclass{article}
\usepackage[english]{babel}
\usepackage[latin1]{inputenc}
\selectlanguage{english}
% additional packages
\usepackage{amsmath}
\begin{document}
etex


% EXAMPLE FOR THE PACKAGE automata.mp

beginfig(0);

  % Setting of some global parameters.

  size := 15;
  incominglength := 25;
  loopsize := 22;

  % The following command defines a grid of nodes with 7 rows and 10 columns.
  % Instead of using "spacing := (36,29)" and thus change the global parameter for the ditance 
  % between nodes, one can use `with spacing((36,29))' to temporarily change the value of this 
  % parameter.

  with spacing((36,29)) matrix.a(7,10);
  
  % Next, one declares the shape of the nodes and possibly their labels. One can use the directives 
  % "btex" and "etex" to include LaTeX-compiled labels in nodes, or just "" to leave them unlabeled. 
  % Moreover, one can use several shapes for the nodes. There are some predefined shorhands for 
  % modified styles of nodes and arrows (e.g., "node_double" is equivalent to "with border(double) node")

  node.a[4][0](btex $s_0$ etex);
  with shape(boxed) node.a[4][2](btex $s_1$ etex);
  with shape(fixedbox) node.a[4][4](btex $s_2$ etex);
  with size(18) with shape(diamond) node.a[4][6](btex $s_3$ etex);
  with size(22) with shape(utriangle) node.a[4][8](btex $s_4$ etex);
  node.a[2][2](btex $s_5$ etex);
  node.a[2][5](btex $s_6$ etex);
  node.a[1][6](btex $s_7$ etex);
  node_double.a[1][8](btex $s_8$ etex);
 
  % One can define arrows of different forms. Each arrow follow a path from an initial node to a 
  % final node and this path must be explicitly given as a MetaPost path expression (e.g., point..point--point).
  % Note that labels of nodes and arrows can be rotated using the MetaPost directive "rotated".
  
  incoming(0, "") (a[4][0]) 120;
  loop.bot(.65, btex \small $a$ etex) (a[4][0]) -165;
  arrow_dash.bot(.5, btex \small $b$ etex) (a[4][0],a[4][2]) a[4][0].c..a[4][2].c;
  arrow.bot(.5, btex \small $a$ etex) (a[4][2],a[4][4]) a[4][2].c..a[4][4].c;
  arrow.bot(.5, btex \small $a$ etex) (a[4][4],a[4][6]) a[4][4].c..a[4][6].c;
  arrow.bot(.5, btex \small $a$ etex) (a[4][6],a[4][8]) a[4][6].c..a[4][8].c;
  arrow.bot(.5, btex \small $c$ etex) (a[4][8],a[4][0]) a[4][8].c..a[5][7].c---a[5][1].c..a[4][0].c;
  arrow_dash.rt(.5, btex \small $a$ etex rotated 90) (a[4][2],a[2][2]) a[4][2].c..a[2][2].c;
  arrow.bot(.5, btex \small $c$ etex rotated 39) (a[2][2],a[4][0]) a[2][2].c..a[4][0].c;
  arrow_dash.bot(.5, btex \small $a$ etex) (a[2][2],a[2][5]) a[2][2].c{dir-30}..a[2][5].c;
  arrow_dash.rt(.5, btex \small $a$ etex rotated 90) (a[4][6],a[1][6]) a[4][6].c..a[1][6].c;
  arrow.top(.5, btex \small $a$ etex) (a[1][6],a[1][8]) a[1][6].c..a[1][8].c;
  arrow.bot(.38, btex \small $c$ etex) (a[2][5],a[4][0]) a[2][5].c..a[1][4].c---a[1][2].c shifted(-10,0)..a[2][1].c shifted (-10,0)---a[4][0].c;
  arrow_bold.bot(.5, btex \small $c$ etex) (a[1][8],a[4][0]) a[1][8].c..a[2][9].c---a[5][9].c..a[6][8].c---a[6][1].c..a[5][0].c---a[4][0].c;
  arrow_bold.top(.38, btex \small $c$ etex) (a[1][8],a[4][0]) a[1][8].c..a[0][7].c---a[0][1].c..a[1][0].c---a[4][0].c;
  
  % Finally, one may want to highligh a given set of nodes. This can done using the macro
  % "enclose"/"enclosure", which draws/returns a path that goes around the given sequence 
  % of vertices (it is always a good idea if the sequence is given in clockwise order and
  % if it defines a convex region!).
  
  with enclosuremargin(20) draw enclosure(a[4][4].c--a[2][5].c--a[1][6].c--a[4][6].c--cycle) withcolor light.light.green withpen pencircle scaled 3;  

endfig;

end.