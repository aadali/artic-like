% statement.mp
% L. Nobre G.
% 2007

% prologues := 3;

input boxes;

verbatimtex \documentclass{article}\usepackage{palatino}\begin{document} etex

beginfig(1);
  picture topl, midd, botl;
  numeric tl, ml, bl;
%	topl = thelabel(btex \mbox{G E T T ' N ' R E A D Y} etex, origin);
%	midd = thelabel(btex \mbox{4D} etex, origin);
%	botl = thelabel(btex \mbox{c l i m b} etex, origin);
%	defaultfont := "putr8r";
%	defaultfont := "bchr8r";
  topl = thelabel( "G E T T ' N ' R E A D Y", origin) yscaled 0.3;
  midd = thelabel( "4D", origin) scaled 1;
  botl = thelabel( "C   L   I   M   B", origin) yscaled 0.6;
  tl = abs( urcorner topl - ulcorner topl );
  ml = abs( urcorner midd - ulcorner midd );
  bl = abs( urcorner botl - ulcorner botl );
  boxjoin( a.se=b.ne; a.sw=b.nw  );
  boxit.topli(topl);
  boxit.middl(midd scaled (tl/ml) );
  boxit.botli(botl xscaled (tl/bl) );
  drawunboxed( topli, middl, botli );
endfig;

beginfig(2);
  picture topl, midd, botl;
  numeric tl, ml, bl, sc;
  sc = 10;
  topl = thelabel( btex o\ n\ e\ \ l\ i\ f\ e etex scaled sc, origin) yscaled 0.3;
  midd = thelabel( btex {\huge 40} etex scaled sc, origin) scaled 1;
  botl = thelabel( btex D\ R\ I\ L\ L\ S etex scaled sc, origin) yscaled 0.6;
  tl = abs( urcorner topl - ulcorner topl );
  ml = abs( urcorner midd - ulcorner midd );
  bl = abs( urcorner botl - ulcorner botl );
  boxjoin( a.se=b.ne; a.sw=b.nw  );
  boxit.topli(topl);
  boxit.middl(midd scaled (tl/ml) );
  boxit.botli(botl xscaled (tl/bl) );
  drawunboxed( topli, middl, botli );
endfig;

verbatimtex \end{document} etex

end;
