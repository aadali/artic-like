% directfonts.mp
% Check out /usr/share/texmf/fonts/tfm/public/
% or /usr/TeX/texmf/fonts/tfm/public

prologues := 3;

beginfig(0);
	numeric xCoord;
	xCoord = 3in;
	defaultfont := "bchr8r";
	label( "QWERTYUIOPasdfghjklZXCVBNM", (xCoord, 4cm) );
	label( "What gives?", (xCoord, 3cm) );
	defaultfont := "cmss17";
	label( "0.123,456,789", (xCoord, 2cm) );
%	defaultfont := "fplmr";
%	label( "áéíóúàèìòùãõçÁÉÍÓÚÀÈÌÒÙÃÕÇ°", (xCoord, 0) );
endfig;

end;

