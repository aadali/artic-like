input boxes;

beginfig(0);
	picture topl, midd, botl;
	numeric tl, ml, bl;
	defaultfont := "bchbo8r";
%	defaultfont := "putro8r";
	topl = thelabel( "Eu Vivo", origin) yscaled 1;
	midd = thelabel( "na", origin) scaled 1;
	botl = thelabel( "Cauda da Europa", origin) yscaled 1;
	tl = abs( urcorner topl - ulcorner topl );
	ml = abs( urcorner midd - ulcorner midd );
	bl = abs( urcorner botl - ulcorner botl );
	boxjoin( a.se=b.ne; a.sw=b.nw  );
	boxit.topli(topl);
	boxit.middl(midd);
	boxit.botli(botl);
	drawunboxed( topli, middl, botli );
endfig;

end;
