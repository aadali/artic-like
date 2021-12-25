
beginfig(1);
	picture topl;
	color one, two, thr, fur;
	string drawcommand;
	drawcommand = "draw topl withcolor one;";
%	defaultfont := "putr8r";
%	defaultfont := "bchr8r";
	one = (0.13333,0.54510,0.13334);
	two = (0.74118,0.71764,0.41961);
	thr = (0.80000,0.20000,0.20000);
	fur = (0.00000,0.85000,1.00000);
	topl = thelabel( "40", origin) scaled 7;
	scantokens( drawcommand );
endfig;

end;
