% defaultcmr.mp
% use dvips -Pcmz 

beginfig(0);
	picture sent[];
	numeric xCoord, i;
	color one, two, thr, fur;
	one = (0.13333,0.54510,0.13334);
	two = (0.74118,0.71764,0.41961);
	thr = (0.80000,0.20000,0.20000);
	fur = (0.00000,0.85000,1.00000);
	xCoord = 3in;

	defaultfont := "cmbx12";
	defaultscale := 6.5;

	sent1 = thelabel( " ", (xCoord, 6in) );
	sent2 = thelabel( "Be", (xCoord, 5in) );
	sent3 = thelabel( "Aware", (xCoord, 4in) );
	sent4 = thelabel( "of", (xCoord, 3in) );
	sent5 = thelabel( "Technical", (xCoord, 2in) );
	sent6 = thelabel( "Arrogance.", (xCoord, 1in) );
	for i = 1 upto 4:
		draw sent[i] withcolor two;
	endfor;
	for i = 5 upto 6:
		draw sent[i] withcolor thr;
	endfor;
endfig;

end;

