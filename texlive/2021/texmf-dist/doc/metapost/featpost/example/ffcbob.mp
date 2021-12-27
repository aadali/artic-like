% ffcbob.mp
% L. Nobre G.
% 2003

input featpost3Dplus2D;

f := 2*(15,0,2);
Spread := 20;

beginfig(1);
	numeric phi, theta, back, numbacks, minnum, maxnum;
	phi = 30;
	theta = 80;
	back = 7;
	numbacks = 30;
	minnum = 3;
	maxnum = 15;
	numeric upl[], dol[], pos[], len[], rai[] rao[];
	upl1=1;  dol1=17;pos1=-5.0; len1=2.0; rai1=1.5; rao1=6.5;
	upl2=2;  dol2=16;pos2=-3.0; len2=0.8; rai2=1.5; rao2=6.3;
	upl3=3;  dol3=15;pos3=-2.2; len3=0.4; rai3=1.5; rao3=4.8;
	upl4=4;  dol4=14;pos4=-1.8; len4=0.4; rai4=1.5; rao4=4.4;
	upl5=5;  dol5=13;pos5=-1.4; len5=0.4; rai5=1.5; rao5=4.2;
	upl6=6;  dol6=11;pos6=-1.0; len6=0.2; rai6=1.5; rao6=3.2;
	upl7=7;  dol7=12;pos7=-1.0; len7=0.2; rai7=4.0; rao7=4.6;
	upl8=8;  dol8=10;pos8=-0.8; len8=0.6; rai8=1.5; rao8=3.2;
	upl9=9;  dol9=7; pos9=-0.2; len9=0.4; rai9=1.5; rao9=3.2;
	upl10=10;dol10=9;pos10=-0.2;len10=0.2;rai10=4.0;rao10=4.4;
	upl11=11;dol11=8;pos11=0.0; len11=0.2;rai11=4.0;rao11=4.6;
	upl12=12;dol12=6;pos12=0.2; len12=0.8;rai12=1.5;rao12=4.8;
	upl13=13;dol13=5;pos13=1.0; len13=0.4;rai13=1.5;rao13=4.4;
	upl14=14;dol14=3;pos14=1.4; len14=0.2;rai14=1.5;rao14=4.4;
	upl15=15;dol15=4;pos15=1.4; len15=0.2;rai15=4.8;rao15=6.3;
	upl16=16;dol16=2;pos16=1.6; len16=0.6;rai16=1.5;rao16=6.3;
	upl17=17;dol17=1;pos17=2.2; len17=2.0;rai17=1.5;rao17=6.5;
	color axe;
	axe = (cosd(phi)*cosd(theta),sind(phi)*cosd(theta),sind(theta));
	numeric i;
	for i=0 upto numbacks:
	    draw rp(back*(1,2*i/numbacks-1,-1))--
		 rp(back*(-1,2*i/numbacks-1,-1))--
		 rp(back*(-1,2*i/numbacks-1,1)) withcolor 0.7background;
	endfor;
	numeric crux;
	crux = cdotprod( f, axe );
	numeric j, k;
	if crux <= pos[minnum+1]:
	    for i=minnum upto maxnum:
		j := maxnum;
		forever:
		    exitunless dol[j] <> i;
		    j := decr( j );
		endfor;
%		rigorousdisc(rai[j],true,axe*pos[j],rao[j],axe*len[j]);
		rigorousdisc(0,true,axe*pos[j],rao[j],axe*len[j]);
	    endfor;
	elseif crux > pos[maxnum]:
	    for i=minnum upto maxnum:
%		rigorousdisc(rai[i],true,axe*pos[i],rao[i],axe*len[i]);
		rigorousdisc(0,true,axe*pos[i],rao[i],axe*len[i]);
	    endfor;
	else:
	    k := minnum+2;
	    forever:
		exitif crux <= pos[k];
		k := incr( k );
	    endfor;
	    for i=minnum upto k-1:
%		rigorousdisc(rai[i],true,axe*pos[i],rao[i],axe*len[i]);
		rigorousdisc(0,true,axe*pos[i],rao[i],axe*len[i]);
	    endfor;
	    for i=k upto maxnum:
		j := maxnum;
		forever:
		    exitunless dol[j] <> i;
		    j := decr( j );
		endfor;
%		rigorousdisc(rai[j],true,axe*pos[j],rao[j],axe*len[j]);
		rigorousdisc(0,true,axe*pos[j],rao[j],axe*len[j]);
	    endfor;
	fi;
endfig;

end.
