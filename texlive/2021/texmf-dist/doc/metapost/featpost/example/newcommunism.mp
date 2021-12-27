% newcommunism.mp
% L. Nobre G.
% 2007

prologues := 2;

beginfig(1);
  numeric u, elhei, elwid, elpos, bawid, bahei, thstp, thdel, i;
  numeric refang, reft, tfin, tini, spikang;
  path el, ba, aux, cutel, horns;
  pair reftpair;
  color yellow, thered;
  yellow = 0.87*(0.95*red+0.7*green);
  thered = 0.66red;
  u = 2mm;
  thdel = 25;
  thstp = 5;
  elhei = 30u;
  elwid = 44u;
  elpos = 26u;
  bawid = 27u;
  bahei = 10u;
  el = (elwid,elpos)..for i=thstp step thstp until 360-thstp:
    (elwid*cosd(i),elpos+elhei*sind(i)).. endfor cycle;
  ba = (bawid,-bahei){right}..{left}(bawid,0)--
    (elwid*cosd(90-thdel),elpos+elhei*sind(90-thdel))--
    (bawid/3.0,0)--(0,elpos+elhei)--(-bawid/3.0,0)--
    (elwid*cosd(90+thdel),elpos+elhei*sind(90+thdel))--
    (-bawid,0){left}..{right}(-bawid,-bahei)--cycle;
  aux = (0,elpos)--
    1.1[(0,elpos),(elwid*cosd(270+2thdel),elpos+elhei*sind(270+2thdel))];
  reftpair = el intersectiontimes aux;
  reft = xpart reftpair;
  refang = angle( direction reft of el );
  tfin = directiontime dir(180-refang) of el;
  tini = directiontime dir(180+refang) of el - length el;
  cutel = subpath (tini,tfin) of el;
  spikang = (90-angle((0,elpos+elhei)-(-bawid,0)));
  horns = cutel{dir(-refang-spikang)}..
    {left}(0,elpos-elhei+bawid/3.0)..
    {dir(refang+spikang)}cycle;
  fill ba withcolor yellow;
  fill horns withcolor yellow;
  picture storepic;
  storepic = currentpicture;
  currentpicture := nullpicture;
  bboxmargin := bahei;
  fill bbox storepic withcolor thered;
  draw storepic;  
endfig;

beginfig(2);
	numeric u;
	color yellow, thered;
	path martelo, foice, monitor, ecran, mar, foi;
	u = 1pt;
%	yellow = 0.87*(0.95*red+0.7*green);
	thered = 0.66red;
	yellow = 0.7blue;
	picture storemf;

	numeric cabosal, cabesal, cabolar, cabelar;
	numeric cabolen, cabeclen, cabecontrol, martelrot;
	pair martelpos;
	cabosal = 2u;
	cabesal = 6u;
	cabolar = 8u;
	cabelar = 10u;
	cabolen = 70u;
	cabeclen = 7u;
	cabecontrol = 4u;
	martelrot = 50;
	martelpos = dir(180)*0.03*cabolen;
	z100 = origin;
	z101 = ( cabesal, 0 );
	z102 = ( cabesal, cabosal+cabelar-cabolen );
	z103 = ( cabesal+cabolar, cabosal+cabelar-cabolen );
	z104 = ( cabesal+cabolar, 0 );
	z105 = ( cabesal+cabolar+cabeclen, 0 );
	z106 = ( cabesal+cabolar+cabecontrol, cabelar );
	z107 = ( cabesal+cabolar, cabelar );
	z108 = ( cabesal+cabolar, cabelar+cabosal );
	z109 = ( cabesal, cabelar+cabosal );
	z110 = ( cabesal, cabelar );
	z111 = ( 0, cabelar );
	martelo = z100--
		  z101--
		  z102--
		  z103--
		  z104--
		  z105..controls z106..
		  z107--
		  z108--
		  z109--
		  z110--
		  z111--cycle;
	mar = martelo 
		 rotated martelrot 
		 shifted martelpos;

	filldraw mar withcolor yellow;
	show mar;	 
	numeric fcabolen, fcabolar, foicecontrol, laminang;
	numeric foicerot, foicaber, flamilar, fcabomar, refang;
	numeric fmaomarg, fmaorais, fmaospac, fprotrud, basangul;
	pair foicepos;
	laminang = 3;
	refang = 10
	basangul = 23;
	fcabomar = 3u;
	fcabolar = 8u - fcabomar;
	flamilar = 4u;
	fcabolen = 22u;
	foicecontrol = 2u;
	fmaomarg = 6u;
	fmaorais = 2u;
	fmaospac = 14u;
	fprotrud = 30u;
	foicaber = 42u;
	foicerot = -50;
	foicepos = dir(270)*0.3*foicaber;
	z200 = origin;
	z201 = ( -fcabomar, 0 );
	z202 = ( -fcabomar, -fcabolen );
	z203 = ( fcabolar, -fcabolen );
	z204 = ( fcabolar+foicecontrol, -fcabolen );
	z205 = ( fcabolar+fmaomarg, fmaorais-fmaospac );
	z206 = ( fcabolar, -fmaospac );
	z207 = ( fcabolar, 0 );
	z208 = ( -fprotrud, foicaber );
	z209 = ( 0, flamilar );
	foice = z200--
		z201--
		z202--
		z203..controls z204..
		z205{dir(-160)}..{left}
		z206--
		z207{dir(basangul)}..tension 0.8..{dir(180+refang+laminang)}
		z208{dir(refang-laminang)}..tension 0.8..{dir(180+basangul)}
		z209--cycle;
	foi = foice
		 rotated foicerot 
		 shifted foicepos;
		 
	filldraw foi  withcolor yellow;
	show foi;
	storemf = currentpicture shifted -(center currentpicture);	
	currentpicture := nullpicture;





	
	numeric monith, monitw;
	pair mll, mlr, mur, mul;
	monith = 18u + ypart ((urcorner storemf)-(llcorner storemf));
	monitw = 1.34monith;
	mll = (-0.5monitw, -0.5monith);
	mlr = (0.5monitw, -0.5monith);
	mur = (0.5monitw, 0.5monith);
	mul = (-0.5monitw, 0.5monith);
	ecran = mll--
		mlr--
		mur--
		mul--cycle;
	
	numeric bordthick, roundfrac, alturpe, borraspac, largupe;
	numeric borralarg, apoiolarg, cornermarg, apoiothick;
	bordthick = 13u;
	roundfrac = 0.3;
	alturpe = 15u;
	borraspac = 0.95u;
	largupe = 50u;
	borralarg = 60u;
	apoiolarg = 156u;
	cornermarg = 2u;
	apoiothick = 18u;
	z301 = (0.5largupe,-0.5monith-bordthick);
	z302 = (0.5monitw+roundfrac*bordthick,-0.5monith-bordthick);
	z303 = (0.5monitw+bordthick,-0.5monith-roundfrac*bordthick);
	z304 = (0.5monitw+bordthick,0.5monith+roundfrac*bordthick);
	z305 = (0.5monitw+roundfrac*bordthick,0.5monith+bordthick);
	z306 = (-0.5monitw-roundfrac*bordthick,0.5monith+bordthick);
	z307 = (-0.5monitw-bordthick,0.5monith+roundfrac*bordthick);
	z308 = (-0.5monitw-bordthick,-0.5monith-roundfrac*bordthick);
	z309 = (-0.5monitw-roundfrac*bordthick,-0.5monith-bordthick);
	z310 = (-0.5largupe,-0.5monith-bordthick);
	z311 = (-0.5largupe,-0.5monith-bordthick-alturpe);
	z312 = (-0.5apoiolarg+cornermarg,
		-0.5monith-bordthick-alturpe-cornermarg);
	z313 = (-0.5apoiolarg,
		-0.5monith-bordthick-alturpe-apoiothick-borraspac);
	z314 = (-0.5borralarg,
		-0.5monith-bordthick-alturpe-apoiothick-borraspac);
	z315 = (0,-0.5monith-bordthick-alturpe-apoiothick);
	z316 = (0.5borralarg,
		-0.5monith-bordthick-alturpe-apoiothick-borraspac);
	z317 = (0.5apoiolarg,
		-0.5monith-bordthick-alturpe-apoiothick-borraspac);
	z318 = (0.5apoiolarg-cornermarg,
		-0.5monith-bordthick-alturpe-cornermarg);
	z319 = (0.5largupe,-0.5monith-bordthick-alturpe);
	monitor = z301--
		  z302{right}..{up}
		  z303--
		  z304{up}..{left}
		  z305--
		  z306{left}..{down}
		  z307--
		  z308{down}..{right}
		  z309--
		  z310--
		  z311..controls z312..
		  z313--
		  z314{dir(73)}...{right}
		  z315{right}...{dir(-73)}
		  z316--
		  z317..controls z318..
		  z319--cycle;
	draw monitor withcolor thered;
	bboxmargin := monitw;
	fill bbox monitor withcolor thered;
	filldraw monitor withcolor yellow;
	filldraw ecran withcolor thered;
	show monitor;
	show ecran;
%	label.urt("New Communism", llcorner monitor);
	draw storemf;


      endfig;


end;


