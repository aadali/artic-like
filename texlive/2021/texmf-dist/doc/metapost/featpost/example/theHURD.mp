% Converted from PostScript(TM) to MetaPost by pstoedit
% MetaPost backend contributed by Scott Pakin <pakin_AT_uiuc.edu>
% pstoedit is Copyright (C) 1993 - 2006 Wolfgang Glunz <wglunz34_AT_pstoedit.net>

% Generate structured PostScript
prologues := 1;

% Display a given string with its *baseline* at a given location
% and with a given rotation angle
vardef showtext(expr origin)(expr angle)(expr string) =
  draw string infont defaultfont scaled defaultscale
    rotated angle shifted origin;
enddef;

% Modified by L. Nobre G's Lex scanner %%%%%%%%%%% 
randomseed := 1010; 
vardef renewcolore = 
  color tmpcolo; 
  tmpcolo = (uniformdeviate(0.8),uniformdeviate(0.8),uniformdeviate(0.8)); 
  drawoptions( withcolor tmpcolo ); 
enddef; 
defaultscale := 0.16; 
beginfig(1);
 path PreeditedPath[];
linecap := butt;
linejoin := mitered;
 PreeditedPath1 = hide(z1=(9.82422,68.9492))z1--hide(z2=(9.82422,12.8594))z2--hide(z3=(17.6836,11.7891))z3
	..controls hide(z4=(21.9688,11.2031))z4 and hide(z5=(26.6367,10.7148))z5..hide(z6=(27.9531,10.7148))z6
	..controls hide(z7=(29.8477,10.7148))z7 and hide(z8=(30.3633,15.4492))z8..hide(z9=(30.3633,32.7773))z9
	--hide(z10=(30.3633,54.8359))z10--hide(z11=(55.8203,54.6562))z11..controls hide(z12=(69.75,54.5586))z12 and hide(z13=(83.5781,55.125))z13..hide(z14=(86.3633,55.9062))z14
	--hide(z15=(91.4531,57.3359))z15--hide(z16=(90.2031,38.0469))z16..controls hide(z17=(89.5195,27.5078))z17 and hide(z18=(88.6602,18.2656))z18..hide(z19=(88.3125,17.6836))z19
	..controls hide(z20=(86.7109,15.0156))z20 and hide(z21=(91.2773,12.7227))z21..hide(z22=(100.387,11.6172))z22
	--hide(z23=(110.746,10.3594))z23--hide(z24=(110.477,65.0195))z24--hide(z25=(110.211,119.676))z25
	--hide(z26=(101.191,122.355))z26..controls hide(z27=(96.1875,123.844))z27 and hide(z28=(91.4531,125.035))z28..hide(z29=(90.5625,125.035))z29
	..controls hide(z30=(89.625,125.035))z30 and hide(z31=(89.4766,115.391))z31..hide(z32=(90.2031,101.992))z32
	--hide(z33=(91.4531,78.9492))z33--hide(z34=(82.7031,77.5547))z34..controls hide(z35=(71.9844,75.8438))z35 and hide(z36=(45.9961,75.8281))z36..hide(z37=(36.9727,77.5273))z37
	--hide(z38=(30.3633,78.7734))z38--hide(z39=(30.3633,99.1367))z39--hide(z40=(30.3633,119.5))z40
	--hide(z41=(21.0742,122.266))z41..controls hide(z42=(15.8945,123.812))z42 and hide(z43=(11.3516,125.035))z43..hide(z44=(10.8047,125.035))z44
	..controls hide(z45=(10.2617,125.035))z45 and hide(z46=(9.82422,99.8516))z46..cycle;
% renewcolore;
 fill PreeditedPath1;
% dotlabels.rt(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46);
 PreeditedPath2 = hide(z47=(126.375,82.7031))z47..controls hide(z48=(126.719,59.3008))z48 and hide(z49=(127.871,37.457))z49..hide(z50=(128.941,34.1172))z50
	..controls hide(z51=(131.738,25.3555))z51 and hide(z52=(141.648,15.3555))z52..hide(z53=(152.367,10.4766))z53
	..controls hide(z54=(159.797,7.09375))z54 and hide(z55=(164.691,6.26562))z55..hide(z56=(176.836,6.32031))z56
	..controls hide(z57=(194.344,6.39844))z57 and hide(z58=(204.965,10.4414))z58..hide(z59=(215.602,21.0742))z59
	..controls hide(z60=(226.223,31.6953))z60 and hide(z61=(226.871,35.0078))z61..hide(z62=(227.145,80.0234))z62
	--hide(z63=(227.391,120.57))z63--hide(z64=(220.422,121.91))z64..controls hide(z65=(216.672,122.633))z65 and hide(z66=(211.602,123.664))z66..hide(z67=(209.438,124.145))z67
	--hide(z68=(205.418,125.035))z68--hide(z69=(205.418,84.668))z69..controls hide(z70=(205.418,37.1523))z70 and hide(z71=(204.617,33.9414))z71..hide(z72=(191.305,28.043))z72
	..controls hide(z73=(175.586,21.0781))z73 and hide(z74=(156.254,24.7109))z74..hide(z75=(150.438,35.7227))z75
	..controls hide(z76=(148.648,39.1133))z76 and hide(z77=(148.25,49.9258))z77..hide(z78=(148.793,80.6484))z78
	--hide(z79=(149.508,121.105))z79--hide(z80=(139.684,122.355))z80..controls hide(z81=(134.324,123.039))z81 and hide(z82=(128.926,123.973))z82..hide(z83=(127.805,124.41))z83
	..controls hide(z84=(126.285,125.004))z84 and hide(z85=(125.914,114.141))z85..cycle;
% renewcolore;
 fill PreeditedPath2;
% dotlabels.rt(47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85);
 PreeditedPath3 = hide(z86=(241.305,122.18))z86..controls hide(z87=(243.297,119.777))z87 and hide(z88=(243.898,107.887))z88..hide(z89=(244.18,65.375))z89
	..controls hide(z90=(244.375,35.7227))z90 and hide(z91=(244.578,11.0898))z91..hide(z92=(244.625,10.4492))z92
	..controls hide(z93=(244.676,9.78906))z93 and hide(z94=(247.039,8.82031))z94..hide(z95=(250.074,8.21484))z95
	..controls hide(z96=(253.109,7.60938))z96 and hide(z97=(257.754,6.67969))z97..hide(z98=(260.793,6.07031))z98
	--hide(z99=(266.152,5))z99--hide(z100=(266.152,26.6133))z100..controls hide(z101=(266.152,38.582))z101 and hide(z102=(266.707,48.2188))z102..hide(z103=(267.402,48.207))z103
	..controls hide(z104=(275.438,48.0781))z104 and hide(z105=(299.418,38.2695))z105..hide(z106=(306.164,32.3555))z106
	..controls hide(z107=(310.332,28.6953))z107 and hide(z108=(316.168,17.3242))z108..hide(z109=(316.168,12.8594))z109
	..controls hide(z110=(316.168,10.3594))z110 and hide(z111=(317.773,8.96484))z111..hide(z112=(321.523,8.21484))z112
	..controls hide(z113=(324.418,7.63672))z113 and hide(z114=(329.473,6.64844))z114..hide(z115=(332.512,6.07031))z115
	--hide(z116=(338.137,5))z116--hide(z117=(337.016,13.2148))z117..controls hide(z118=(335.402,25.0039))z118 and hide(z119=(330.18,34.0078))z119..hide(z120=(319.023,44.2305))z120
	--hide(z121=(309.199,53.2305))z121--hide(z122=(314.379,55.3711))z122..controls hide(z123=(321.926,58.4883))z123 and hide(z124=(333.449,71.6289))z124..hide(z125=(335.691,79.6641))z125
	..controls hide(z126=(339.273,92.5273))z126 and hide(z127=(338.176,101.074))z127..hide(z128=(332.043,108.066))z128
	..controls hide(z129=(324.508,116.656))z129 and hide(z130=(312.594,119.887))z130..hide(z131=(288.121,119.969))z131
	..controls hide(z132=(277.582,120.004))z132 and hide(z133=(263.32,121.176))z133..hide(z134=(256.684,122.547))z134
	..controls hide(z135=(240.07,125.984))z135 and hide(z136=(238.184,125.941))z136..cycle;
% renewcolore;
 fill PreeditedPath3;
% dotlabels.rt(86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136);
 PreeditedPath4 = hide(z137=(308.305,101.707))z137..controls hide(z138=(316.125,99.1367))z138 and hide(z139=(317.668,93.5977))z139..hide(z140=(313.555,82.8789))z140
	..controls hide(z141=(307.645,67.4766))z141 and hide(z142=(288.836,60.0312))z142..hide(z143=(272.762,66.7305))z143
	--hide(z144=(266.152,69.4844))z144--hide(z145=(266.152,85.3828))z145..controls hide(z146=(266.152,94.1328))z146 and hide(z147=(266.551,101.672))z147..hide(z148=(267.043,102.156))z148
	..controls hide(z149=(268.93,104.012))z149 and hide(z150=(302.41,103.648))z150..cycle;
% renewcolore;
unfill PreeditedPath4;
% dotlabels.rt(137,138,139,140,141,142,143,144,145,146,147,148,149,150);
 PreeditedPath5 = hide(z151=(354.57,65.6445))z151--hide(z152=(354.57,7.14453))z152--hide(z153=(376.898,7.14453))z153
	..controls hide(z154=(401.727,7.14453))z154 and hide(z155=(412.688,9.31641))z155..hide(z156=(426.914,17.0547))z156
	..controls hide(z157=(438.168,23.1797))z157 and hide(z158=(445.777,31.6133))z158..hide(z159=(451.641,44.4766))z159
	..controls hide(z160=(466.125,76.2461))z160 and hide(z161=(444.777,112.316))z161..hide(z162=(406.73,120.355))z162
	..controls hide(z163=(399.598,121.859))z163 and hide(z164=(385.117,123.309))z164..hide(z165=(374.219,123.605))z165
	--hide(z166=(354.57,124.145))z166--cycle;
% renewcolore;
fill PreeditedPath5;
% dotlabels.rt(151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166);
 PreeditedPath6 = hide(z167=(413.695,95.0625))z167..controls hide(z168=(428.523,87.418))z168 and hide(z169=(436.469,71.9844))z169..hide(z170=(433.215,57.1602))z170
	..controls hide(z171=(431.301,48.418))z171 and hide(z172=(419.41,36.4961))z172..hide(z173=(409.23,33.1016))z173
	..controls hide(z174=(399.941,30.0078))z174 and hide(z175=(383.602,28.1562))z175..hide(z176=(379.398,29.7227))z176
	..controls hide(z177=(377.191,30.5469))z177 and hide(z178=(376.766,36.6172))z178..hide(z179=(376.988,64.3047))z179
	..controls hide(z180=(377.137,82.7031))z180 and hide(z181=(377.793,98.4492))z181..hide(z182=(378.445,99.1367))z182
	..controls hide(z183=(381.016,101.84))z183 and hide(z184=(406.195,98.9297))z184..cycle;
% renewcolore;
unfill PreeditedPath6;
% dotlabels.rt(167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184);

      endfig;
end
