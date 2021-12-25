% Converted from PostScript(TM) to MetaPost by pstoedit
% MetaPost backend contributed by Scott Pakin <pakin@uiuc.edu>
% pstoedit is Copyright (C) 1993 - 2003 Wolfgang Glunz <wglunz@pstoedit.net>

% Generate structured PostScript
prologues := 1;

% Display a given string with its *baseline* at a given location
% and with a given rotation angle
vardef showtext(expr origin)(expr angle)(expr string) =
  draw string infont defaultfont scaled defaultscale
    rotated angle shifted origin;
enddef;

% Modified by L. Nobre G's Lex scanner %%%%%%%%%%% 

beginfig(1);
 path PreeditedPath[];
 color bluecolor;
 bluecolor = (0.607845,0.61569,0.823535);

 linecap := butt;
 linejoin := mitered;

 PreeditedPath1 = hide(z1=(26.4846,113.11))z1..controls hide(z2=(47.8035,130.71))z2 and hide(z3=(67.1843,138.71))z3..hide(z4=(82.689,138.71))z4
	..controls hide(z5=(98.1936,138.71))z5 and hide(z6=(123.389,117.91))z6..hide(z7=(144.708,117.91))z7
	..controls hide(z8=(164.088,117.91))z8 and hide(z9=(173.719,138.71))z9..hide(z10=(204.788,138.71))z10
	..controls hide(z11=(233.859,138.71))z11 and hide(z12=(241.615,116.244))z12..hide(z13=(290.064,117.91))z13
	..controls hide(z14=(336.578,119.51))z14 and hide(z15=(342.304,140.311))z15..hide(z16=(363.711,140.311))z16
	..controls hide(z17=(385.03,140.311))z17 and hide(z18=(408.181,129.21))z18..hide(z19=(419.916,111.51))z19
	..controls hide(z20=(400.535,125.91))z20 and hide(z21=(376.046,119.227))z21..hide(z22=(354.021,109.91))z22
	..controls hide(z23=(328.981,99.3196))z23 and hide(z24=(295.878,85.9104))z24..hide(z25=(268.745,95.5105))z25
	..controls hide(z26=(241.612,105.11))z26 and hide(z27=(214.479,109.91))z27..hide(z28=(191.222,101.91))z28
	..controls hide(z29=(167.965,93.9104))z29 and hide(z30=(127.265,85.9104))z30..hide(z31=(100.132,101.91))z31
	..controls hide(z32=(72.9985,117.91))z32 and hide(z33=(55.5559,121.11))z33..cycle;

 filldraw PreeditedPath1 withcolor bluecolor;

 PreeditedPath3 = hide(z67=(101.492,205.833))z67..controls hide(z68=(114.075,219.402))z68 and hide(z69=(125.514,225.569))z69..hide(z70=(134.666,225.569))z70
	..controls hide(z71=(143.817,225.569))z71 and hide(z72=(158.688,209.533))z72..hide(z73=(171.271,209.533))z73
	..controls hide(z74=(182.71,209.533))z74 and hide(z75=(188.394,225.569))z75..hide(z76=(206.733,225.569))z76
	..controls hide(z77=(223.892,225.569))z77 and hide(z78=(228.469,208.249))z78..hide(z79=(257.065,209.533))z79
	..controls hide(z80=(284.52,210.767))z80 and hide(z81=(287.899,226.803))z81..hide(z82=(300.535,226.803))z82
	..controls hide(z83=(313.118,226.803))z83 and hide(z84=(326.782,218.245))z84..hide(z85=(333.708,204.599))z85
	..controls hide(z86=(322.269,215.701))z86 and hide(z87=(307.815,210.548))z87..hide(z88=(294.815,203.366))z88
	..controls hide(z89=(280.036,195.201))z89 and hide(z90=(260.497,184.863))z90..hide(z91=(244.482,192.264))z91
	..controls hide(z92=(228.467,199.665))z92 and hide(z93=(212.452,203.366))z93..hide(z94=(198.725,197.198))z94
	..controls hide(z95=(184.998,191.031))z95 and hide(z96=(160.976,184.863))z96..hide(z97=(144.961,197.198))z97
	..controls hide(z98=(128.946,209.533))z98 and hide(z99=(118.651,212))z99..cycle;

 filldraw PreeditedPath3 withcolor bluecolor;

 PreeditedPath4 = hide(z166=(217.85,207.797))z166
	..controls hide(z167=(236.119,207.797))z167 and hide(z168=(249.691,188.745))z168..hide(z169=(249.691,162.646))z169
	..controls hide(z170=(249.691,137.069))z170 and hide(z171=(236.38,119.583))z171..hide(z172=(217.328,119.583))z172
	..controls hide(z173=(200.103,119.583))z173 and hide(z174=(188.097,137.069))z174..hide(z175=(188.097,162.646))z175
	..controls hide(z176=(188.097,189.267))z176 and hide(z177=(200.103,207.797))z177..cycle;

 clip currentpicture to PreeditedPath4;
 picture gthrough;
 gthrough = currentpicture;
 filldraw PreeditedPath1 withcolor bluecolor;
 filldraw PreeditedPath3 withcolor bluecolor;

 PreeditedPath5 = hide(z133=(249.169,235.984))z133--hide(z134=(249.169,214.322))z134..controls hide(z135=(238.207,231.025))z135 and hide(z136=(226.201,238.333))z136..hide(z137=(210.02,238.333))z137
	..controls hide(z138=(176.614,238.333))z138 and hide(z139=(151.558,205.709))z139..hide(z140=(151.558,162.124))z140
	..controls hide(z141=(151.558,118.278))z141 and hide(z142=(174.526,89.0466))z142..hide(z143=(209.237,89.0466))z143
	..controls hide(z144=(225.679,89.0466))z144 and hide(z145=(235.597,94.0056))z145..hide(z146=(249.169,109.404))z146
	--hide(z147=(249.169,90.3516))z147..controls hide(z148=(249.169,73.3872))z148 and hide(z149=(236.38,61.6428))z149..hide(z150=(218.111,61.6428))z150
	..controls hide(z151=(204.278,61.6428))z151 and hide(z152=(195.144,67.3845))z152..hide(z153=(192.273,78.0852))z153
	--hide(z154=(154.429,78.0852))z154..controls hide(z155=(154.69,67.1235))z155 and hide(z156=(158.866,59.5549))z156..hide(z157=(168.262,51.7251))z157
	..controls hide(z158=(179.484,42.5906))z158 and hide(z159=(195.144,38.1536))z159..hide(z160=(216.806,38.1536))z160
	..controls hide(z161=(259.347,38.1536))z161 and hide(z162=(283.88,57.2061))z162..hide(z163=(283.88,90.3516))z163
	--hide(z164=(283.88,235.984))z164--hide(z165=(249.169,235.984))z165--cycle;

 filldraw PreeditedPath5;

 unfill PreeditedPath4;	
 draw gthrough;	

 PreeditedPath6 = hide(z178=(58.4002,162.71))z178..controls hide(z179=(76,180.311))z179 and hide(z180=(92,188.311))z180..hide(z181=(104.8,188.311))z181
	..controls hide(z182=(117.6,188.311))z182 and hide(z183=(138.4,167.51))z183..hide(z184=(156,167.51))z184
	..controls hide(z185=(172,167.51))z185 and hide(z186=(179.95,188.311))z186..hide(z187=(205.6,188.311))z187
	..controls hide(z188=(229.6,188.311))z188 and hide(z189=(236.003,165.844))z189..hide(z190=(276,167.51))z190
	..controls hide(z191=(314.4,169.11))z191 and hide(z192=(319.127,189.91))z192..hide(z193=(336.8,189.91))z193
	..controls hide(z194=(354.4,189.91))z194 and hide(z195=(373.512,178.81))z195..hide(z196=(383.2,161.11))z196
	..controls hide(z197=(367.2,175.51))z197 and hide(z198=(346.983,168.826))z198..hide(z199=(328.8,159.51))z199
	..controls hide(z200=(308.128,148.92))z200 and hide(z201=(280.8,135.51))z201..hide(z202=(258.4,145.11))z202
	..controls hide(z203=(236,154.71))z203 and hide(z204=(213.6,159.51))z204..hide(z205=(194.4,151.51))z205
	..controls hide(z206=(175.2,143.51))z206 and hide(z207=(141.6,135.51))z207..hide(z208=(119.2,151.51))z208
	..controls hide(z209=(96.8001,167.51))z209 and hide(z210=(82.4002,170.71))z210..cycle;

 filldraw PreeditedPath6 withcolor bluecolor;

 PreeditedPath8 = hide(z244=(76,61.8328))z244..controls hide(z245=(93.6001,75.4011))z245 and hide(z246=(109.6,81.5686))z246..hide(z247=(122.4,81.5686))z247
	..controls hide(z248=(135.2,81.5686))z248 and hide(z249=(156,65.5332))z249..hide(z250=(173.6,65.5332))z250
	..controls hide(z251=(189.6,65.5332))z251 and hide(z252=(197.55,81.5686))z252..hide(z253=(223.2,81.5686))z253
	..controls hide(z254=(247.2,81.5686))z254 and hide(z255=(253.603,64.2485))z255..hide(z256=(293.6,65.5332))z256
	..controls hide(z257=(332,66.7668))z257 and hide(z258=(336.728,82.802))z258..hide(z259=(354.4,82.802))z259
	..controls hide(z260=(372,82.802))z260 and hide(z261=(391.112,74.2446))z261..hide(z262=(400.8,60.5994))z262
	..controls hide(z263=(384.8,71.7007))z263 and hide(z264=(364.583,66.5479))z264..hide(z265=(346.4,59.366))z265
	..controls hide(z266=(325.728,51.2012))z266 and hide(z267=(298.4,40.8635))z267..hide(z268=(276,48.2644))z268
	..controls hide(z269=(253.6,55.6653))z269 and hide(z270=(231.2,59.366))z270..hide(z271=(212,53.1985))z271
	..controls hide(z272=(192.8,47.031))z272 and hide(z273=(159.2,40.8635))z273..hide(z274=(136.8,53.1985))z274
	..controls hide(z275=(114.4,65.5332))z275 and hide(z276=(100,68.0002))z276..cycle;

 filldraw PreeditedPath8 withcolor bluecolor;

endfig;
end
