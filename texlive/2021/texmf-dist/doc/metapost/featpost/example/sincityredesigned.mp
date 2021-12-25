
prologues := 1;

beginfig(1);
  linejoin := mitered;
  numeric controlparam;
  controlparam = 5mm;
  path PreeditedPath[];
  PreeditedPath1 = hide(z1=(304.883,266.453))z1..controls
  hide(z2=(303.009,265.609))z2 and
  hide(z3=(297.112,264.877))z3..hide(z4=(291.817,264.83))z4
  --hide(z5=(282.167,264.745))z5--hide(z6=(280.885,249.372))z6..controls
  hide(z7=(280.174,240.832))z7 and
  hide(z8=(279.111,234))z8..hide(z9=(278.494,234))z9
  ..controls hide(z10=(277.871,234))z10 and
  hide(z11=(274.31,236.621))z11..hide(z12=(270.381,239.971))z12
  ..controls hide(z13=(262.387,246.786))z13 and
  hide(z14=(259.108,247.209))z14..hide(z15=(245.785,243.138))z15
  ..controls hide(z16=(240.661,241.572))z16 and
  hide(z17=(234.133,240.4))z17..hide(z18=(231.096,240.5))z18
  ..controls hide(z19=(221.716,240.809))z19 and 
  hide(z20=(220.916,238.612))z20..hide(z21=(219.534,208.721))z21
  ..controls hide(z22=(218.847,193.861))z22 and
  hide(z23=(218.823,178.971))z23..hide(z24=(219.482,175.585))z24
  ..controls hide(z25=(220.592,169.876))z25 and 
  hide(z26=(221.19,169.381))z26..hide(z27=(227.851,168.665))z27
  ..controls hide(z28=(235.366,167.856))z28 and 
  hide(z29=(247.664,172.576))z29..hide(z30=(247.664,176.269))z30
  ..controls hide(z31=(247.664,180.709))z31 and 
  hide(z32=(255.027,204.963))z32..hide(z33=(256.375,204.963))z33
  ..controls hide(z34=(257.214,204.963))z34 and 
  hide(z35=(262.182,202.207))z35..hide(z36=(267.307,198.899))z36
  ..controls hide(z37=(280.629,190.298))z37 and 
  hide(z38=(287.926,188.002))z38..hide(z39=(297.368,189.437))z39
  ..controls hide(z40=(310.542,191.44))z40 and 
  hide(z41=(314.953,195.911))z41..hide(z42=(318.083,210.429))z42
  ..controls hide(z43=(323.752,236.733))z43 and 
  hide(z44=(324.666,242.082))z44..hide(z45=(325.697,255.009))z45
  --hide(z46=(326.746,268.161))z46
  --hide(z47=(317.523,268.075))z47..controls hide(z48=(312.398,268.028))z48
  and hide(z49=(306.781,267.307))z49..cycle;
  PreeditedPath2 = hide(z50=(182.759,231.97))z50..controls
  hide(z51=(174.744,229.708))z51 and 
  hide(z52=(174.156,228.876))z52..hide(z53=(173.193,218.457))z53
  ..controls hide(z54=(172.721,213.357))z54 and 
  hide(z55=(171.267,202.743))z55..hide(z56=(169.941,194.715))z56
  ..controls hide(z57=(163.9,158.163))z57 and 
  hide(z58=(164.703,149.053))z58..hide(z59=(173.707,151.991))z59
  ..controls hide(z60=(176.323,152.845))z60 and 
  hide(z61=(182.588,153.88))z61..hide(z62=(187.542,154.276))z62
  ..controls hide(z63=(201.233,155.371))z63 and 
  hide(z64=(204.848,159.872))z64..hide(z65=(208.362,180.197))z65
  ..controls hide(z66=(209.514,186.858))z66 and 
  hide(z67=(211.473,197.07))z67..hide(z68=(212.786,203.255))z68
  ..controls hide(z69=(219.311,234))z69 and 
  hide(z70=(219.311,234))z70..hide(z71=(201.889,233.805))z71
  ..controls hide(z72=(194.715,233.725))z72 and 
  hide(z73=(186.02,232.891))z73..cycle;
  PreeditedPath3 = hide(z74=(127.248,220.173))z74..controls
  hide(z75=(107.435,213.404))z75 and 
  hide(z76=(79.562,199.048))z76..hide(z77=(74.5039,193.007))z77
  ..controls hide(z78=(68.4971,185.833))z78 and 
  hide(z79=(61.489,170.106))z79..hide(z80=(61.489,163.8))z80
  ..controls hide(z81=(61.489,159.558))z81 and 
  hide(z82=(62.5991,158.602))z82..hide(z83=(71.3101,155.345))z83
  ..controls hide(z84=(76.7759,153.302))z84 and 
  hide(z85=(85.1082,150.878))z85..hide(z86=(90.0984,149.879))z86
  ..controls hide(z87=(95.0518,148.889))z87 and 
  hide(z88=(99.0657,147.551))z88..hide(z89=(99.0657,146.89))z89
  ..controls hide(z90=(99.0657,146.226))z90 and 
  hide(z91=(94.2832,143.342))z91..hide(z92=(88.3049,140.401))z92
  ..controls hide(z93=(76.178,134.435))z93 and 
  hide(z94=(63.4753,131.962))z94..hide(z95=(58.5852,134.617))z95
  ..controls hide(z96=(56.0232,136.008))z96 and 
  hide(z97=(54.4861,135.947))z97..hide(z98=(52.4363,134.375))z98
  ..controls hide(z99=(49.1121,131.825))z99 and 
  hide(z100=(51.0422,117.683))z100..hide(z101=(55.3677,112.901))z101
  ..controls hide(z102=(59.9397,107.845))z102 and 
  hide(z103=(71.7371,105.477))z103..hide(z104=(81.4729,107.659))z104
  ..controls hide(z105=(100.944,112.023))z105 and 
  hide(z106=(124.681,122.025))z106..hide(z107=(142.791,133.498))z107
  ..controls hide(z108=(155.26,141.396))z108 and 
  hide(z109=(156.223,142.747))z109..hide(z110=(157.815,154.577))z110
  ..controls hide(z111=(159.443,166.668))z111 and 
  hide(z112=(156.455,169.328))z112..hide(z113=(135.959,174.038))z113
  ..controls hide(z114=(126.906,176.118))z114 and 
  hide(z115=(119.187,178.418))z115..hide(z116=(118.724,179.172))z116
  ..controls hide(z117=(117.77,180.728))z117 and 
  hide(z118=(147.232,198.131))z118..hide(z119=(150.819,198.131))z119
  ..controls hide(z120=(152.271,198.131))z120 and 
  hide(z121=(154.38,195.911))z121..hide(z122=(155.687,193.007))z122
  ..controls hide(z123=(156.993,190.104))z123 and 
  hide(z124=(158.437,187.883))z124..hide(z125=(159.018,187.883))z125
  ..controls hide(z126=(159.582,187.883))z126 and 
  hide(z127=(162.09,191.641))z127..hide(z128=(164.603,196.252))z128
  ..controls hide(z129=(170.373,206.842))z129 and 
  hide(z130=(170.682,217.263))z130..hide(z131=(165.337,221.03))z131
  ..controls hide(z132=(159.717,224.992))z132 and 
  hide(z133=(140.058,224.549))z133..cycle;
  PreeditedPath4 = hide(z134=(373.29,181.905))z134..controls
  hide(z135=(367.765,173.877))z135 and 
  hide(z136=(363.069,167.387))z136..hide(z137=(362.785,167.387))z137
  ..controls hide(z138=(361.491,167.387))z138 and 
  hide(z139=(348.78,175.029))z139..hide(z140=(343.826,178.785))z140
  --hide(z141=(338.361,182.93))z141
  --hide(z142=(331.529,176.098))z142..controls
  hide(z143=(327.771,172.34))z143 and 
  hide(z144=(323.39,167.267))z144..hide(z145=(321.793,164.825))z145
  ..controls hide(z146=(318.889,160.384))z146 and 
  hide(z147=(318.889,160.384))z147..hide(z148=(327.429,151.909))z148
  ..controls hide(z149=(340.46,138.978))z149 and 
  hide(z150=(341.004,137.325))z150..hide(z151=(336.781,123.49))z151
  ..controls hide(z152=(334.853,117.171))z152 and 
  hide(z153=(333.626,109.929))z153..hide(z154=(334.091,107.606))z154
  ..controls hide(z155=(334.924,103.437))z155 and 
  hide(z156=(335.286,103.338))z156..hide(z157=(349.292,103.421))z157
  ..controls hide(z158=(357.149,103.468))z158 and 
  hide(z159=(363.91,103.854))z159..hide(z160=(364.237,104.275))z160
  ..controls hide(z161=(364.569,104.701))z161 and 
  hide(z162=(368.408,113.413))z162..hide(z163=(372.856,123.832))z163
  ..controls hide(z164=(380.732,142.279))z164 and 
  hide(z165=(388.471,155.431))z165..hide(z166=(400.498,170.803))z166
  ..controls hide(z167=(403.838,175.073))z167 and 
  hide(z168=(406.511,179.628))z168..hide(z169=(406.511,181.051))z169
  ..controls hide(z170=(406.511,183.073))z170 and 
  hide(z171=(394.042,191.439))z171..hide(z172=(384.307,195.948))z172
  ..controls hide(z173=(383.74,196.211))z173 and 
  hide(z174=(378.815,189.933))z174..cycle;
  PreeditedPath5 = hide(z175=(306.591,170.416))z175..controls
  hide(z176=(284.045,163.442))z176 and 
  hide(z177=(261.654,155.676))z177..hide(z178=(250.226,150.864))z178
  ..controls hide(z179=(242.882,147.771))z179 and 
  hide(z180=(235.916,145.182))z180..hide(z181=(234.939,145.182))z181
  ..controls hide(z182=(233.842,145.182))z182 and 
  hide(z183=(232.98,137.496))z183..hide(z184=(232.719,125.369))z184
  --hide(z185=(232.292,105.556))z185
  --hide(z186=(238.697,106.666))z186..controls
  hide(z187=(242.284,107.288))z187 and 
  hide(z188=(246.53,108.19))z188..hide(z189=(248.347,108.716))z189
  ..controls hide(z190=(251.392,109.597))z190 and 
  hide(z191=(251.532,108.801))z191..hide(z192=(250.615,95.8203))z192
  ..controls hide(z193=(249.493,79.9358))z193 and 
  hide(z194=(250.483,78.3792))z194..hide(z195=(260.816,79.7881))z195
  ..controls hide(z196=(270.245,81.074))z196 and 
  hide(z197=(273.79,85.0598))z197..hide(z198=(284.336,106.239))z198
  ..controls hide(z199=(292.162,121.953))z199 and 
  hide(z200=(293.65,123.879))z200..hide(z201=(301.04,127.846))z201
  ..controls hide(z202=(308.217,131.699))z202 and 
  hide(z203=(309.449,133.226))z203..hide(z204=(311.715,141.083))z204
  ..controls hide(z205=(313.144,146.036))z205 and 
  hide(z206=(314.277,151.74))z206..hide(z207=(314.277,153.979))z207
  ..controls hide(z208=(314.277,156.237))z208 and 
  hide(z209=(314.957,160.384))z209..hide(z210=(315.831,163.458))z210
  ..controls hide(z211=(317.785,170.333))z211 and 
  hide(z212=(314.448,172.845))z212..cycle;
  PreeditedPath6 = hide(z213=(200.352,136.887))z213..controls
  hide(z214=(194.544,135.348))z214 and 
  hide(z215=(189.243,133.54))z215..hide(z216=(188.612,132.885))z216
  ..controls hide(z217=(187.443,131.669))z217 and 
  hide(z218=(185.57,95.9912))z218..hide(z219=(185.617,75.8364))z219
  ..controls hide(z220=(185.647,63.0261))z220 and 
  hide(z221=(188.351,57.2344))z221..hide(z222=(193.52,58.9067))z222
  ..controls hide(z223=(195.607,59.5823))z223 and 
  hide(z224=(199.327,60.5022))z224..hide(z225=(201.889,60.9766))z225
  ..controls hide(z226=(204.405,61.4424))z226 and 
  hide(z227=(208.467,64.3926))z227..hide(z228=(210.827,67.467))z228
  ..controls hide(z229=(218.335,77.2478))z229 and 
  hide(z230=(226.481,103.677))z230..hide(z231=(227.989,123.149))z231
  ..controls hide(z232=(228.823,133.909))z232 and 
  hide(z233=(228.631,134.752))z233..hide(z234=(224.777,137.26))z234
  ..controls hide(z235=(219.588,140.637))z235 and 
  hide(z236=(214.187,140.554))z236..cycle;
  PreeditedPath7 = hide(z237=(147.232,125.456))z237..controls
  hide(z238=(124.174,119.222))z238 and 
  hide(z239=(97.5039,98.374))z239..hide(z240=(88.1265,79.2524))z240
  --hide(z241=(83.3518,69.5166))z241
  --hide(z242=(88.6506,59.0977))z242..controls
  hide(z243=(97.6846,41.3342))z243 and 
  hide(z244=(99.9316,39.4353))z244..hide(z245=(114.267,37.449))z245
  ..controls hide(z246=(125.814,35.8491))z246 and 
  hide(z247=(128.102,36.0742))z247..hide(z248=(140.741,40.0542))z248
  ..controls hide(z249=(160.896,46.4009))z249 and 
  hide(z250=(167.846,50.8127))z250..hide(z251=(170.471,58.927))z251
  ..controls hide(z252=(171.687,62.6846))z252 and 
  hide(z253=(174.603,67.5432))z253..hide(z254=(176.952,69.7241))z254
  ..controls hide(z255=(181.734,74.165))z255 and 
  hide(z256=(181.908,75.1035))z256..hide(z257=(178.752,79.4233))z257
  ..controls hide(z258=(176.002,83.186))z258 and 
  hide(z259=(173.877,82.7844))z259..hide(z260=(152.869,74.531))z260
  ..controls hide(z261=(142.108,70.3037))z261 and 
  hide(z262=(134.704,68.4663))z262..hide(z263=(128.188,68.4065))z263
  --hide(z264=(118.879,68.321))z264
  --hide(z265=(123.902,73.616))z265..controls
  hide(z266=(135.407,85.7429))z266 and 
  hide(z267=(158.517,97.1145))z267..hide(z268=(164.654,93.6689))z268
  ..controls hide(z269=(168.903,91.2834))z269 and 
  hide(z270=(170.716,93.4292))z270..hide(z271=(173.831,104.531))z271
  ..controls hide(z272=(178.624,121.612))z272 and 
  hide(z273=(178.563,124.468))z273..hide(z274=(173.365,126.432))z274
  ..controls hide(z275=(167.23,128.75))z275 and 
  hide(z276=(158.163,128.411))z276..cycle;
  numeric i, j;
  for i=1 upto 7:
    PreeditedPath[7+i] = point 0 of PreeditedPath[i]
    for j=controlparam step controlparam until arclength PreeditedPath[i]:
      --point (arctime j of PreeditedPath[i]) of PreeditedPath[i]
    endfor --cycle;
    fill PreeditedPath[7+i];
  endfor;
  picture sincitytext;
  sincitytext = currentpicture;
  currentpicture := nullpicture;
  bboxmargin := 1cm;
  fill bbox sincitytext;
  draw sincitytext withcolor (0.99,0.05,0.05);
endfig;

beginfig(2);
  pair shif;
  shif = (87mm,61.5mm);
  color madred;
  madred = (0.99,0.05,0.05);
  linejoin := mitered;
  numeric controlparam, devparam, ranparam;
  controlparam = 0.95mm;
  devparam = 1mm;
  ranparam = 1.5mm;
  pickup pencircle scaled 2mm;
  numeric i, j, t;
  pair aux, dev, tdr;
  for i=1 upto 7:
    aux := point 0 of PreeditedPath[i];
    tdr := unitvector( direction 0 of PreeditedPath[i] );
    dev := normaldeviate*ranparam*(ypart tdr,-xpart tdr);
    dev := dev+devparam*(ypart tdr,-xpart tdr);
    PreeditedPath[14+i] = (dev+aux)
    for j=controlparam step controlparam until arclength PreeditedPath[i]:
      hide( t := arctime j of PreeditedPath[i];
	aux := point t of PreeditedPath[i];
	tdr := unitvector( direction t of PreeditedPath[i] );
	dev := normaldeviate*ranparam*(ypart tdr,-xpart tdr) ;
      	dev := dev+devparam*(ypart tdr,-xpart tdr) )
      --(dev+aux)
    endfor --cycle;
  endfor;
  for i=1 upto 3:
    fill PreeditedPath[14+i] withcolor madred;
    fill PreeditedPath[i] withcolor black;
    draw PreeditedPath[i] withcolor madred;
  endfor;
  for i=4 upto 7:
    fill PreeditedPath[14+i] shifted shif withcolor madred;
    fill PreeditedPath[i] shifted shif withcolor black;
    draw PreeditedPath[i] shifted shif withcolor madred;
  endfor;
  picture sincitytext;
  sincitytext = currentpicture rotated -14.5;
  currentpicture := nullpicture;
  bboxmargin := 1cm;
  fill bbox sincitytext;
  draw sincitytext;  
endfig;

end
