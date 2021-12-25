% photoreverse.mp
% L. Nobre G. 
% 2004

input featpost3Dplus2D;

% bash: export TEX=latex
% tcsh: setenv TEX latex

% Converted from PostScript(TM) to MetaPost by pstoedit
% MetaPost backend contributed by Scott Pakin <pakin@uiuc.edu>
% pstoedit is Copyright (C) 1993 - 2003 Wolfgang Glunz <wglunz@pstoedit.net>

% MetaPost source tuned by PstoeditMpostPreeditor (lex scanner)
% and edited by L. Nobre G. <nobre@lince.cii.fc.ul.pt>

% Generate structured PostScript
prologues := 1;

% Display a given string with its *baseline* at a given location
% and with a given rotation angle
vardef showtext(expr origin)(expr angle)(expr string) =
  draw string infont defaultfont scaled defaultscale
    rotated angle shifted origin;
enddef;

beginfig(1);
  path PreeditedPath[];
  linecap := butt;
  linejoin := mitered;
  PreeditedPath1 =  hide(z1=(267.637,309.195))z1..controls hide(z2=(265.005,308.373))z2 and hide(z3=(262.069,307.351))z3..hide(z4=(261.056,306.906))z4
  ..controls hide(z5=(260.051,306.463))z5 and hide(z6=(257.254,305.728))z6..hide(z7=(254.841,305.272))z7
  ..controls hide(z8=(250.692,304.488))z8 and hide(z9=(204.457,292.788))z9..hide(z10=(182.447,286.951))z10
  ..controls hide(z11=(176.816,285.458))z11 and hide(z12=(168.096,283.167))z12..hide(z13=(163.069,281.859))z13
  ..controls hide(z14=(151.222,278.777))z14 and hide(z15=(141.515,276.205))z15..hide(z16=(134.184,274.205))z16
  ..controls hide(z17=(130.967,273.328))z17 and hide(z18=(125.044,271.795))z18..hide(z19=(121.022,270.8))z19
  ..controls hide(z20=(112.027,268.574))z20 and hide(z21=(85.2788,261.351))z21..hide(z22=(82.2654,260.334))z22
  ..controls hide(z23=(81.0576,259.926))z23 and hide(z24=(74.8799,258.292))z24..hide(z25=(68.5181,256.697))z25
  ..controls hide(z26=(53.8198,253.011))z26 and hide(z27=(54.844,253.994))z27..hide(z28=(56.676,245.334))z28
  ..controls hide(z29=(57.3293,242.246))z29 and hide(z30=(56.3792,241.622))z30..hide(z31=(47.3118,239.178))z31
  ..controls hide(z32=(35.9773,236.124))z32 and hide(z33=(35.1968,235.536))z33..hide(z34=(38.4636,232.521))z34
  ..controls hide(z35=(40.325,230.802))z35 and hide(z36=(46.4053,219.375))z36..hide(z37=(48.1816,214.256))z37
  ..controls hide(z38=(50.9392,206.31))z38 and hide(z39=(75.4822,148.371))z39..hide(z40=(81.2634,136.159))z40
  ..controls hide(z41=(83.6174,131.186))z41 and hide(z42=(85.5562,126.827))z42..hide(z43=(85.5562,126.506))z43
  ..controls hide(z44=(85.5562,124.811))z44 and hide(z45=(87.6768,124.222))z45..hide(z46=(90.6018,125.105))z46
  ..controls hide(z47=(94.1848,126.187))z47 and hide(z48=(94.4399,127.076))z48..hide(z49=(91.2598,127.398))z49
  ..controls hide(z50=(86.803,127.85))z50 and hide(z51=(81.7302,144.202))z51..hide(z52=(66.9248,205.847))z52
  ..controls hide(z53=(57.2651,246.065))z53 and hide(z54=(57.3979,243.248))z54..hide(z55=(65.0811,244.906))z55
  ..controls hide(z56=(67.3086,245.386))z56 and hide(z57=(71.2236,246.393))z57..hide(z58=(73.8562,247.162))z58
  ..controls hide(z59=(76.4578,247.923))z59 and hide(z60=(81.7537,249.412))z60..hide(z61=(85.5562,250.453))z61
  ..controls hide(z62=(89.3586,251.494))z62 and hide(z63=(93.6614,252.702))z63..hide(z64=(95.0623,253.122))z64
  ..controls hide(z65=(96.4841,253.549))z65 and hide(z66=(99.0842,254.15))z66..hide(z67=(100.912,254.475))z67
  ..controls hide(z68=(102.711,254.795))z68 and hide(z69=(106.689,255.777))z69..hide(z70=(109.687,256.642))z70
  ..controls hide(z71=(115.705,258.378))z71 and hide(z72=(115.391,258.29))z72..hide(z73=(129.431,262.147))z73
  ..controls hide(z74=(134.842,263.633))z74 and hide(z75=(139.964,265.139))z75..hide(z76=(140.765,265.48))z76
  ..controls hide(z77=(141.57,265.822))z77 and hide(z78=(145.519,266.957))z78..hide(z79=(149.541,268.003))z79
  ..controls hide(z80=(153.562,269.049))z80 and hide(z81=(158.005,270.236))z81..hide(z82=(159.412,270.641))z82
  ..controls hide(z83=(164.668,272.153))z83 and hide(z84=(171.332,273.79))z84..hide(z85=(180.984,275.938))z85
  ..controls hide(z86=(185.445,276.931))z86 and hide(z87=(190.982,278.282))z87..hide(z88=(193.416,278.972))z88
  ..controls hide(z89=(195.835,279.657))z89 and hide(z90=(203.214,281.625))z90..hide(z91=(209.869,283.359))z91
  ..controls hide(z92=(216.523,285.094))z92 and hide(z93=(224.23,287.209))z93..hide(z94=(227.053,288.076))z94
  ..controls hide(z95=(229.905,288.952))z95 and hide(z96=(233.467,289.972))z96..hide(z97=(235.097,290.379))z97
  ..controls hide(z98=(236.722,290.785))z98 and hide(z99=(239.484,291.549))z99..hide(z100=(241.312,292.098))z100
  ..controls hide(z101=(243.136,292.645))z101 and hide(z102=(247.309,293.664))z102..hide(z103=(250.672,294.384))z103
  ..controls hide(z104=(255.718,295.463))z104 and hide(z105=(258.712,296.887))z105..hide(z106=(255.937,296.887))z106
  ..controls hide(z107=(255.547,296.887))z107 and hide(z108=(255.206,297.399))z108..hide(z109=(255.206,297.984))z109
  ..controls hide(z110=(255.206,298.587))z110 and hide(z111=(255.864,299.081))z111..hide(z112=(256.669,299.081))z112
  ..controls hide(z113=(257.546,299.081))z113 and hide(z114=(258.131,298.613))z114..hide(z115=(258.131,297.911))z115
  ..controls hide(z116=(258.131,296.975))z116 and hide(z117=(258.643,296.837))z117..hide(z118=(260.69,297.219))z118
  ..controls hide(z119=(265.886,298.191))z119 and hide(z120=(265.944,297.765))z120..hide(z121=(265.075,264.932))z121
  ..controls hide(z122=(264.736,252.135))z122 and hide(z123=(264.811,245.7))z123..hide(z124=(265.297,245.7))z124
  ..controls hide(z125=(266.883,245.7))z125 and hide(z126=(267.399,252.793))z126..hide(z127=(267.917,281.75))z127
  ..controls hide(z128=(268.149,294.694))z128 and hide(z129=(268.149,294.694))z129..hide(z130=(265.364,299.227))z130
  ..controls hide(z131=(262.174,304.419))z131 and hide(z132=(261.87,306.04))z132..hide(z133=(263.981,306.6))z133
  ..controls hide(z134=(268.294,307.744))z134 and hide(z135=(275.854,303.103))z135..hide(z136=(280.281,296.595))z136
  ..controls hide(z137=(283.962,291.184))z137 and hide(z138=(284.721,290.268))z138..hide(z139=(289.063,286.003))z139
  ..controls hide(z140=(291.403,283.705))z140 and hide(z141=(292.5,282.033))z141..hide(z142=(292.5,280.763))z142
  ..controls hide(z143=(292.5,279.62))z143 and hide(z144=(293.377,278.063))z144..hide(z145=(294.767,276.741))z145
  ..controls hide(z146=(296.017,275.551))z146 and hide(z147=(298.525,272.208))z147..hide(z148=(300.361,269.283))z148
  ..controls hide(z149=(302.197,266.358))z149 and hide(z150=(303.934,263.981))z150..hide(z151=(304.236,263.981))z151
  ..controls hide(z152=(304.545,263.981))z152 and hide(z153=(305.201,262.957))z153..hide(z154=(305.736,261.641))z154
  ..controls hide(z155=(306.263,260.343))z155 and hide(z156=(308.249,257.363))z156..hide(z157=(310.196,254.95))z157
  ..controls hide(z158=(312.143,252.537))z158 and hide(z159=(313.706,250.225))z159..hide(z160=(313.706,249.758))z160
  ..controls hide(z161=(313.706,249.294))z161 and hide(z162=(314.657,248))z162..hide(z163=(315.827,246.87))z163
  ..controls hide(z164=(316.997,245.74))z164 and hide(z165=(318.128,243.888))z165..hide(z166=(318.35,242.738))z166
  ..controls hide(z167=(318.613,241.372))z167 and hide(z168=(319.52,240.264))z168..hide(z169=(320.982,239.521))z169
  ..controls hide(z170=(322.225,238.889))z170 and hide(z171=(323.698,237.449))z171..hide(z172=(324.309,236.267))z172
  ..controls hide(z173=(324.914,235.097))z173 and hide(z174=(325.652,234.031))z174..hide(z175=(325.955,233.89))z175
  ..controls hide(z176=(326.254,233.75))z176 and hide(z177=(327.334,231.989))z177..hide(z178=(328.331,230.015))z178
  ..controls hide(z179=(329.365,227.967))z179 and hide(z180=(331.096,225.887))z180..hide(z181=(332.316,225.225))z181
  ..controls hide(z182=(333.504,224.581))z182 and hide(z183=(335.098,222.921))z183..hide(z184=(335.863,221.532))z184
  ..controls hide(z185=(336.635,220.13))z185 and hide(z186=(338.471,217.547))z186..hide(z187=(339.995,215.719))z187
  ..controls hide(z188=(341.518,213.89))z188 and hide(z189=(343.403,211.128))z189..hide(z190=(344.236,209.503))z190
  ..controls hide(z191=(345.069,207.878))z191 and hide(z192=(347.734,204.531))z192..hide(z193=(350.232,201.971))z193
  ..controls hide(z194=(352.945,199.192))z194 and hide(z195=(354.903,196.494))z195..hide(z196=(355.168,195.17))z196
  ..controls hide(z197=(355.411,193.957))z197 and hide(z198=(356.943,191.441))z198..hide(z199=(358.596,189.54))z199
  ..controls hide(z200=(361.841,185.811))z200 and hide(z201=(364.128,182.683))z201..hide(z202=(365.551,180.034))z202
  ..controls hide(z203=(366.846,177.621))z203 and hide(z204=(365.23,174.037))z204..hide(z205=(362.846,174.037))z205
  ..controls hide(z206=(361.091,174.037))z206 and hide(z207=(358.527,172.096))z207..hide(z208=(359.273,171.332))z208
  ..controls hide(z209=(359.507,171.093))z209 and hide(z210=(360.799,171.399))z210..hide(z211=(362.115,172.005))z211
  ..controls hide(z212=(368.257,174.837))z212 and hide(z213=(375.377,173.789))z213..hide(z214=(371.109,170.681))z214
  ..controls hide(z215=(370.133,169.97))z215 and hide(z216=(369.915,169.138))z216..hide(z217=(370.268,167.456))z217
  ..controls hide(z218=(370.58,165.969))z218 and hide(z219=(370.316,164.385))z219..hide(z220=(369.501,162.849))z220
  ..controls hide(z221=(368.801,161.533))z221 and hide(z222=(368.321,159.646))z222..hide(z223=(368.404,158.535))z223
  ..controls hide(z224=(368.485,157.438))z224 and hide(z225=(368.095,156.301))z225..hide(z226=(367.526,155.976))z226
  ..controls hide(z227=(366.853,155.591))z227 and hide(z228=(366.377,153.635))z228..hide(z229=(366.137,150.272))z229
  ..controls hide(z230=(365.933,147.42))z230 and hide(z231=(365.431,143.858))z231..hide(z232=(365.003,142.228))z232
  ..controls hide(z233=(364.584,140.633))z233 and hide(z234=(363.853,137.475))z234..hide(z235=(363.395,135.281))z235
  ..controls hide(z236=(362.936,133.087))z236 and hide(z237=(362.088,130.162))z237..hide(z238=(361.53,128.846))z238
  ..controls hide(z239=(360.967,127.519))z239 and hide(z240=(360.506,125.117))z240..hide(z241=(360.506,123.508))z241
  ..controls hide(z242=(360.506,121.899))z242 and hide(z243=(360.111,119.859))z243..hide(z244=(359.629,118.974))z244
  ..controls hide(z245=(359.15,118.097))z245 and hide(z246=(358.834,116.269))z246..hide(z247=(358.934,114.952))z247
  ..controls hide(z248=(359.034,113.636))z248 and hide(z249=(358.602,111.825))z249..hide(z250=(357.983,110.967))z250
  ..controls hide(z251=(357.35,110.088))z251 and hide(z252=(356.85,108.298))z252..hide(z253=(356.85,106.909))z253
  ..controls hide(z254=(356.85,105.519))z254 and hide(z255=(356.382,103.567))z255..hide(z256=(355.79,102.485))z256
  ..controls hide(z257=(355.209,101.424))z257 and hide(z258=(354.401,97.4023))z258..hide(z259=(354.005,93.5999))z259
  ..controls hide(z260=(352.968,83.6548))z260 and hide(z261=(352.529,81.2454))z261..hide(z262=(351.429,79.4868))z262
  ..controls hide(z263=(350.069,77.3098))z263 and hide(z264=(348.615,69.103))z264..hide(z265=(349.163,66.6899))z265
  ..controls hide(z266=(349.581,64.8442))z266 and hide(z267=(349.391,64.6501))z267..hide(z268=(346.247,63.7119))z268
  ..controls hide(z269=(341.811,62.3882))z269 and hide(z270=(340.625,63.3262))z270..hide(z271=(340.914,67.9331))z271
  ..controls hide(z272=(341.263,73.5166))z272 and hide(z273=(339.121,72.5398))z273..hide(z274=(338.347,66.7629))z274
  ..controls hide(z275=(337.642,61.498))z275 and hide(z276=(335.318,60.0554))z276..hide(z277=(334.459,64.3499))z277
  ..controls hide(z278=(334.225,65.5198))z278 and hide(z279=(333.412,67.0278))z279..hide(z280=(332.646,67.7136))z280
  ..controls hide(z281=(331.914,68.3679))z281 and hide(z282=(331.256,69.5457))z282..hide(z283=(331.256,70.2))z283
  ..controls hide(z284=(331.256,70.8691))z284 and hide(z285=(330.353,72.7593))z285..hide(z286=(329.3,74.2949))z286
  ..controls hide(z287=(325.789,79.4136))z287 and hide(z288=(320.082,90.9722))z288..hide(z289=(318.72,95.7205))z289
  ..controls hide(z290=(317.42,100.254))z290 and hide(z291=(314.317,105.3))z291..hide(z292=(312.829,105.3))z292
  ..controls hide(z293=(308.461,105.3))z293 and hide(z294=(294.694,129.504))z294..hide(z295=(294.694,137.182))z295
  ..controls hide(z296=(294.694,141.643))z296 and hide(z297=(291.904,145.559))z297..hide(z298=(290.891,142.521))z298
  ..controls hide(z299=(290.542,141.474))z299 and hide(z300=(287.381,143.252))z300..hide(z301=(287.381,144.495))z301
  ..controls hide(z302=(287.381,144.918))z302 and hide(z303=(286.723,145.888))z303..hide(z304=(285.992,146.542))z304
  ..controls hide(z305=(285.21,147.242))z305 and hide(z306=(284.219,148.809))z306..hide(z307=(283.725,150.125))z307
  ..controls hide(z308=(283.246,151.402))z308 and hide(z309=(282.133,154.294))z309..hide(z310=(281.275,156.487))z310
  ..controls hide(z311=(280.406,158.711))z311 and hide(z312=(278.15,162.52))z312..hide(z313=(276.229,165.006))z313
  ..controls hide(z314=(274.309,167.493))z314 and hide(z315=(272.756,170.223))z315..hide(z316=(272.756,171.112))z316
  ..controls hide(z317=(272.756,172))z317 and hide(z318=(271.953,173.672))z318..hide(z319=(270.965,174.842))z319
  ..controls hide(z320=(269.404,176.689))z320 and hide(z321=(269.276,177.255))z321..hide(z322=(269.969,179.229))z322
  ..controls hide(z323=(270.893,181.862))z323 and hide(z324=(269.866,185.685))z324..hide(z325=(267.573,188.15))z325
  ..controls hide(z326=(266.139,189.692))z326 and hide(z327=(265.841,220.618))z327..hide(z328=(267.229,223.909))z328
  ..controls hide(z329=(267.537,224.64))z329 and hide(z330=(267.426,225.664))z330..hide(z331=(266.979,226.212))z331
  ..controls hide(z332=(266.536,226.755))z332 and hide(z333=(265.939,230.782))z333..hide(z334=(265.65,235.17))z334
  ..controls hide(z335=(264.577,251.473))z335 and hide(z336=(263.579,226.029))z336..hide(z337=(264.07,194.878))z337
  ..controls hide(z338=(264.101,192.884))z338 and hide(z339=(263.729,188.479))z339..hide(z340=(263.25,185.189))z340
  ..controls hide(z341=(262.42,179.485))z341 and hide(z342=(262.452,179.086))z342..hide(z343=(263.835,177.876))z343
  ..controls hide(z344=(264.65,177.164))z344 and hide(z345=(266.512,174.33))z345..hide(z346=(268.04,171.478))z346
  ..controls hide(z347=(269.567,168.626))z347 and hide(z348=(271.381,165.501))z348..hide(z349=(272.135,164.421))z349
  ..controls hide(z350=(272.9,163.324))z350 and hide(z351=(273.487,162.293))z351..hide(z352=(273.487,162.045))z352
  ..controls hide(z353=(273.487,161.804))z353 and hide(z354=(275.511,158.242))z354..hide(z355=(277.976,154.147))z355
  ..controls hide(z356=(282.728,146.25))z356 and hide(z357=(284.237,143.655))z357..hide(z358=(287.245,138.206))z358
  ..controls hide(z359=(289.99,133.234))z359 and hide(z360=(299.06,117.685))z360..hide(z361=(301.023,114.587))z361
  ..controls hide(z362=(301.95,113.124))z362 and hide(z363=(302.737,111.747))z363..hide(z364=(302.737,111.589))z364
  ..controls hide(z365=(302.737,111.183))z365 and hide(z366=(314.213,91.2598))z366..hide(z367=(321.733,78.6091))z367
  ..controls hide(z368=(325.081,72.9785))z368 and hide(z369=(329.837,64.9495))z369..hide(z370=(332.304,60.7668))z370
  ..controls hide(z371=(338.083,50.968))z371 and hide(z372=(337.658,51.1626))z372..hide(z373=(346.101,54.4392))z373
  ..controls hide(z374=(354.37,57.6489))z374 and hide(z375=(353.838,56.9644))z375..hide(z376=(355.457,66.4705))z376
  ..controls hide(z377=(356.179,70.7117))z377 and hide(z378=(357.29,77.0364))z378..hide(z379=(357.91,80.4373))z379
  ..controls hide(z380=(358.537,83.8743))z380 and hide(z381=(359.518,89.2683))z381..hide(z382=(360.104,92.5029))z382
  ..controls hide(z383=(360.686,95.7136))z383 and hide(z384=(362.512,105.446))z384..hide(z385=(364.153,114.075))z385
  ..controls hide(z386=(367.6,132.21))z386 and hide(z387=(370.829,149.516))z387..hide(z388=(373.021,161.606))z388
  ..controls hide(z389=(373.857,166.213))z389 and hide(z390=(374.711,170.875))z390..hide(z391=(374.912,171.917))z391
  ..controls hide(z392=(375.191,173.367))z392 and hide(z393=(374.51,174.842))z393..hide(z394=(372.046,178.132))z394
  ..controls hide(z395=(364.981,187.565))z395 and hide(z396=(360.559,193.664))z396..hide(z397=(357.253,198.534))z397
  ..controls hide(z398=(355.351,201.336))z398 and hide(z399=(351.495,206.797))z399..hide(z400=(348.733,210.6))z400
  ..controls hide(z401=(345.967,214.408))z401 and hide(z402=(340.482,222.007))z402..hide(z403=(336.594,227.419))z403
  ..controls hide(z404=(332.706,232.83))z404 and hide(z405=(327.32,240.264))z405..hide(z406=(324.675,243.872))z406
  ..controls hide(z407=(321.994,247.528))z407 and hide(z408=(318.404,252.436))z408..hide(z409=(316.598,254.914))z409
  ..controls hide(z410=(307.005,268.076))z410 and hide(z411=(302.032,274.894))z411..hide(z412=(295.05,284.456))z412
  ..controls hide(z413=(290.779,290.306))z413 and hide(z414=(285.377,297.845))z414..hide(z415=(283.001,301.275))z415
  ..controls hide(z416=(275.704,311.805))z416 and hide(z417=(275.811,311.75))z417..cycle;
  PreeditedPath2 =  hide(z418=(254.987,177.913))z418..controls hide(z419=(254.744,177.67))z419 and hide(z420=(253.305,177.229))z420..hide(z421=(251.769,176.926))z421
  ..controls hide(z422=(250.239,176.624))z422 and hide(z423=(246.065,175.346))z423..hide(z424=(242.482,174.084))z424
  ..controls hide(z425=(235.17,171.507))z425 and hide(z426=(231.263,171.562))z426..hide(z427=(228.296,174.285))z427
  ..controls hide(z428=(226.59,175.85))z428 and hide(z429=(225.265,175.207))z429..hide(z430=(225.222,172.794))z430
  ..controls hide(z431=(225.181,170.533))z431 and hide(z432=(221.496,168.977))z432..hide(z433=(216.377,169.058))z433
  ..controls hide(z434=(212.647,169.117))z434 and hide(z435=(210.978,168.823))z435..hide(z436=(209.978,167.931))z436
  ..controls hide(z437=(209.247,167.279))z437 and hide(z438=(208.138,166.725))z438..hide(z439=(207.565,166.725))z439
  ..controls hide(z440=(206.982,166.725))z440 and hide(z441=(205.847,166.11))z441..hide(z442=(205.042,165.359))z442
  ..controls hide(z443=(202.136,162.644))z443 and hide(z444=(198.973,161.521))z444..hide(z445=(194.147,161.49))z445
  ..controls hide(z446=(190.929,161.469))z446 and hide(z447=(188.598,161.014))z447..hide(z448=(186.688,160.034))z448
  ..controls hide(z449=(185.149,159.244))z449 and hide(z450=(181.789,157.939))z450..hide(z451=(179.156,157.109))z451
  ..controls hide(z452=(176.524,156.279))z452 and hide(z453=(173.783,155.284))z453..hide(z454=(173.014,154.879))z454
  ..controls hide(z455=(172.251,154.477))z455 and hide(z456=(169.942,153.653))z456..hide(z457=(167.895,153.051))z457
  ..controls hide(z458=(165.847,152.448))z458 and hide(z459=(163.522,151.492))z459..hide(z460=(162.74,150.93))z460
  ..controls hide(z461=(161.973,150.38))z461 and hide(z462=(160.18,149.906))z462..hide(z463=(158.864,149.906))z463
  ..controls hide(z464=(157.548,149.906))z464 and hide(z465=(155.467,149.416))z465..hide(z466=(154.367,148.846))z466
  ..controls hide(z467=(153.27,148.278))z467 and hide(z468=(152.149,147.955))z468..hide(z469=(151.954,148.151))z469
  ..controls hide(z470=(151.753,148.351))z470 and hide(z471=(150.016,147.838))z471..hide(z472=(148.115,147.018))z472
  ..controls hide(z473=(146.213,146.197))z473 and hide(z474=(143.748,145.515))z474..hide(z475=(142.667,145.51))z475
  ..controls hide(z476=(128.554,145.452))z476 and hide(z477=(126.589,144.869))z477..hide(z478=(124.219,140.034))z478
  ..controls hide(z479=(123.528,138.624))z479 and hide(z480=(122.192,137.646))z480..hide(z481=(119.925,136.89))z481
  ..controls hide(z482=(118.097,136.281))z482 and hide(z483=(115.887,135.452))z483..hide(z484=(114.952,135.025))z484
  ..controls hide(z485=(114.015,134.597))z485 and hide(z486=(112.283,134.322))z486..hide(z487=(111.04,134.404))z487
  ..controls hide(z488=(109.836,134.483))z488 and hide(z489=(107.75,133.91))z489..hide(z490=(106.506,133.16))z490
  ..controls hide(z491=(105.256,132.406))z491 and hide(z492=(103.033,131.403))z492..hide(z493=(101.644,130.967))z493
  ..controls hide(z494=(100.254,130.53))z494 and hide(z495=(98.2485,129.46))z495..hide(z496=(97.2561,128.627))z496
  ..controls hide(z497=(95.428,127.091))z497 and hide(z498=(95.428,127.091))z498..hide(z499=(98.426,127.436))z499
  ..controls hide(z500=(101.129,127.747))z500 and hide(z501=(115.684,132.079))z501..hide(z502=(128.334,136.337))z502
  ..controls hide(z503=(131.163,137.289))z503 and hide(z504=(139.669,139.925))z504..hide(z505=(147.347,142.228))z505
  ..controls hide(z506=(155.025,144.531))z506 and hide(z507=(164.185,147.312))z507..hide(z508=(167.822,148.444))z508
  ..controls hide(z509=(171.478,149.581))z509 and hide(z510=(176.841,151.206))z510..hide(z511=(179.887,152.1))z511
  ..controls hide(z512=(182.911,152.987))z512 and hide(z513=(191.441,155.611))z513..hide(z514=(198.9,157.95))z514
  ..controls hide(z515=(206.359,160.288))z515 and hide(z516=(215.873,163.225))z516..hide(z517=(220.106,164.495))z517
  ..controls hide(z518=(224.347,165.767))z518 and hide(z519=(230.894,167.764))z519..hide(z520=(234.731,168.955))z520
  ..controls hide(z521=(238.534,170.136))z521 and hide(z522=(244.657,171.975))z522..hide(z523=(248.259,173.019))z523
  ..controls hide(z524=(255.864,175.222))z524 and hide(z525=(259.346,176.934))z525..hide(z526=(257.985,177.8))z526
  ..controls hide(z527=(256.888,178.498))z527 and hide(z528=(255.62,178.546))z528..cycle;
  dotlabels.rt(371,375,416,44,26);
  draw PreeditedPath1 withcolor green;
  draw PreeditedPath2 withcolor red;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEATPOST specific %%
  PhotoMarks := 5;
  PhotoPair1 := z26; PhotoPoint1 := (7.0,-0.6,5.6);
  PhotoPair2 := z44; PhotoPoint2 := (7.0,0,0);
  PhotoPair3 := z371; PhotoPoint3:= (0,6.4,0);
  PhotoPair4 := z375; PhotoPoint4:= (-0.6,6.4,0);
  PhotoPair5 := z416; PhotoPoint5:= (-0.6,-0.6,5.6);
  color pp[];
  pp1 = (7.0,-0.6,7.0);
  pp5 = (-0.6,-0.6,7.0);
  pp2 = (8.4,0,0);
  pp3 = (8.4,6.4,0);
  path pt[];
  f := (15,20,15);                          %%%%%%%%%%%%%%%%%%%%%%%
  viewcentr := black;                       %%
  Spread := 30;                             %%     Parameters
  ShiftV := (250,100);                      %%
% SphericalDistortion := true;              %%%%%%%%%%%%%%%%%%%%%%%

  pt1 = rp(PhotoPoint1)--rp(pp1)--rp(pp5)--rp(PhotoPoint5);
  pt2 = rp(PhotoPoint2)--rp(pp2)--rp(pp3)--rp(PhotoPoint3);
  draw pt1 withpen pencircle scaled 3pt withcolor blue;
  draw pt2 withpen pencircle scaled 3pt withcolor blue;

  photoreverse(400,5.5,2.5);               %% Changes parameters

  pt3 = rp(PhotoPoint1)--rp(pp1)--rp(pp5)--rp(PhotoPoint5);
  pt4 = rp(PhotoPoint2)--rp(pp2)--rp(pp3)--rp(PhotoPoint3);
  draw pt3 withpen pencircle scaled 1pt;
  draw pt4 withpen pencircle scaled 1pt;
endfig;
end
