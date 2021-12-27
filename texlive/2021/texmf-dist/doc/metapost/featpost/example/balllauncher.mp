% balllauncher.mp
% L. Nobre G.
% 2014

prologues := 1;

verbatimtex \documentclass{article}\begin{document} etex

beginfig(9);
  ahlength := 8mm;
  ahangle := 24;
  numeric tlen, thig, plen, phig, llen, lhig, pray, sray, wray, tang;
  numeric u, gmar, lang, hlen, rlen, pmal, pmah, aray, alim, aste, alen;
  pen basepen, outerpen;
  picture mydash, currpic[];
  u = 2.5mm;
  basepen = pencircle scaled 0.35u;
  outerpen = pencircle scaled 1.1u;
  mydash = dashpattern( off 1.5u on 1.2u off 1.5u);
  pickup basepen;
  gmar = 0.65u;       % margem de bordadura
  lang = 45;          % angulo do lancador
  tang = 45;          % angulo da tabua (10)
  hlen = 70u;         % comprimento da linha horizontal
  tlen = 120u;        % comprimento da tabua
  thig = 6u;          % altura da tabua
  plen = 60u;         % comprimento da placa
  phig = 35u;         % altura da placa
  pmal = 20u;         % parte da placa que fica em cima da tabua
  pmah = 15u;         % parte da placa que fica abaixo da horizontal
  llen = 60u;         % comprimento da lancador (na inclinacao zero)
  lhig = 15u;         % altura do lancador (na inclinacao zero)
  pray = 2u;          % raio dos cantos arredondados da placa
  sray = 5.7u;        % raio da marca da esfera no lancador
  wray = 1u;          % raio do peso de marcacao do angulo
  rlen = 20u;         % comprimento da cordinha que segura o peso de marcacao
  alim = 45;          % angulo limite do tansferidor
  aste = 5;           % divisao da escala do transferidor
  alen = 0.9u;        % comprimento das marcas da escala do transferidor
  aray = lhig-4*gmar; % raio do transferidor
  path t, p, a, lo, li, sp;
  li = (aray+3*gmar-llen,-gmar)--(aray+gmar,-gmar)--
       (aray+gmar,lhig-3*gmar)--(aray+3*gmar-llen,lhig-3*gmar)--cycle;
  lo = (aray+2*gmar-llen,-2*gmar)--(aray+2*gmar,-2*gmar)--
       (aray+2*gmar,lhig-2*gmar)--(aray+2*gmar-llen,lhig-2*gmar)--cycle;
  t = origin--(tlen,0)--(tlen,-thig)--(0,-thig)--cycle;
  p = origin--(pmal,0)---(pmal,phig-pmah-pray)...(pmal-pray,phig-pmah)---
      (pmal-plen+pray,phig-pmah)...(pmal-plen,phig-pmah-pray)---
      (pmal-plen,pray-pmah)...(pmal-plen+pray,-pmah)---
      (-pray,-pmah)...(0,pray-pmah)---cycle;
  a = buildcycle( fullcircle rotated (-10), dir(90+alim)--origin--right );
  numeric auxa, auxb, aucc;
  auxa = aray+3*gmar-llen+sray;
  auxb = 0.5*lhig-2*gmar;
  z0 = (auxa,auxb);
  sp = fullcircle scaled (2*sray) shifted z0;
  unfill lo;
  draw lo;
  fill li;
  unfill sp;
  unfill a scaled (2*aray);
  numeric i;
  for i=aste step aste until alim-aste+90:
    draw (aray-gmar)*dir(i)--(aray-gmar-alen)*dir(i);
  endfor;
  z1 = rlen*dir(90-lang+tang);
  draw origin--z1;
  fill fullcircle scaled (2*wray) shifted z1;
  currpic0 = currentpicture shifted -z0 rotated 180;
  pair haxis, vaxis, pos;
  numeric velh, velv, halfaccel, maxtime;
  halfaccel = -1.866;                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  maxtime = 20;                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  haxis = dir(180+tang-lang);
  vaxis = dir(270+tang-lang);
  velh = tlen/maxtime;
  velv = velh*sind(lang-tang)/cosd(lang-tang);
  for i=1 upto maxtime:
    pos := haxis*velh*i+vaxis*(velv+halfaccel*i)*i;
    fill fullcircle scaled (2*sray) shifted (z0+pos);
  endfor;
  drawdblarrow z0--(z0+pos) withcolor white withpen outerpen;
  drawdblarrow z0--(z0+pos) withcolor blue withpen basepen;
  picture tmppicture;
  tmppicture = currentpicture;
  currentpicture := nullpicture;
  draw t;
  draw p;
  z999 = (-0.5*(plen-pmal),sray)-z0;
  draw tmppicture rotatedaround (z0,180+lang) shifted (z999);
  undraw origin--(-hlen,0) dashed mydash withpen outerpen;
  draw origin--(-hlen,0) dashed mydash withpen basepen;
  tmppicture := currentpicture;
  currentpicture := nullpicture;
  draw tmppicture rotated (-tang);
endfig;

beginfig(0);
  draw currpic0;
  for i=0 step 90 until 270:
    draw origin--(u*dir(i));
  endfor;
  label(btex Launch etex,(0,0.5sray));
  label(btex \begin{tabular}{c} Position \\ of Ball \end{tabular} etex,
      (0,-0.5sray));
  basepen := pencircle scaled 0.35u;
  outerpen := pencircle scaled 1.1u;
  ahlength := 4mm;
  drawdblarrow origin--(gmar+sray,0) withcolor white withpen outerpen;
  drawdblarrow origin--(gmar+sray,0) withcolor blue withpen basepen;
endfig;

verbatimtex \end{document} etex

end.