% xraycamera.mp 
% L. Nobre G.
% IYP (2005)

input featpost3Dplus2D;

verbatimtex
\documentclass{article} 
\usepackage{beton} 
\usepackage{concmath} 
\usepackage{ccfonts} 
\begin{document} 
etex

beginfig(1);
  f := (10,7,8);
  Spread := 45;
  numeric w, h, l, d, t, tr, tl, ma, mb, we, ar;
  w = 3;
  h = 2;
  l = 4;
  d = 0.75;
  t = 1.2;
  tr= 0.15;
  tl= 0.25;
  ma= 0.3;
  mb= 0.2;
  ar= 2;
  we= 0.6*(h-t);
  pen thin, thick, aver;
  thin = pencircle scaled 0.5pt;
  thick= pencircle scaled 1.5pt;
  aver = pencircle scaled 1.0pt;
  ahlength:=2mm;
  ahangle:=30;
  color p[];
  p1 = (0,w/2,-ma-mb);
  p2 = (0,w/2,-ma);
  p3 = (0,w/2,0);
  p4 = (0,w/2,h);
  p5 = (0,-w/2,-ma-mb);
  p6 = (0,-w/2,-ma);
  p7 = (0,-w/2,0);
  p8 = (0,-w/2,h);
  p9 = (0,-w/2-ma-mb,0);
  p10= (0,-w/2-ma,0);
  p11= (0,-w/2-ma-mb,h);
  p12= (0,-w/2-ma,h);
  p13= (-l,-w/2-ma-mb,h);
  p14= (-l,-w/2-ma,h);
  p15= (-l,-w/2,h);
  p16= (-l,-w/2,0);
  p17= (-l,w/2,0);
  p18= (-l,w/2,h);
  p19= (-l+d,0,0);
  p20= (-l+d,0,t);
  p21= (-l-ma,0,t);
  p22= (-l-ma-mb,0,t);
  p23= (-l-ma,0,0);
  p24= (-l-ma-mb,0,0);
  p25= (-l+d,w/2+ma,0);
  p26= (-l+d,w/2+ma+mb,0);
  p27= (-l,w/2+ma,0);
  p28= (-l,w/2+ma+mb,0);
  p29= 0.5[p2,p6];
  p30= 0.5[p10,p12];
  p31= 0.5[p12,p14];
  p32= 0.5[p21,p23];
  p33= 0.5[p25,p27];
  p34= (0,0,h);
  p35= (0,0,2*t-h);
  drawoptions( withpen thick );
  draw rp(p8)--rp(p4)--rp(p3)--rp(p7)--rp(p8)--rp(p15)--rp(p16);
  draw rp(p15)--rp(p18)--rp(p17)--rp(p3)--rp(p4)--rp(p18);
  path pat[];
  pat1 = rp(p8)--rp(p4);
  pat2 = rp(p7)--rp(p16)--rp(p17);
  pat3 = rp(p4)--rp(p18);
  pat4 = pat2 cutbefore pat1 cutafter pat3;
  draw pat4;
  drawoptions( withpen thin );
  drawdblarrow rp(p2)--rp(p6);
  drawdblarrow rp(p10)--rp(p12);
  drawdblarrow rp(p14)--rp(p12);
  drawdblarrow rp(p21)--rp(p23);
  drawdblarrow rp(p25)--rp(p27);
  draw rp(p1)--rp(p3);
  draw rp(p5)--rp(p7);
  draw rp(p7)--rp(p9);
  draw rp(p8)--rp(p11);
  draw rp(p13)--rp(p15);
  draw rp(p20)--rp(p22);
  draw rp(p19)--rp(p24);
  draw rp(p17)--rp(p28);
  draw rp(p19)--rp(p26);
  drawoptions( withpen thick );
  rigorousdisc( 0, true, p20, tr, (-tl,0,0) );
  draw rp(p19)--rp(p20) withpen thin;
  drawoptions( withpen aver );
  draw rp(p34)--rp(p20)--rp(p35);
  fill ellipticpath( (0,0,t), (h-t)*blue, we*green ) withcolor 0.5white;
  draw pat2 dashed evenly withpen thick;
  draw rp(p20)--rp(p35)--rp(p34) dashed evenly;
  draw rp((0,0,t)-we*green)--rp((0,0,t)+we*green) dashed evenly;
  angline( p34, p35, p20, ar, btex $\alpha$ etex, llft );
  label.bot(btex $w$ etex, rp(p29));
  label.lft(btex $h$ etex, rp(p30));
  label.top(btex $l$ etex, rp(p31));
  label.lft(btex $t$ etex, rp(p32));
  label.rt(btex $d$ etex, rp(p33));
endfig;

verbatimtex 
\end{document} 
etex

end;
