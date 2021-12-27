% casadopessoal.mp
% L. Nobre G.
% 2012

prologues := 1;

numeric u, raio, cc, cp, espa, cpc;
u = 50mm;
raio = 0.4u;
espa = u-2*raio;
cpc = u;
cc = 2.16u;
cp = 1.6u;

color pgreen, cblue, fgreen;
pgreen = (139,194,52)/255;
cblue = (12,84,130)/255;
fgreen = (150,200,50)/255;
picture C, P;

def boxpath( expr w, h, r ) =
  ((w,h-r)...(w-r,h)---(r-w,h)...(-w,h-r)---(-w,r-h)...(r-w,-h)---(w-r,-h)...(w,r-h)---cycle)    
enddef;

def boxdraw( expr W, H, R, coox, cooy, col ) =
    fill boxpath( W/2, H/2, R ) shifted (coox,cooy) withcolor col
enddef;

beginfig(2);
  numeric altur;
  altur = 2*(cp+espa)+raio;
  boxdraw(u,cp,raio,-0.5*(cc+espa),0.5*(cp-altur),pgreen);
  boxdraw(cc,u,raio,0,cp+0.5*(espa-altur),pgreen);
  boxdraw(u,cp,raio,-0.5*(cc+espa),0.5*(altur-cp)-espa-raio,pgreen);
  boxdraw(cc,u,raio,0,0.5*(altur-u),pgreen);
  boxdraw(u,cp,raio,0.5*(cc+espa),0.5*(altur-cp)-espa-raio,pgreen);
  P = currentpicture;
endfig;

beginfig(1);
  boxdraw(u,altur-2*(raio+espa),raio,-0.5*(cc+espa),0,cblue);
  boxdraw(cc,u,raio,0,0.5*(u-altur),cblue);
  boxdraw(u,cpc,raio,0.5*(cc+espa),0.5*(altur-cpc)-espa-raio,cblue);
  boxdraw(cc,u,raio,0,0.5*(altur-u),cblue);
  boxdraw(u,cpc,raio,0.5*(cc+espa),0.5*(cpc-altur)+espa+raio,cblue);
  C = currentpicture;
endfig;

beginfig(3);
  draw C shifted (-raio-0.5*(cc+3*espa),0);
  draw P shifted (raio+0.5*(cc+3*espa),0);
endfig;

beginfig(4);
  numeric altur;
  altur = 2*(cp+espa)+raio;
  boxdraw(u,cp,raio,-0.5*(cc+espa),0.5*(cp-altur),fgreen);
  boxdraw(cp,u,raio,0.5*(cp-cc),cp+0.5*(espa-altur),fgreen);
  boxdraw(u,cp,raio,-0.5*(cc+espa),0.5*(altur-cp)-espa-raio,fgreen);
  boxdraw(cc,u,raio,0,0.5*(altur-u),fgreen);
endfig;

end.