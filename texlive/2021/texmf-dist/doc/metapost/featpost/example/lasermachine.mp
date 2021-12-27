% lasermachine.mp
% L. Nobre G.
% 2006

input featpost3Dplus2D;

prologues := 1;

beginfig(1);
  numeric devi;
  devi = 0.1;
  path hart, it;
  numeric u, abot, atop, delta, i;
  u = 1in;
  abot = 60;
  atop = 145;
  delta = 23;
  z1 = ( 0, 2u);
  z2 = ( -1u, 1.25u );
  z3 = (  1u, 1.25u );
  hart = origin{dir(abot-delta)}
      ...z3...{dir(-atop)}z1{dir(atop)}
      ...z2...{dir(-abot)}cycle;
  it = hart yscaled 0.85;
  for i=2 downto -3: 
    draw lasermachine( it, devi*u*i, 0.95 ) withcolor uniformdeviate(1)*red;
  endfor;  
endfig;
end.
