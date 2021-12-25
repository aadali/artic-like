% cylinderborders.mp
% L. Nobre G.
% 2014

prologues := 1;

input featpost3Dplus2D;

beginfig(1);
  f := (3,5,2);
  numeric margi, ray, dist, len, marj;
  margi = 0.15;
  marj = 0.37;
  ray = 0.2;
  dist = 0.6;
  len = 0.25;
  color pos, pa, pb, pc, pd, qa, qb, qc, qd, ra, P, Q;
  pos = (marj,dist,0);
  numeric pqwid, pqlen, theta;
  pqwid = (dist ++ ray) + 2*margi;
  pqlen = len + 4*margi;
  theta = 2*angle(dist,ray);
  pa = (2margi+marj,dist+margi,0);
  pb = (2margi+marj,-margi,0);
  pc = (-len-margi,-margi,0);
  pd = (-len-margi,dist+margi,0);
  qa = (margi+marj,dist+margi*cosd(theta),-margi*sind(theta));
  qb = (margi+marj,dist-(dist+margi)*cosd(theta),(dist+margi)*sind(theta));
  qc = (-len-2margi,dist-(dist+margi)*cosd(theta),(dist+margi)*sind(theta));
  qd = (-len-2margi,dist+margi*cosd(theta),-margi*sind(theta));
  ra = (margi+marj,dist,0);
  path lowplan, higplan;
  lowplan = rp(pa)--rp(pb)--rp(pc)--rp(pd)--cycle;
  higplan = rp(qa)--rp(qb)--rp(qc)--rp(qd)--cycle;
  unfill higplan;
  draw higplan;
  unfill lowplan;
  draw lowplan;
  draw rp(pos)--rp(-len*red);
  rigorousdisc( 0, true, (0,0,ray), ray, -len*red );
  draw rp(ra)--rp(qb);
  draw rp(pos) withpen pencircle scaled 3pt;
  P = (0,ray*sind(theta),ray*(1+cosd(theta)));
  Q = P-len*red;
  draw rp(black)--rp(pos)--rp(P)--rp(Q)--rp(pos);
  draw rp(ra)--rp((-margi-len,dist,0)) dashed withdots withcolor 0.5white;
  draw rp(P)--rp(Q) withcolor blue withpen pencircle scaled 2pt;
endfig;

end.