% joinedemptylines.mp
% L. Nobre G.
% 2001

input featpost3Dplus2D;
 
numeric altura;
altura = 1;

f := (2,3,1.8+altura);
Spread := 125;
LightSource := 0.7*(4,-1,4);
ShadowOn := true;

def helixone( expr TheVal ) = 
	begingroup
		numeric par[], cred, cblu, param;
		par1 = 3*360;
		par2 = 0.3;
		par3 = 0.8;
		param = par1*TheVal;
		cred = par2*cosd( param );
		cblu = par2*sind( param );
		( (cred,cblu,par3*(-1+TheVal)+altura) )
	endgroup 
enddef;

def helixtwo( expr TheVal ) = 
	begingroup
		numeric par[], cred, cblu, param;
		par1 = 3*360;
		par2 = 0.1;
		par3 = 0.8/3;
		param = par1*TheVal;
		cred = par2*cosd( param );
		cblu = par2*sind( param );
		( (cred+0.2,cblu,par3*TheVal+altura) )
	endgroup 
enddef;

beginfig(1); 
  emptyline( false, 0.91, black, 0.5background, 360, 0.88, 2, helixone );
  emptyline( true, 0.91, black, 0.5background, 120, 0.88, 2, helixtwo );
  drawarrow rp((0.5,0,altura))--rp((0.32,0,altura));
endfig;

end;
