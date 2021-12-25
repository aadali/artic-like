% emptylines.mp
% L. Nobre G.
% 2001

input featpost3Dplus2D;
 

def closec( expr TheVal ) = 
	begingroup
		numeric par[], cred, cblu, param;
		par1 = 360;
		par2 = 0.3;
		param = 360*TheVal;
		cred = par2*cosd( param );
		cblu = par2*sind( param );
		( (cred,cblu,1) )
	endgroup 
enddef;

f := (4,3,2);
Spread := 170;
ShadowOn := true;

beginfig(1); 
	     emptyline( true, 1, black, 0.5background, 360, 0.88, 2, closec );
endfig;

end;
