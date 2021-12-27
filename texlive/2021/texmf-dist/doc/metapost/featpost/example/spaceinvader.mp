% spaceinvader.mp
% L. Nobre G.
% 2012

prologues := 1;

vardef renewcolore =
  color tmpcolo;
  tmpcolo = (uniformdeviate(0.8),uniformdeviate(0.8),uniformdeviate(0.8));
  drawoptions( withcolor tmpcolo );
enddef;

primarydef A xor B = ((not A) and B) or (A and not B) enddef;

def strai( expr FlagVer, FlagForward, IndI, IndJ ) =
  begingroup
    path tmppath;
    if FlagVer:
      tmppath = (IndI*u,(IndJ-1)*u+mar)---(IndI*u,IndJ*u-mar);
    else:
      tmppath = ((IndI-1)*u+mar,IndJ*u)---(IndI*u-mar,IndJ*u);
    fi;
    if not FlagForward:
      tmppath := reverse tmppath;
    fi;
    ( tmppath )
  endgroup
enddef;
      
beginfig(1);
  string st[];
  st16 = "..............M.";
  st15 = ".............M..";
  st14 = "MM..........M...";
  st13 = "..MMM......M....";
  st12 = ".....M....M.....";
  st11 = "MMM...M..M...MMM";
  st10 = "..MMM.MMMM.MMM..";
  st9  = "....MMMMMMMM....";
  st8  = "....M..MM..M....";
  st7  = "....MMMMMMMM....";
  st6  = "....MMM..MMM....";
  st5  = "...MMMMMM.MMMMM.";
  st4  = "..MM..MMMM..M.MM";
  st3  = "MMM..MM..MMM...M";
  st2  = "M...MM.....MMM.M";
  st1  = "...MM........M..";
  numeric i, j, hv[][][];
  boolean b[][];
  numeric refi, refj, pn, u, mar, numa, numb, numc;
  boolean startfound, forward, statever, pathfinished;
  path frontier[];
  picture borderpic;
  color cola, colb, colc;
  pen thepen;
  u = 6mm;
  thepen = pencircle scaled 0.2u;
  mar = 0.3u;
  startfound = false;
  pn = 0;
  for j=1 upto 16:
    for i=1 upto 16:
      if substring (i-1,i) of st[j] = ".":
	b[i][j] = false;
      else:
	b[i][j] = true;
      fi;
    endfor;
  endfor;
  for i=0 upto 17:
    b[i][0] := false;
    b[i][17] := false;
    b[0][i] := false;
    b[17][i] := false;
  endfor;
  for j=0 upto 16:
    for i=0 upto 16:
      if b[i][j] xor b[i+1][j]:
	hv[1][i][j] = 1;
	draw strai( true, true, i, j ); 
      else:
	hv[1][i][j] = 0;
      fi;
      if b[i][j] xor b[i][j+1]:
	hv[0][i][j] = 1;
	draw strai( false, true, i, j );
      else:
	hv[0][i][j] = 0;
      fi; 
    endfor;
  endfor;
  for i=0 upto 17:
    hv[0][17][i] := 0;
    hv[1][i][17] := 0;
  endfor;
  borderpic = currentpicture;
endfig;

def trioone( expr flstat, flforw ) =
  begingroup
    color tmpcolor;
    if flstat and flforw:
      tmpcolor = (0,0,0);
    elseif flstat and (not flforw):
      tmpcolor = (0,1,-1);
    elseif (not flstat) and flforw:
      tmpcolor = (1,0,1);
    else:
      tmpcolor = (1,-1,0);
    fi;
    ( tmpcolor )
  endgroup
enddef;

def triotwo( expr flstat, flforw ) =
  begingroup
    color tmpcolor;
    if flstat and flforw:
      tmpcolor = (1,0,1);
    elseif flstat and (not flforw):
      tmpcolor = (1,0,-1);
    elseif (not flstat) and flforw:
      tmpcolor = (0,1,0);
    else:
      tmpcolor = (0,-1,0);
    fi;
    ( tmpcolor )
  endgroup
enddef;

def triothr( expr flstat, flforw ) =
  begingroup
    color tmpcolor;
    if flstat and flforw:
      tmpcolor = (0,1,0);
    elseif flstat and (not flforw):
      tmpcolor = (0,0,-1);
    elseif (not flstat) and flforw:
      tmpcolor = (1,0,0);
    else:
      tmpcolor = (1,-1,1);
    fi;
    ( tmpcolor )
  endgroup
enddef;

beginfig(2);
  % draw borderpic;
  forever:
    startfound := false;
    for j=0 upto 16:
      for i=1 upto 16:
	if (not startfound) and (hv[0][i][j] = 1):
	  startfound := true;
	  refi := i;
	  refj := j;
	fi;
      endfor;
    endfor;
    exitunless startfound;
    statever := false;
    forward := true;
    frontier[incr(pn)] = strai( statever, forward, refi, refj );
    hv[0][refi][refj] := 2;
    pathfinished := false;
    forever:
      cola := trioone(statever,forward);
      colb := triotwo(statever,forward);
      colc := triothr(statever,forward);
      numa := hv[redpart cola][refi+greenpart cola][refj+bluepart cola];
      numb := hv[redpart colb][refi+greenpart colb][refj+bluepart colb];
      numc := hv[redpart colc][refi+greenpart colc][refj+bluepart colc];
      if numa = 1:
	refi := refi+greenpart cola;
	refj := refj+bluepart cola;
	forward := statever xor forward;
	statever := not statever;
	frontier[pn] := frontier[pn]..strai( statever, forward, refi, refj );
	hv[redpart cola][refi][refj] := 2;
      elseif numa = 2:
	frontier[pn] := frontier[pn]..cycle;
	pathfinished := true;
      elseif numb = 1:
	refi := refi+greenpart colb;
	refj := refj+bluepart colb;
	frontier[pn] := frontier[pn]..strai(statever,forward,refi,refj);
	hv[redpart colb][refi][refj] := 2;
      elseif numb = 2:
	frontier[pn] := frontier[pn]..cycle;
	pathfinished := true;
      elseif numc = 1:
	refi := refi+greenpart colc;
	refj := refj+bluepart colc;
	forward := not (statever xor forward);
	statever := not statever;
	frontier[pn] := frontier[pn]..strai( statever, forward, refi, refj );
	hv[redpart colc][refi][refj] := 2;
      elseif numc = 2:
	frontier[pn] := frontier[pn]..cycle;
	pathfinished := true;	
      fi;
      exitif pathfinished;
    endfor;
    renewcolore;
    fill frontier[pn] withpen thepen;
  endfor;
  show pn;
endfig;

end.