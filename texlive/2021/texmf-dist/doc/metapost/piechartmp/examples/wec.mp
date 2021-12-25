%H
%H piechartMP example - world energy consumption 2000
%H Jens-Uwe Morawski, 2002
%H


%I The piechartMP documentation includes many useful
%I examples which not require MetaPost experience.
%I In contrast this example employs many MetaPost
%I commands and definitions. Thus it is not meant
%I for piechart "beginners".
%I The chief purpose of this example is to show
%I how you can draw pie-charts with different data-sets
%I using ResetSegments and how to save them in
%I MetaPost pictures.
%I The file wec.pdf in the same directory is the
%I result of this MP-file.

input piechartmp ;

%D Helvetica-Bold
% This should work on most systems. If not, comment the next line.
defaultfont := "phvb8r" ;

%D 
pair NameShift ; NameShift := (-3.5cm,0) ;
Scale := 8mm ;

%D define some colors
color oil_color, gas_color, coal_color, nuclear_color, water_color ;
oil_color     := (0.4, 0  , 0  ) ;
gas_color     := (0  , 0.5, 0  ) ;
coal_color    := (0.5, 0.3, 0  ) ;
nuclear_color := (0.8, 0.8, 0.2) ;
water_color   := (0.3, 0.3, 1.0) ;


ChartCounter := 0 ;
picture Charts[] ;

%D
def AChart (expr Size, Name, Water, Nuclear, Coal, Gas, Oil) =
	begingroup;
	save s , Radius ; string s ;

	% the 'Size' is the absolute energy consumption
	% that means the 'Size' is proportional to the area
	% of the pie-chart and not the radius.
	Radius := (sqrt Size)*Scale ;

	
	s := Name & "; " & (decimal Size) & " billion t MCU" ;
	
	ResetSegments ; % not documented in piechartMP doc

		if Water>0:
			Segment(Water, "", water_color) ;
		fi;
		if Nuclear>0:
			Segment(Nuclear, "", nuclear_color) ;
		fi;
		if Coal>0:
			Segment(Coal, "", coal_color) ;
		fi;
		if Gas>0:
			Segment(Gas, "", gas_color) ;
		fi;
		if Oil>0:
			Segment(Oil, "", oil_color) ;
		fi;

		PieChart(Radius, 0.3, 65, 90, 0) ;
		label.lft (s, NameShift) ;
		ChartCounter := ChartCounter + 1 ;
		Charts[ChartCounter] := currentpicture ;
		currentpicture := nullpicture ;
	endgroup;
enddef;

def ShipoutCharts =
	begingroup;
	save h, t ; numeric h, t ;
	h := 0 ; t := 0 ;
	for i=1 upto ChartCounter:
		t := ypart(ulcorner currentpicture) ;
		h := -1*(ypart(lrcorner Charts[i])) ;
		h := h + t + 5mm ;
		draw Charts[i] shifted (0, (h));
	endfor;
	endgroup;
enddef;

beginfig(1);
AChart (3.8, "North-America", 2, 6, 25, 32, 35 ) ;
AChart (3.4, "Asia/Oceania", 2, 5, 46, 10, 38 ) ;
AChart (2.6, "Europe", 4, 16, 18, 20, 42 ) ;
AChart (1.3, "CIS", 2, 6, 18, 50, 24 ) ;
AChart (0.5, "Middle East", 0, 0, 0, 35, 65 ) ;
AChart (0.5, "South-America", 15, 0, 5, 15, 65 ) ;
AChart (0.4, "Africa", 0, 0, 37, 18, 45 ) ;
ShipoutCharts ;

x := 5mm+(xpart(lrcorner currentpicture)) ;
y := ypart(ulcorner currentpicture) ;

fill fullcircle scaled 5mm shifted (x,y) withcolor oil_color;
label.rt ("Oil", (x+5mm,y)) ;
fill fullcircle scaled 5mm shifted (x,(y-1cm)) withcolor gas_color;
label.rt ("Gas", (x+5mm,y-1cm)) ;
fill fullcircle scaled 5mm shifted (x,(y-2cm)) withcolor coal_color;
label.rt ("Coal", (x+5mm,y-2cm)) ;
fill fullcircle scaled 5mm shifted (x,(y-3cm)) withcolor nuclear_color;
label.rt ("Nuclear Power", (x+5mm,y-3cm)) ;
fill fullcircle scaled 5mm shifted (x,(y-4cm)) withcolor water_color;
label.rt ("Water Power", (x+5mm,y-4cm)) ;

%D 
defaultfont:="ptmbi8r" ; defaultscale:= 2;
label.top("World Energy Consumption 2000", (0, y+5mm)) ;

endfig;
end