
%% This is the documentation of the minim-hatching pattern library. To typeset 
%% this file, say:
%%
%%     luatex --fmt=minim-mp minim-hatching-doc.mp

input minim-hatching

verbatimtex
% work around a bug: pdf.setnames will do nothing if this is not called at least once.
\pdfextension names{} % (TODO: remove after bug is fixed).
\font \titlf = {Latin Modern Roman Dunhill} at 12pt
\font \mainf = {Latin Modern Roman Unslanted} at 10pt \mainf
\parindent 0em \parskip 1em
\rightskip 0em plus 2em
\startmetadata
    author   {Esger Renkema}
    title    {minim-hatching}
    date     {2021-06-01}
    version  {2021/1.0}
    keywords {MetaPost; tiling patterns}
stopmetadata
etex;

def makepage = path page; page = unitsquare xscaled 210mm yscaled 297mm;
    setbounds currentpicture to page shifted (-7.1cm, 4cm-297mm); enddef;
def xyscaled expr p = xscaled xpart p yscaled ypart p enddef;

pickup pencircle scaled 0.4pt; save p; picture p;
save w, h, o; (w, h, o) = (2cm, 1cm, 3mm);

def table(expr n)(text i)(text t) =
    for d = t : x := 0;
        for s = i :
            if d = "=" : draw label(s) shifted (x+o, y+o);
            else: filldraw unitsquare xscaled w yscaled h shifted (x,y)
                withpattern(scantokens(s&" "&d&n)); fi
            x := x + w;
        endfor
        patname(d, n)
        y := y - h;
    endfor enddef;
def double(expr l, r) = y:=y-h;
    filldraw unitsquare xscaled w yscaled h shifted (0,y)
        withpattern(scantokens(l));
    filldraw unitsquare xscaled w yscaled h shifted (w,y)
        withpattern(scantokens(r));
    draw label(l) shifted(2w+o, y+o);
    draw label(r) shifted(3w+o, y+o); enddef;

vardef label(expr t)=
    save p; picture p; p := scantokens ("btex "&t&"\strut etex;");
    setbounds p to unitsquare xyscaled urcorner p; p enddef;
def patname(expr n, s) =
    if n <> "=" : draw label(n&s) shifted(x+o, y+o); fi enddef;
def titl(expr t) =
    p := label("\titlf "&t);
    draw p shifted ((-2o, y+o) - lrcorner p); enddef;

beginfig(1) y := 0;
    titl("Introduction");
    p := btex \vtop{\hsize=9cm\noindent
        This is a small library of monochrome patterns for the minim-mp metapost processor.
        Load the patterns by saying \par
        \strut \qquad input minim-hatching;\par
        at the top of your metapost file or instance.\par
        I recommend printing this file, as patterns do not always display well on a computer screen.
        } etex shifted (0,o);
    draw p; y := y - ypart (urcorner p - lrcorner p) - h;
    titl("Line patterns");
    table(" lines")("", "dense")
        ("=", "horizontal", "vertical", "upward", "downward", "broken");
    titl("Zigzag patterns");
    table("zigzag")("horizontal", "vertical")
        ("=", "" );
    titl("Grid patterns");
    table(" grid")("", "dense")
        ("=", "", "upward", "downward", "diagonal", "horizontal", "vertical");
makepage endfig;

verbatimtex \embedfile
    global file {EUPL-1.2-EN.txt}
    desc {Licence for this package (EUPL).}
    name LICENCE mimetype text/plain
    relation Unspecified uncompressed etex;

beginfig(1) y := 0;
    titl("Dot patterns");
    table(" dots")("sparse", "", "dense")
        ("=", "", "square", "diamond", "heavy", "square heavy", "diamond heavy");
    titl("Cloth patterns");
    double("knit", "purl")
    double("weave", "dense weave")
    double("stitches", "crosses")
    double("chords", "tea cloth")
    y:=y-h; titl("Circle patterns");
    table("")("", "small")
        ("=", "circles", "square circles", "hexagons", "scales");
    titl("Stone patterns");
    double("bricks", "tiles");
    y:=y-2h; titl("Licence");
    p := btex
        \vtop{\hsize=9cm\noindent
        This package may be distributed under the terms of the European Union Public 
        Licence (EUPL) version 1.2 or later. An english version of this licence has 
        been included as an attachment to this file; copies in other languages can be 
        obtained at
        $$\hbox to \hsize{\hss\hyperlink
            url {https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12}\relax
            https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12\endlink\hss}$$
        } etex shifted (0,y+o);
    draw p;
makepage endfig;

bye;


