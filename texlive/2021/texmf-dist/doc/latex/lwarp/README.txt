
LaTeX lwarp package v0.902   README.txt

Files included are:

lwarp.dtx: The documented source code.
lwarp.ins: The documentation driver.
lwarp.pdf: The documentation.
lwarp_tutorial.txt: More documentation.  A sample tutorial.
lwarpmk.lua: A utility program which compiles print or HTML versions.
lwarp_mathjax.txt: A script inserted when MathJax is used.
lwarp_baseline_marker.png: A tiny transparent image used to
    improve SVG math width and baselines, to be located in
    tex/latex/lwarp and thus be found by kpsewhich.

Derived by compiling lwarp.ins:
lwarp.sty: The lwarp package.
lwarp-*.sty: Numerous compatibility files for other packages.


NOTE FOR TEX DISTRIBUTION MAINTAINERS:
For a TeX distribution, an executable called "lwarpmk" should
be created in each bin directory, which is a link to or which
calls the "lwarpmk.lua" script found in the scripts directory.
See the installation instructions in the documentation.


The documentation includes a file listing "tutorial.tex" which may be
copied/pasted from the documentation into an editor, or it may
be copied from the file lwarp_tutorial.tex, which is found in the
documentation directory.


License:
This work may be distributed and/or modified under the
conditions of the LaTeX Project Public License, either version 1.3
of this license or (at your option) any later version.
The latest version of this license is in
  http://www.latex-project.org/lppl.txt
and version 1.3 or later is part of all distributions of LaTeX
version 2005/12/01 or later.


Copyright 2016-2021 Brian Dunn

Homepage: http://BDTechConcepts.com
Email: bd@BDTechConcepts.com
