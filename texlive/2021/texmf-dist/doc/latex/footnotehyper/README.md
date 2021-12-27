<!-- -->

    Source:  footnotehyper.dtx (v1.1e 2021/08/13)
    Author:  Jean-Francois Burnol
    Info:    hyperref aware footnote.sty
    License: LPPL 1.3c
    Copyright (C) 2016-2021 Jean-Francois Burnol <jfbu at free dot fr>.

ABSTRACT
========

The `footnote` package by Mark Wooding (`1997/01/28` `1.13`)
allows to gather footnotes (`\begin{savenotes}`) and later insert
them (after `\end{savenotes}`) at the bottom of the page, even
if the intervening material consists of tabulars, minipages or
framed contents for example. One can also use the
`\savenotes/\spewnotes` syntax.

Also, `footnote.sty` provides a `footnote` environment which
allows to insert verbatim material.

Earlier releases of the present `footnotehyper` package added
patches for `hyperref` compatibility and some bugfixes, addressing
in particular the incompatibility with `color/xcolor`, and with
`babel-frenchb`, and also fixing the `footnote` environment with
optional argument `[NUM]`. Since `v0.99` all macros are defined
internally and the `footnote` package is not loaded at all.

The same user interface is kept. Since `v1.0` it is possible to
use `footnotehyper` also in absence of `hyperref` or when the latter is
loaded with its `hyperfootnotes=false` option. The order of loading of
`footnotehyper` and `hyperref` is inconsequential.

INSTALLATION
============

The recommended way is to first extract the package (.sty)
and driver (.tex) files from footnotehyper.dtx via

    tex footnotehyper.dtx

and then produce the documentation via

    latex footnotehyper.tex (twice)
    dvipdfmx footnotehyper.dvi

It is also possible to produce simultaneously the package
and the documentation via one of these two routes:

    pdflatex footnotehyper.dtx (twice)

or

    latex footnotehyper.dtx (twice)
    dvips
    ps2pdf

The method using the extracted file footnotehyper.tex produces
the smallest pdf file and is the officially preferred one as
it allows to set options in footnotehyper.tex to customize the
footnotehyper.pdf file:

- scrdoc class options (paper size, font size, ...)
- with or without source code,
- with dvipdfmx or with latex+dvips or pdflatex.

Installation:

    footnotehyper.sty    -> TDS:tex/latex/footnotehyper/footnotehyper.sty
    footnotehyper.dtx    -> TDS:source/latex/footnotehyper/footnotehyper.dtx
    footnotehyper.pdf    -> TDS:doc/latex/footnotehyper/footnotehyper.pdf
    README.md            -> TDS:doc/latex/footnotehyper/README.md

The other files may be discarded.

LICENSE
=======

This Work may be distributed and/or modified under the conditions
of the LaTeX Project Public License, version 1.3c. This version of
this license is in:

> <http://www.latex-project.org/lppl/lppl-1-3c.txt>

and the latest version of this license is in:

> <http://www.latex-project.org/lppl.txt>

Version 1.3 or later is part of all distributions of
LaTeX version 2005/12/01 or later.

The Author of this Work is:

- Jean-Francois Burnol `<jfbu at free dot fr>`

This Work consists of the main source file footnotehyper.dtx
and the derived files footnotehyper.sty, footnotehyper.tex,
footnotehyper.pdf, footnotehyper.dvi.
