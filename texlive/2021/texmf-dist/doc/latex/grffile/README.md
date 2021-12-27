# grffile legacy LaTeX package

By default this package just loads the standard graphicx package
Which can include files with spaces or multiple dots.

The original code, by Heiko Oberdiek, is still available as an option.
This extended the original LaTeX graphics code to support spaces and
multiple dots in filenames.

If required use

\usepackage{grffile}[=v1]

to use the original version of this package (grffile-2017-06-30.sty).

tex grffile.dtx

will extract grffile.sty and grffile-2017-06-30.sty