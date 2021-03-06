Awesome Boxes is all about drawing admonition blocks around text to
inform or alert your readers about something particular. The specific
aim of this package is to use [FontAwesome 5](https://fontawesome.com)
icons to ease the illustration of these blocks.

The idea of admonition blocks comes from the ones you can easily do with
[AsciiDoc](http://asciidoctor.org/docs/user-manual/#admonition).

For more details, please refer to the
[awesomebox.pdf](https://github.com/milouse/latex-awesomebox/blob/master/awesomebox.pdf)
document.

Requirements
============

The following LaTeX packages are required (they should be already
included in your LaTeX distribution):

-   `array`
-   `fontawesome5`
-   `ifthen`
-   `xcolor`
-   `xparse`

Installation
============

Download the `awesomebox.sty` file and put it in the same folder of the
document your are composing.

For system wide installation, please refer to the documentation of your
LaTeX distribution.

Compatibility
=============

This repository also hosts the package `awesomebox-compat`, which
depends on the `fontawesome` package, instead of `fontawesome5`. Apart
from that, it has the exact same features set.

This can be usefull for you if you want to use it with an old LaTeX
distribution, which does not embed `fontawesome5` yet (like Overleaf
system). To use it, download the file `awesomebox-compat.sty` in this
repository, put it near your tex files and just replace your
`\usepackage{awesomebox}` instruction by
`\usepackage{awesomebox-compat}`.

Be aware that icon names changes between FontAwesome and FontAwesome5
and thus using the compatibility package may break your current files.
FontAwesome also requires you to use XeLaTeX or LuaTeX: you cannot use
PDFLaTeX with the compatibility package.

License
=======

Awesome Box is released under the
[WTFPL](http://www.wtfpl.net/txt/copying/). A copy of this license is
distributed in this package.
