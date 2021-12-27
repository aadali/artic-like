The `unicode-bidi` package
========

The experimental `unicode-bidi` package allows you to mix non-RTL script with RTL script without any markup.

Usage
-----

To use the `unicode-bidi` package, you just need to load it

````latex
\usepackage{unicode-bidi}
````

Example
-------

Here is a sample `xepersian` document:

````latex
\documentclass{article}
\usepackage{xepersian}
\usepackage{unicode-bidi}
\settextfont{Yas}
\begin{document}
این یک آزمایش است (یک آزمایش)
One, two, and three (numbers) or [numbers]
و آزمایش بعدی,
This is a [real test].
و این هم یک [آزمایش دیگر] است.

نام نخستین شرکت هواپیمایی ایران
Iran Air or \hbox{Homa airlines}
بوده است که به بسیاری از نقاط دنیا پرواز داشته است.


من نام
\hbox{من و you and we ما می‌شویم.}
را دوست می‌دارم.

Joanne 123 and you
نام من در زبان انگلیسی می‌باشد.

\begin{latin}
This is a test as you can see

\end{latin}

نام من به زبان انگلیسی 
Joanne Patterson
هست.
\end{document}
````


    
    
Reporting bugs
-------
Please use [the issue tracker][1] to report bugs.


  [1]: https://github.com/vafa/unicode-bidi/issues