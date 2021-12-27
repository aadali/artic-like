# In order to compile the sources and to get a pdf file,
# use `make' or `./readme.sh' or:

# 1. run metapost THREE TIMES:
mpost drv-guide.mp
mpost drv-guide.mp
mpost drv-guide.mp

# 2. run latex (three times):
latex drv-guide.tex
latex drv-guide.tex
latex drv-guide.tex

# 3. run dvips (once):
dvips drv-guide.dvi -o

# 4. run ps2pdf (once):
ps2pdf drv-guide.ps
