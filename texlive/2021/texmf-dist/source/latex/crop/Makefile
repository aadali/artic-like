NAME=crop
VERSION=`cat VERSION`
ARCHIVE_NAME=$(NAME).zip
ARCHIVE_CONTENTS=$(NAME).dtx $(NAME).pdf Makefile README $(NAME).ins

all: $(NAME).sty $(NAME).pdf VERSION

archive: $(ARCHIVE_CONTENTS)
	rm -rf $(NAME)/
	mkdir $(NAME)/
	cp $(ARCHIVE_CONTENTS) $(NAME)/
	zip -r $(ARCHIVE_NAME) $(NAME)

$(NAME).dtx: $(NAME).dtx.in
	nancy $(NAME).dtx.in . > $(NAME).dtx

$(NAME).pdf: $(NAME.dtx)
	latexmk $(NAME).dtx

$(NAME).sty: $(NAME).ins
	tex $(NAME).ins

$(NAME).ins:
	pdflatex $(NAME).dtx

release: archive
	git diff --exit-code && \
	git tag -a -m "Release tag" "v$(VERSION)" && \
	git push && git push --tags
