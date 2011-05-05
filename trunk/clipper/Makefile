CC=g++
PREFIX?=/opt
CFLAGS?=-O3

PKG=ea-utils
SRC=fastq-clipper.c fastq-mcf.c fastq-multx.c fastq-join.c
BIN=$(SRC:.c=)

%: %.c
	$(CC) $(CFLAGS) -o $@ $<

all: $(BIN)

install: $(BIN) $(PREFIX)/bin/fastq-clipper $(PREFIX)/bin/fastq-mcf $(PREFIX)/bin/fastq-multx $(PREFIX)/bin/fastq-join

$(PREFIX)/bin:
	mkdir -p $(PREFIX)/bin

$(PREFIX)/bin/%: % $(PREFIX)/bin
	cp $< $@

dist: $(PKG).tar.gz

$(PKG).tar.gz: Makefile $(SRC) $(PKG).spec
	perl -i -pe '$$r=`svnversion | cut -f 2 -d :`; $$r=~ s/M$$//; s/%define rel \w+/%define rel $$r/' ea-utils.spec
	tar -cvzf $(PKG).tar.gz $^
