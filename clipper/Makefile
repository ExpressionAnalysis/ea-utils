CC=g++
PREFIX?=/opt
CFLAGS?=-O3

PKG=ea-utils
REL := $(shell svnversion 2>/dev/null | perl -ne 'print $$1 if /:(\d+)/' )
ALWAYS := $(shell grep ${REL} $(PKG).spec || touch $(PKG).spex && echo 'OK' )
VER := $(shell grep '%define ver' ${PKG}.spex | perl -ne 'print $$1 if / (\S+) *$$/')

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

dist: $(PKG).${VER}-${REL}.tar.gz

$(PKG).spec: $(PKG).spex
	perl -pe 's/%RELEASE%/${REL}/' $(PKG).spex > $(PKG).spec

$(PKG).tar.gz: Makefile $(SRC) $(PKG).spec README CHANGES
	tar -cvzf $(PKG).tar.gz $^

$(PKG).${VER}-${REL}.tar.gz: $(PKG).tar.gz
	cp $< $@
