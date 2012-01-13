CC=g++
PREFIX?=/opt
CFLAGS?=-O3

PKG=ea-utils
REL := $(shell svnversion 2>/dev/null | perl -ne 'print $$1 if /:(\d+)/' )
VER := $(shell grep '%define ver' ${PKG}.spex | perl -ne 'print $$1 if / (\S+) *$$/')

SRC=fastq-clipper.c fastq-mcf.c fastq-multx.c fastq-join.c
BIN=$(SRC:.c=)

all: $(BIN)

install: $(BIN) $(PREFIX)/bin/fastq-clipper $(PREFIX)/bin/fastq-mcf $(PREFIX)/bin/fastq-multx $(PREFIX)/bin/fastq-join

$(PREFIX)/bin:
	mkdir -p $(PREFIX)/bin

$(PREFIX)/bin/%: % $(PREFIX)/bin
	cp $< $@

dist: getrel $(PKG).${VER}-${REL}.tar.gz

# these shenanigans are done to ensure than the release in the spec file is the same as the subversion release
# a less verbose way should be possible

getrel:
	grep "${REL}" $(PKG).spec || touch $(PKG).spex

.PHONY: getrel

$(PKG).spec: $(PKG).spex
	perl -pe 's/%RELEASE%/${REL}/' $(PKG).spex > $(PKG).spec

$(PKG).tar.gz: Makefile $(SRC) $(PKG).spec sam-stats.cpp fastq-stats README CHANGES
	tar -cvzf $(PKG).tar.gz $^

$(PKG).${VER}-${REL}.tar.gz: $(PKG).tar.gz
	cp $< $@

%: %.c
	$(CC) $(CFLAGS) fastq-lib.cpp -o $@ $<

sam-stats: sam-stats.cpp
	$(CC) $(CFLAGS) -lbamtools -o $@ $<
