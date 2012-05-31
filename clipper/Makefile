CC=g++
PREFIX?=/opt
CFLAGS?=-O3
CPPFLAGS?=-O3

PKG=ea-utils
REL := $(shell svnversion 2>/dev/null | perl -ne 'print $$1 if /:(\d+)/' )
VER := $(shell grep '%define ver' ${PKG}.spec | perl -ne 'print $$1 if / (\S+) *$$/')

SRC=fastq-clipper.c fastq-mcf.c fastq-multx.c fastq-join.c fastq-stats.cpp
BIN=fastq-mcf fastq-multx fastq-join  fastq-stats fastq-clipper

all: $(BIN)

install: $(BIN) $(PREFIX)/bin/fastq-clipper $(PREFIX)/bin/fastq-mcf $(PREFIX)/bin/fastq-multx $(PREFIX)/bin/fastq-join $(PREFIX)/bin/fastq-stats

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

$(PKG).tar.gz: Makefile $(SRC) $(PKG).spec fastq-lib.cpp fastq-lib.h sam-stats.cpp fastq-stats.cpp README CHANGES
	tar -cvzf $(PKG).tar.gz $^

$(PKG).${VER}-${REL}.tar.gz: $(PKG).tar.gz
	cp $< $@

%: %.c fastq-lib.cpp
	$(CC) $(CFLAGS) fastq-lib.cpp -o $@ $<

%: %.cpp fastq-lib.cpp
	$(CC) $(CFLAGS) fastq-lib.cpp -o $@ $<

# this requires bamtools and google's sparsehash headers
sam-stats: sam-stats.cpp
	$(CC) $(CFLAGS) -I. fastq-lib.cpp -lbamtools -o $@ $<
