CC=g++
PREFIX?=/opt
CFLAGS?=-O3 -I. 
CPPFLAGS?=-O3 -I.

PKG=ea-utils
REL := $(shell svnversion 2>/dev/null | perl -ne 'print $$1 if /:(\d+)/' )
VER := $(shell grep '%define ver' ${PKG}.spec | perl -ne 'print $$1 if / (\S+) *$$/')

SRC=fastq-clipper.c fastq-mcf.c fastq-multx.c fastq-join.c fastq-stats.cpp
BIN=fastq-mcf fastq-multx fastq-join fastq-stats fastq-clipper sam-stats varcall

all: $(BIN)

debug: 
	CPPFLAGS=-g ${MAKE} $(MFLAGS) varcall

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

.PHONY: getrel debug

$(PKG).spec: $(PKG).spex
	perl -pe 's/%RELEASE%/${REL}/' $(PKG).spex > $(PKG).spec

$(PKG).tar.gz: Makefile $(SRC) $(PKG).spec fastq-lib.cpp fastq-lib.h sam-stats.cpp fastq-stats.cpp varcall.cpp utils.h README CHANGES google sparsehash samtools/*.c
	rm -rf $(PKG).${VER}-${REL}
	mkdir $(PKG).${VER}-${REL}
	mkdir $(PKG).${VER}-${REL}/tidx
	mkdir $(PKG).${VER}-${REL}/samtools
	cp -nr $^ $(PKG).${VER}-${REL}
	cp -nr tidx/*.cpp tidx/*.h $(PKG).${VER}-${REL}/tidx
	cp -nr samtools/*.c samtools/*.h samtools/Makefile $(PKG).${VER}-${REL}/samtools
	tar --exclude=".svn" -cvzf $(PKG).tar.gz $(PKG).${VER}-${REL}
	rm -rf $(PKG).${VER}-${REL}

disttest: $(PKG).tar.gz
	tar -xzvf $(PKG).tar.gz
	cd $(PKG).${VER}-${REL} && make
	rm -rf $(PKG).${VER}-${REL}

$(PKG).${VER}-${REL}.tar.gz: $(PKG).tar.gz
	cp $< $@

%: %.c fastq-lib.cpp
	$(CC) $(CFLAGS) fastq-lib.cpp -o $@ $<

%: %.cpp fastq-lib.cpp
	$(CC) $(CFLAGS) fastq-lib.cpp -o $@ $<

# why the libbam.a doesn't work?  not sure... *.o works
sam-stats: sam-stats.cpp samtools/libbam.a samtools/bam.h
	$(CC) $(CFLAGS) samtools/*.o -lz fastq-lib.cpp $< -o $@

samtools/libbam.a: samtools/*.c
	cd samtools && make libbam.a

varcall: varcall.cpp fastq-lib.cpp tidx/tidx-lib.cpp
	$(CC) $(CFLAGS) fastq-lib.cpp tidx/tidx-lib.cpp -o $@ $<
