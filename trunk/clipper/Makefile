CC=g++
PREFIX?=/opt
CFLAGS?=-O3

PKG=ea-utils
REL := $(shell svnversion 2>/dev/null | perl -ne 'print $$1 if /:(\d+)/' )
VER := $(shell grep '%define ver' ${PKG}.spec | perl -ne 'print $$1 if / (\S+) *$$/')

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
	perl -i -pe 's/%define rel \w*/%define rel ${REL}/' ea-utils.spec
	tar -cvzf $(PKG).tar.gz $^
	[ -n "${REL}" ] && cp $(PKG).tar.gz $(PKG).${VER}-${REL}.tar.gz
