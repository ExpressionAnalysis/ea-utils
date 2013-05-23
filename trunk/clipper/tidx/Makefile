# set these if you want to override

PREFIX?=/usr
CPPFLAGS?=-O3 -I.
REL := $(shell svnversion 2>/dev/null | perl -ne 'print $$1 if /:(\d+)/' )
VER := $(shell perl -e 'do "perl/lib/Text/Tidx.pm"; print "$$Text::Tidx::VERSION\n"' )

all: tidx perl/Tidx.bs

debug:
	CPPFLAGS="-I. -g" ${MAKE} $(MFLAGS) all

# no fPIC for linking main module
tidx: tidx.cpp tidx-lib.cpp fastq-lib.cpp utils.cpp tidx.h

perl/Tidx.bs perl/Makefile: perl/Makefile.PL
	cd perl; perl Makefile.PL; ${MAKE} $(MFLAGS) INSTALL_BASE=${PREFIX}

${PREFIX}/bin/tidx: tidx
	cp tidx ${PREFIX}/bin

install: ${PREFIX}/bin/tidx
	cd perl; ${MAKE} $(MFLAGS) install

clean:
	-$(RM) -f tidx *.a perl/*.bs Makefile.old
	cd perl; ${MAKE} $(MFLAGS) clean

test: all
	cd perl; make test

dist: tidx.cpp tidx-lib.cpp fastq-lib.cpp utils.cpp tidx.h
	tar --exclude=".svn" -cvhzf tidx.$(VER)-$(REL).tar.gz $^ sparsehash
