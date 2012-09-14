# set these if you want to override

PREFIX?=/usr
CPPFLAGS?=-O3
REL := $(shell svnversion 2>/dev/null | perl -ne 'print $$1 if /:(\d+)/' )
VER := $(shell perl -e 'do "perl/lib/Text/Tidx.pm"; print "$$Text::Tidx::VERSION\n"' )

all: tidx tidx.a perl/Tidx.bs

debug:
	CPPFLAGS=-g	${MAKE} $(MFLAGS) all

# no fPIC for linking main module
tidx: tidx.cpp tidx-lib.cpp fastq-lib.cpp utils.cpp tidx.h

# fPIC needed for shared linkage of perl module
tidx.a: tidx-lib.cpp fastq-lib.cpp utils.cpp tidx.h
	g++ -O3 -fPIC -c tidx-lib.cpp fastq-lib.cpp utils.cpp
	ar rvs tidx.a tidx-lib.o fastq-lib.o utils.o

perl/Tidx.bs perl/Makefile: perl/Makefile.PL tidx.a
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
	tar -cvzf tidx.$(VER)-$(REL).tar.gz $^ sparsehash
