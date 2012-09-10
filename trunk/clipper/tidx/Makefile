# set these if you want to override

PREFIX?=/usr
CPPFLAGS?=-O3

all: tidx tidx.a perl/Tidx.bs

debug:
	CPPFLAGS=-g	${MAKE} $(MFLAGS) all

tidx: tidx.cpp tidx-lib.cpp utils.cpp fastq-lib.cpp tidx.h

tidx.a: tidx-lib.cpp fastq-lib.cpp utils.cpp tidx.h
	g++ -O3 -fPIC tidx-lib.cpp fastq-lib.cpp utils.cpp -shared -o tidx.a

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
	cd perl; perl test.pl > test.tmp; diff test.tmp test.ok; rm test.tmp
