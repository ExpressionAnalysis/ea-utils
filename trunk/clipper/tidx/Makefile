# set these if you want to override

PREFIX?=/usr
CPPFLAGS?=-O3
PERL_INC?=$(shell perl -e 'use Config; print $$Config{archlib};' )

all: tidx tidx_wrap.so perl/Tidx.bs

debug:
	CPPFLAGS=-g	${MAKE} $(MFLAGS) all

tidx: tidx.cpp tidx-lib.cpp utils.cpp fastq-lib.cpp tidx.h

tidx_wrap.cxx: tidx.i
	mkdir -p perl/blib/lib/Text/Tidx 
	swig -c++ -perl -outdir perl/blib/lib/Text/Tidx tidx.i

tidx_wrap.so: tidx_wrap.cxx tidx-lib.cpp fastq-lib.cpp utils.cpp tidx.h
	g++ -O3 -fPIC -I${PERL_INC}/CORE  tidx_wrap.cxx tidx-lib.cpp fastq-lib.cpp utils.cpp -shared -o tidx_wrap.so

perl/Tidx.bs perl/Makefile: perl/Makefile.PL tidx_wrap.so
	cd perl; perl Makefile.PL; ${MAKE} $(MFLAGS) INSTALL_BASE=${PREFIX}

${PREFIX}/bin/tidx: tidx
	cp tidx ${PREFIX}/bin

${PREFIX}/lib/tidx_wrap.so: tidx_wrap.so
	cp tidx_wrap.so ${PREFIX}/lib

install: ${PREFIX}/bin/tidx ${PREFIX}/lib/tidx_wrap.so perl/Makefile
	cd perl; ${MAKE} $(MFLAGS) install

clean:
	-$(RM) -f tidx *.so perl/*.bs tidx_wrap.cxx Makefile.old
	cd perl; ${MAKE} $(MFLAGS) clean

test: all
	cd perl; perl test.pl > test.tmp; diff test.tmp test.ok; rm test.tmp
