use 5.010001;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
$CC = "g++";
$EXE_EXT = $^O =~ /win32/i ? ".exe" : "";

WriteMakefile(
    NAME              => 'Text::Tidx',
    VERSION_FROM      => 'lib/Text/Tidx.pm', # finds $VERSION
    H                 => ['sparsehash'],
    OPTIMIZE          => '-O3',
    PREREQ_PM         => {}, # e.g., Module::Name => 1.1
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM  => 'lib/Text/Tidx.pm', # retrieve abstract from module
       AUTHOR         => 'A. U. Thor <earonesty@>') : ()),
#    LIBS              => ['-ltidx'], # e.g., '-lm'
    CC                => $CC,
    LD                => '$(CC)',
    XSOPT             => '-C++',
    TYPEMAPS          => ['perlobject.map'],
#    'dist'            => {  },
    DEFINE            => '', # e.g., '-DHAVE_SOMETHING'
    INC               => '-I.', # e.g., '-I. -I/usr/include/other'
	# Un-comment this if you add C files to link with later:
    OBJECT            => '$(O_FILES)', # link all the C files too
    EXE_FILES         => ['tidx'.$EXE_EXT],
);


sub MY::postamble 
{
# todo... this probably doesn't work on windows/mig
return <<CMD
tidx: fastq-lib.o tidx-lib.o tidx.o utils.o
	\$(CC) -O3 fastq-lib.o tidx-lib.o tidx.o utils.o -o tidx$EXE_EXT
CMD
}

