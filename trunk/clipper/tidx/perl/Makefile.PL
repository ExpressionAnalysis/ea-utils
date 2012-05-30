# File : Makefile.PL
use ExtUtils::MakeMaker;
WriteMakefile(
    'NAME'    => 'Text::Tidx::SWIG',            # Name of package
    'OBJECT'  => '../tidx_wrap.so',              # Object files
    'PM'      => {
            'Text/Tidx.pm'=>'$(INST_LIBDIR)/Tidx.pm', 
            'blib/lib/Text/Tidx/SWIG.pm'=>'$(INST_LIBDIR)/Tidx/SWIG.pm',
        },
);
