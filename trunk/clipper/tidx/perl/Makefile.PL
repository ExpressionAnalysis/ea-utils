# File : Makefile.PL
use ExtUtils::MakeMaker;
WriteMakefile(
    'NAME'          => 'Text::Tidx::SWIG',            # Name of package buildiung
    'DISTNAME'      => 'Text-Tidx',                   # Name of package installing
    'OBJECT'        => 'tidx_wrap.so',             # Object files while building
# wrapper, and built-swig thing
# Note: there's probably a way to get swig to do this without another layer....
    'PM'            => {
            'Text/Tidx.pm'=>'$(INST_LIBDIR)/Tidx.pm', 
            'blib/lib/Text/Tidx/SWIG.pm'=>'$(INST_LIBDIR)/Tidx/SWIG.pm',
        },
);
