#ifdef __cplusplus
extern "C" {
#endif

#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "ppport.h"

#ifdef __cplusplus
}
#endif

#include "tidx.h"

MODULE = Text::Tidx		PACKAGE = Text::Tidx		

tidx *
tidx::new(char * path)
 
void
tidx::DESTROY()

const char * 
tidx::lookup_c(const char *chr, int pos, const char *msep);

const char * 
tidx::lookup_cr(const char *chr, int beg, int end, const char *msep);
