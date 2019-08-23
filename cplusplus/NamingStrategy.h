/*
* These macros define the naming strategy needed for a fortran77
* routine to call a C routine, and whether to build so they may be
* called from C or f77.  For the f77 call C interface, ADD_ assumes that
* f77 calls expect C routines to have an underscore postfixed to the name
* (Suns, and the Intel expect this).  NOCHANGE indicates that f77 expects
* the name called by fortran to be identical to that compiled by C
* (RS6K's do this).  UPCASE says it expects C routines called by fortran
* to be in all upcase (CRAY wants this).  The variable F77_CALL_C is always
* set to one of these values.  
*/
#define ADD_     0
#define NOCHANGE 1
#define UPCASE   2
#define F77ISF2C 3
#define C_CALL   4

#ifdef _WIN32 /* _WIN32 is usually defined by compilers targeting 32 or 64 bit Windows systems */
#define F77_CALL_C UPCASE
#else
#define F77_CALL_C NOCHANGE
#endif

#if (F77_CALL_C == UPCASE)
/*
* These defines set up the naming scheme required to have a fortran 77
* routine call a C routine for the following Fortran to C interface:
*           FORTRAN CALL               C DECLARATION
*           call dgebs2d(...)          void DGEBS2D(...)
*/
/*
* Support routines
*/
#define save_series_                   SAVE_SERIES
#define sparseassembler_kill_                   SPARSEASSEMBLER_KILL
#define sparseassembler_init_                   SPARSEASSEMBLER_INIT
#define sparseassembler_addelemmatrix_                   SPARSEASSEMBLER_ADDELEMMATRIX
#define sparseassembler_addelemmatrixS_                   SPARSEASSEMBLER_ADDELEMMATRIXS
#define sparseassembler_addnonsymmatrix_                   SPARSEASSEMBLER_ADDNONSYMMATRIX
#define sparseassembler_getsparse_                   SPARSEASSEMBLER_GETSPARSE
#define sparseassembler_getnz_                   SPARSEASSEMBLER_GETNZ



#endif