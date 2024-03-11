#ifndef CT_CONFIG_H
#define CT_CONFIG_H

//---------------------------- Version Flags ------------------//
// Cantera version -> this will be a double-quoted string value
#define CANTERA_VERSION "3.0.0"

// Just the major + minor version (that is, 2.2 instead of 2.2.0)
#define CANTERA_SHORT_VERSION "3.0"

//------------------------ Fortran settings -------------------//

// define types doublereal, integer, and ftnlen to match the
// corresponding Fortran data types on your system. The defaults
// are OK for most systems

typedef double doublereal;   // Fortran double precision
typedef int integer;      // Fortran integer
typedef int ftnlen;       // Fortran hidden string length type

// Fortran compilers pass character strings in argument lists by
// adding a hidden argument with the length of the string. Some
// compilers add the hidden length argument immediately after the
// CHARACTER variable being passed, while others put all of the hidden
// length arguments at the end of the argument list. Define this if
// the lengths are at the end of the argument list. This is usually the
// case for most unix Fortran compilers, but is (by default) false for
// Visual Fortran under Windows.
#define STRING_LEN_AT_END

// Define this if Fortran adds a trailing underscore to names in object files.
// For linux and most unix systems, this is the case.
#define FTN_TRAILING_UNDERSCORE 1


#define CT_SUNDIALS_VERSION 53

//-------- LAPACK / BLAS ---------

#define LAPACK_FTN_STRING_LEN_AT_END 1
#define LAPACK_FTN_TRAILING_UNDERSCORE 1
/* #undef CT_USE_LAPACK */

/* #undef CT_USE_SYSTEM_EIGEN */
/* #undef CT_USE_SYSTEM_EIGEN_PREFIXED */
/* #undef CT_USE_SYSTEM_FMT */
/* #undef CT_USE_SYSTEM_YAMLCPP */

//-------------- Optional Cantera Capabilities ----------------------

// Enable Sundials to use an external BLAS/LAPACK library if it was
// built to use this option
#define CT_SUNDIALS_USE_LAPACK 0

// Enable export/import of HDF data via C++ HighFive
/* #undef CT_USE_HDF5 */
/* #undef CT_USE_SYSTEM_HIGHFIVE */
/* #undef CT_USE_HIGHFIVE_BOOLEAN */

#endif
