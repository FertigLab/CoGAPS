# get version of CoGAPS from DESCRIPTION file
AC_INIT(CoGAPS, m4_esyscmd_s([awk -e '/^Version:/ {print $2}' Repo/DESCRIPTION]))

# get C++ compiler from R configuration
CXX=`"${R_HOME}/bin/R" CMD config CXX`

# Switch to a C++ compiler, and check if it works.
AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX

# macro to check if compile flag throws error
AC_DEFUN([AX_CHECK_COMPILE_FLAG],
[AC_PREREQ(2.64)dnl for _AC_LANG_PREFIX and AS_VAR_IF
AS_VAR_PUSHDEF([CACHEVAR],[ax_cv_check_[]_AC_LANG_ABBREV[]flags_$4_$1])dnl
AC_CACHE_CHECK([whether _AC_LANG compiler accepts $1], CACHEVAR, [
    ax_check_save_flags=$[]_AC_LANG_PREFIX[]FLAGS
    _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $4 $1"
    AC_COMPILE_IFELSE([m4_default([$5],[AC_LANG_PROGRAM()])],
        [AS_VAR_SET(CACHEVAR,[yes])],
        [AS_VAR_SET(CACHEVAR,[no])])
    _AC_LANG_PREFIX[]FLAGS=$ax_check_save_flags])
AS_VAR_IF(CACHEVAR,yes,
    [m4_default([$2], :)],
    [m4_default([$3], :)])
AS_VAR_POPDEF([CACHEVAR])dnl
])dnl AX_CHECK_COMPILE_FLAGS

# Check if compiling debug version
AC_ARG_ENABLE(debug, [AC_HELP_STRING([--enable-debug],
    [build debug version of CoGAPS])], [build_debug=yes], [build_debug=no])

# Check if specific version of SIMD instructions requested
AC_ARG_ENABLE(simd, [AC_HELP_STRING([--enable-simd],
    [specify simd instruction set (sse, avx, avx512)])],
    [simd_version=$enableval], [simd_version=sse])

GAPS_CPP_FLAGS=
GAPS_CXX_FLAGS=

if test "x$build_debug" = "xyes" ; then
    echo "Building Debug Version of CoGAPS"
    GAPS_CPP_FLAGS+=" -DGAPS_DEBUG "
fi

if test "x$simd_version" = "xsse" ; then
    echo "Compiling with SSE instructions"
    AX_CHECK_COMPILE_FLAG([-msse4.2], [GAPS_CXX_FLAGS+=" -msse4.2 "])
elif test "x$simd_version" = "xyes" ; then
    echo "Compiling with SSE instructions"
    AX_CHECK_COMPILE_FLAG([-msse4.2], [GAPS_CXX_FLAGS+=" -msse4.2 "])
elif test "x$simd_version" = "xavx" ; then
    echo "Compiling with AVX instructions"
    AX_CHECK_COMPILE_FLAG([-mavx], [GAPS_CXX_FLAGS+=" -mavx "])    
elif test "x$simd_version" = "xavx512" ; then
    echo "Compiling with AVX-512 instructions"
    AX_CHECK_COMPILE_FLAG([-mavx512f], [GAPS_CXX_FLAGS+=" -mavx512f "])
elif test "x$simd_version" = "xno" ; then
    echo "Compiling without SIMD instructions"
else
    echo "Error: Invalid SIMD type"
    exit -1
fi

AC_SUBST(GAPS_CPP_FLAGS)
AC_SUBST(GAPS_CXX_FLAGS)

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

