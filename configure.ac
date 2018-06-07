AC_INIT(CoGAPS, m4_esyscmd_s([awk -e '/^Version:/ {print $2}' Repo/DESCRIPTION]))

dnl Switch to a C++ compiler, and check if it works.
AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX

dnl Check if compiling debug version
AC_ARG_ENABLE(debug, [AC_HELP_STRING([--enable-debug],
    [build debug version of CoGAPS])], [build_debug=yes], [build_debug=no])

dnl Check if specific version of SIMD instructions requested
AC_ARG_ENABLE(simd, [AC_HELP_STRING([--enable-simd], [specify simd instruction set (sse, avx, avx512)])], [simd_version=$enableval], [simd_version=sse])

dnl Check if compiler supports openmp
dnl TODO

GAPS_CPP_FLAGS=
GAPS_CXX_FLAGS=

if test "x$build_debug" = "xyes" ; then
    echo "Building Debug Version of CoGAPS"
    GAPS_CPP_FLAGS+=" -DGAPS_DEBUG "
fi

if test "x$simd_version" = "xsse" ; then
    echo "Using SSE"
    GAPS_CXX_FLAGS+=" -msse4.2 "
elif test "x$simd_version" = "xavx" ; then
    echo "Using AVX"
    GAPS_CXX_FLAGS+=" -mavx "
elif test "x$simd_version" = "xavx512" ; then
    echo "Using AVX-512"
    GAPS_CXX_FLAGS+=" -mavx512f "
else
    echo "Error: Invalid SIMD type"
    exit -1
fi

AC_SUBST(GAPS_CPP_FLAGS)
AC_SUBST(GAPS_CXX_FLAGS)

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

