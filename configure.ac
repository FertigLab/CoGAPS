# get version of CoGAPS from DESCRIPTION file
AC_INIT(CoGAPS, m4_esyscmd_s([awk -e '/^Version:/ {print $2}' DESCRIPTION]))

# get C++ compiler from R configuration
CXX=`"${R_HOME}/bin/R" CMD config CXX`

# Switch to a C++ compiler, and check if it works.
AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX

# Check if compiling debug version
AC_ARG_ENABLE(debug, [AC_HELP_STRING([--enable-debug],
    [build debug version of CoGAPS])], [build_debug=yes], [build_debug=no])

# Check if compiler warnings should be turned on
AC_ARG_ENABLE(warnings, [AC_HELP_STRING([--enable-warnings],
    [compile CoGAPS with warning messages])], [warnings=yes], [warnings=no])

# Check if specific version of SIMD instructions requested
AC_ARG_ENABLE(simd, [AC_HELP_STRING([--enable-simd],
    [specify simd instruction set (sse, avx)])],
    [simd_version=$enableval], [simd_version=sse])

# default CoGAPS specific flags
GAPS_CPP_FLAGS=" -DBOOST_MATH_PROMOTE_DOUBLE_POLICY=0 -D__GAPS_R_BUILD__ "
GAPS_CXX_FLAGS=
GAPS_LIBS=

# get compiler info
AX_COMPILER_VENDOR
AX_COMPILER_VERSION

# set openmp flags
AX_OPENMP
GAPS_CXX_FLAGS+=" $OPENMP_CXXFLAGS "
GAPS_LIBS+=" $OPENMP_CXXFLAGS "

echo "building on $ax_cv_cxx_compiler_vendor compiler version $ax_cv_cxx_compiler_version"

# set compile flags for debug build
if test "x$build_debug" = "xyes" ; then
    echo "Building Debug Version of CoGAPS"
    GAPS_CPP_FLAGS+=" -DGAPS_DEBUG "
fi

# set compile flags if warnings enabled
if test "x$warnings" = "xyes" ; then
    if test "x$ax_cv_cxx_compiler_vendor" = "xgnu" || test "x$ax_cv_cxx_compiler_vendor" = "xclang"; then
        GAPS_CXX_FLAGS+=" -Wall -Wextra -Werror -Wno-unused-parameter "
    fi
fi

# set compile flags for SIMD
AX_EXT
echo "Supported SIMD flags: " $SIMD_FLAGS

if test "x$simd_version" != "xavx" && test "x$simd_version" != "xsse" && test "x$simd_version" != "xno" && test "x$simd_version" != "xyes" ; then
    echo "Invalid SIMD Type, Defaulting to AVX"
    simd_version=avx
fi

avx_supported=false
if test "x$simd_version" = "xavx" || test "x$simd_version" = "xyes" ; then
    for s in $SIMD_FLAGS; do
        if test "x$s" = "x-mavx"; then
            avx_supported=true
            echo "Building with AVX instructions"
            GAPS_CXX_FLAGS+=" -mavx "
        fi
    done
    if test "x$avx_supported" = "xfalse" ; then
        echo "Note: AVX not supported, trying for SSE"
        simd_version=sse
    fi
fi

sse_supported=false
if test "x$simd_version" = "xsse" ; then
    for s in $SIMD_FLAGS; do
        if test "x$s" = "x-msse4.2"; then
            sse_supported=true
            echo "Building with SSE instructions"
            GAPS_CXX_FLAGS+=" -msse4.2 "
        fi
    done
    if test "x$sse_supported" = "xfalse" ; then
        echo "Note: SSE not supported"
        simd_version=no
    fi
fi

if test "x$simd_version" = "xno" ; then
    echo "Building without SIMD instructions"
fi

# export variables containing flags
AC_SUBST(GAPS_CPP_FLAGS)
AC_SUBST(GAPS_CXX_FLAGS)
AC_SUBST(GAPS_LIBS)

# create makefile, output configure file
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

