# include macros
m4_builtin([include], m4/ax_check_compile_flag.m4)
m4_builtin([include], m4/ax_compiler_vendor.m4)
m4_builtin([include], m4/ax_compiler_version.m4)
m4_builtin([include], m4/ax_openmp.m4)

# get version of CoGAPS from DESCRIPTION file
AC_INIT(CoGAPS, m4_esyscmd_s([awk -e '/^Version:/ {print $2}' Repo/DESCRIPTION]))

# get C++ compiler from R configuration
CXX=`"${R_HOME}/bin/R" CMD config CXX`

# Switch to a C++ compiler, and check if it works.
AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX

# Check if compiling debug version
AC_ARG_ENABLE(debug, [AC_HELP_STRING([--enable-debug],
    [build debug version of CoGAPS])], [build_debug=yes], [build_debug=no])

# Check if compiling debug version with thread sanitizer
AC_ARG_ENABLE(debug_tsan, [AC_HELP_STRING([--enable-debug-tsan],
    [build debug version of CoGAPS])], [build_debug_tsan=yes], [build_debug_tsan=no])

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
AX_OPENMP

# set openmp flags
GAPS_CXX_FLAGS+=" $OPENMP_CXXFLAGS "
GAPS_LIBS+=" $OPENMP_CXXFLAGS "

echo "building on $ax_cv_cxx_compiler_vendor compiler version $ax_cv_cxx_compiler_version"

# set compile flags for debug build
if test "x$build_debug" = "xyes" ; then
    echo "Building Debug Version of CoGAPS"
    GAPS_CPP_FLAGS+=" -DGAPS_DEBUG -fno-omit-frame-pointer "
    if test "x$ax_cv_cxx_compiler_vendor" = "xclang" ; then
        GAPS_CPP_FLAGS+=" -fsanitize=address -fsanitize=undefined "
        GAPS_LIBS+=" -fsanitize=address -fsanitize=undefined "
    elif test "x$ax_cv_cxx_compiler_vendor" = "xgnu" ; then
        GAPS_CPP_FLAGS+=" -fsanitize=address -fsanitize=undefined "
        GAPS_LIBS+=" -fsanitize=address -fsanitize=undefined "
    fi
elif test "x$build_debug_tsan" = "xyes" ; then
    echo "Building Debug Version of CoGAPS with thread sanitizer"
    GAPS_CPP_FLAGS+=" -DGAPS_DEBUG -fno-omit-frame-pointer "
    if test "x$ax_cv_cxx_compiler_vendor" = "xclang" ; then
        GAPS_CPP_FLAGS+=" -fsanitize=thread "
        GAPS_LIBS+=" -fsanitize=thread "
    fi
fi

# set compile flags if warnings enabled
if test "x$warnings" = "xyes" ; then
    AX_CHECK_COMPILE_FLAG([-Wall], [GAPS_CXX_FLAGS+=" -Wall "])
    AX_CHECK_COMPILE_FLAG([-Wextra], [GAPS_CXX_FLAGS+=" -Wextra "])
    AX_CHECK_COMPILE_FLAG([-Werror], [GAPS_CXX_FLAGS+=" -Werror "])
    AX_CHECK_COMPILE_FLAG([-Wno-unused-parameter], [GAPS_CXX_FLAGS+=" -Wno-unused-parameter "])
    AX_CHECK_COMPILE_FLAG([-Wno-unused-function], [GAPS_CXX_FLAGS+=" -Wno-unused-function "])
fi

# set compile flags for SIMD
if test "x$simd_version" = "xsse" ; then
    echo "Compiling with SSE instructions"
    AX_CHECK_COMPILE_FLAG([-msse4.2], [GAPS_CXX_FLAGS+=" -msse4.2 "])
elif test "x$simd_version" = "xyes" ; then
    echo "Compiling with SSE instructions"
    AX_CHECK_COMPILE_FLAG([-msse4.2], [GAPS_CXX_FLAGS+=" -msse4.2 "])
elif test "x$simd_version" = "xavx" ; then
    echo "Compiling with AVX instructions"
    AX_CHECK_COMPILE_FLAG([-mavx], [GAPS_CXX_FLAGS+=" -mavx "])    
elif test "x$simd_version" = "xno" ; then
    echo "Compiling without SIMD instructions"
else
    echo "Error: Invalid SIMD type"
    exit -1
fi

# export variables containing flags
AC_SUBST(GAPS_CPP_FLAGS)
AC_SUBST(GAPS_CXX_FLAGS)
AC_SUBST(GAPS_LIBS)

# create makefile, output configure file
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

