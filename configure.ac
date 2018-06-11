# include macros
m4_builtin([include], ../ax_check_compile_flag.m4)
m4_builtin([include], ../ax_compiler_vendor.m4)
m4_builtin([include], ../ax_compiler_version.m4)
m4_builtin([include], ../ax_openmp.m4)
m4_builtin([include], ../ax_ext.m4)
m4_builtin([include], ../ax_gcc_x86_cpuid.m4)
m4_builtin([include], ../ax_gcc_x86_avx_xgetbv.m4)

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

# Check if specific version of SIMD instructions requested
AC_ARG_ENABLE(simd, [AC_HELP_STRING([--enable-simd],
    [specify simd instruction set (sse, avx)])],
    [simd_version=$enableval], [simd_version=sse])

# default CoGAPS specific flags
GAPS_CPP_FLAGS=
GAPS_CXX_FLAGS=

# get compiler info
AX_COMPILER_VENDOR
AX_COMPILER_VERSION
AX_OPENMP
AX_EXT

echo $SIMD_FLAGS
echo $CPUEXT_FLAGS
echo $OPENMP_CXXFLAGS
echo $ax_cv_cxx_compiler_vendor
echo $ax_cv_cxx_compiler_version

if test "x$build_debug" = "xyes" ; then
    echo "Building Debug Version of CoGAPS"
    GAPS_CPP_FLAGS+=" -DGAPS_DEBUG "
    GAPS_CXX_FLAGS+=" -Wall -Wextra -Werror "
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

