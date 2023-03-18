# get version of CoGAPS from DESCRIPTION file
AC_INIT(CoGAPS, m4_esyscmd_s([awk '/^Version:/ {print $2}' DESCRIPTION]))

# get C++ compiler from R configuration
CXX=`"${R_HOME}/bin/R" CMD config CXX`

# Switch to a C++ compiler, and check if it works.
AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX

# Check if compiling debug version
AC_ARG_ENABLE(debug, [AC_HELP_STRING([--enable-debug],
    [build debug version of CoGAPS])], [build_debug=yes], [build_debug=no])

# Check if running C++ tests (this adds to compilation time)
AC_ARG_ENABLE(cpp-tests, [AC_HELP_STRING([--enable-cpp-tests],
    [turn on C++ unit tests])], [cpp_tests=yes], [cpp_tests=no])

# Check if compiler warnings should be turned on
AC_ARG_ENABLE(warnings, [AC_HELP_STRING([--enable-warnings],
    [compile CoGAPS with warning messages])], [warnings=yes], [warnings=no])

# Use SIMD unless requested not to
AC_ARG_ENABLE(simd, [AC_HELP_STRING([--enable-simd],
    [compile with SIMD support if available])],
    [use_simd=$enableval], [use_simd=yes])

# default CoGAPS specific flags
GAPS_CPP_FLAGS=" -DBOOST_MATH_PROMOTE_DOUBLE_POLICY=0 -DGAPS_DISABLE_CHECKPOINTS -D__GAPS_R_BUILD__"
GAPS_CXX_FLAGS=
GAPS_LIBS=

# get compiler info
AX_COMPILER_VENDOR
AX_COMPILER_VERSION

# set openmp flags, disable only if requested
AC_ARG_ENABLE(openmp, [AC_HELP_STRING([--enable-openmp],
    [compile with openMP support if available])],
    [use_openmp=$enableval], [use_openmp=yes])

AX_OPENMP
if test "x$use_openmp" != "xno" ; then
    GAPS_CXX_FLAGS+=" $OPENMP_CXXFLAGS "
    GAPS_LIBS+=" $OPENMP_CXXFLAGS "
fi

echo "building on $ax_cv_cxx_compiler_vendor compiler version $ax_cv_cxx_compiler_version"

# set compile flags for debug build
if test "x$build_debug" = "xyes" ; then
    echo "Building Debug Version of CoGAPS"
    GAPS_CPP_FLAGS+=" -DGAPS_DEBUG "
fi

# set compile flags if warnings enabled
if test "x$warnings" = "xyes" ; then
    if test "x$ax_cv_cxx_compiler_vendor" = "xgnu" || test "x$ax_cv_cxx_compiler_vendor" = "xclang"; then
        GAPS_CXX_FLAGS+=" -Wall -Wextra -Werror -Wno-unused-parameter -Wno-unknown-pragmas"
    fi
fi

if test "x$use_simd" = "xno" ; then
    echo "Building without SIMD instructions"
fi

if test "x$use_simd" = "xsse" ; then
    echo "Disabling AVX instructions"
    GAPS_CXX_FLAGS+=" -mno-avx "
fi

GAPS_SOURCE_FILES+=" Cogaps.o"
GAPS_SOURCE_FILES+=" GapsParameters.o"
GAPS_SOURCE_FILES+=" GapsResult.o"
GAPS_SOURCE_FILES+=" GapsRunner.o"
GAPS_SOURCE_FILES+=" GapsStatistics.o"
GAPS_SOURCE_FILES+=" RcppExports.o"
GAPS_SOURCE_FILES+=" test-runner.o"
GAPS_SOURCE_FILES+=" atomic/Atom.o"
GAPS_SOURCE_FILES+=" atomic/ConcurrentAtom.o"
GAPS_SOURCE_FILES+=" atomic/AtomicDomain.o"
GAPS_SOURCE_FILES+=" atomic/ConcurrentAtomicDomain.o"
GAPS_SOURCE_FILES+=" atomic/ProposalQueue.o"
GAPS_SOURCE_FILES+=" data_structures/HashSets.o"
GAPS_SOURCE_FILES+=" data_structures/HybridMatrix.o"
GAPS_SOURCE_FILES+=" data_structures/HybridVector.o"
GAPS_SOURCE_FILES+=" data_structures/Matrix.o"
GAPS_SOURCE_FILES+=" data_structures/SparseIterator.o"
GAPS_SOURCE_FILES+=" data_structures/SparseMatrix.o"
GAPS_SOURCE_FILES+=" data_structures/SparseVector.o"
GAPS_SOURCE_FILES+=" data_structures/Vector.o"
GAPS_SOURCE_FILES+=" file_parser/CharacterDelimitedParser.o"
GAPS_SOURCE_FILES+=" file_parser/FileParser.o"
GAPS_SOURCE_FILES+=" file_parser/MatrixElement.o"
GAPS_SOURCE_FILES+=" file_parser/MtxParser.o"
GAPS_SOURCE_FILES+=" gibbs_sampler/AlphaParameters.o"
GAPS_SOURCE_FILES+=" gibbs_sampler/DenseNormalModel.o"
GAPS_SOURCE_FILES+=" gibbs_sampler/SparseNormalModel.o"
GAPS_SOURCE_FILES+=" math/Math.o"
GAPS_SOURCE_FILES+=" math/MatrixMath.o"
GAPS_SOURCE_FILES+=" math/Random.o"
GAPS_SOURCE_FILES+=" math/VectorMath.o"

# add c++ tests to source list
if test "x$cpp_tests" = "xyes" ; then
    echo "Enabling C++ Unit Tests"
    GAPS_CPP_FLAGS+=" -DGAPS_CPP_UNIT_TESTS "
    GAPS_SOURCE_FILES+=" cpp_tests/testVector.o"
fi

# export variables containing flags
AC_SUBST(GAPS_CPP_FLAGS)
AC_SUBST(GAPS_CXX_FLAGS)
AC_SUBST(GAPS_LIBS)
AC_SUBST(GAPS_SOURCE_FILES)

# create makefile, output configure file
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

