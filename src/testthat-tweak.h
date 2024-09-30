//(c) Alexander Favorov favorov@sensi.org
//MIT Artistic license
//We want to use testthat (written for version "3.2.1.1")
//and catch via testthat and to use tags to run tests separately

#ifndef TESTTHAT_TWEAK_HPP
#define TESTTHAT_TWEAK_HPP

//add this after #include <testthat.h>
//if you want to use context/test_that/expect_true style and tags
#undef context
#define context(__X__,__Y__) CATCH_TEST_CASE(__X__,__Y__)

//AF: add this after #include <testthat.h>
//if you want to use TEST_CASE/SECTION/CHECK style
#define TEST_CASE CATCH_TEST_CASE
#define SECTION CATCH_SECTION 
#define CHECK CATCH_CHECK
#define REQUIRE CATCH_CHECK

#endif //TESTTHAT_TWEAK_HPP