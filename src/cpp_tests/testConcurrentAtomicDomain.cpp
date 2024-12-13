#include <testthat.h>
#include "../testthat-tweak.h"
#include "../atomic/ConcurrentAtomicDomain.h"
#include "../math/Random.h"
#include "../utils/GapsPrint.h"

TEST_CASE("AtomicDomain")
{
    GapsRandomState randState(123);
    GapsRng rng(&randState);

    SECTION("Construction")
    {
        ConcurrentAtomicDomain domain(10);
        REQUIRE(domain.size() == 0);
    }
}
