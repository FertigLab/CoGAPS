#include "catch.h"
#include "../AtomicSupport.h"

TEST_CASE("Test AtomicSupport.h")
{
    unsigned nrow = 100, ncol = 500;
    uint64_t maxAtoms = 100;

    SECTION("Make and Convert proposals")
    {
        AtomicSupport domain('A', nrow, ncol, 1.0, 0.5);
        domain.setMaxNumAtoms(maxAtoms);
        REQUIRE(domain.alpha() == 1.0);
        REQUIRE(domain.lambda() == 0.5);
        
        gaps::random::setSeed(1);
        AtomicProposal prop = domain.makeProposal();
        
        REQUIRE(prop.label == 'A');
        REQUIRE(prop.type == 'B');
        REQUIRE(prop.nChanges == 1);
        REQUIRE(prop.pos1 == 0xf21db672ab2f52a4);
        REQUIRE(prop.delta1 == 1.0);
        REQUIRE(prop.pos2 == 0);
        REQUIRE(prop.delta2 == 0.0);
        
        domain.acceptProposal(prop);

        REQUIRE(domain.numAtoms() == 1);
        REQUIRE(domain.totalMass() == prop.delta1);

        for (unsigned i = 0; i < 10000; ++i)
        {
            prop = domain.makeProposal();

            REQUIRE(prop.label == 'A');
            bool cond1 = prop.type == 'B' && prop.nChanges == 1;
            bool cond2 = prop.type == 'D' && prop.nChanges == 1;
            bool cond3 = prop.type == 'M' && prop.nChanges == 2;
            bool cond4 = prop.type == 'E' && prop.nChanges == 2;
            bool cond = cond1 || cond2 || cond3 || cond4;

            REQUIRE(cond);
            
            MatrixChange change = domain.getMatrixChange(prop);
            
            REQUIRE(change.label == 'A');
            REQUIRE(change.nChanges == prop.nChanges);
            cond = change.row1 >= 0 && change.row2 < nrow;
            REQUIRE(cond);
            cond = change.col1 >= 0 && change.col2 < ncol;
            REQUIRE(cond);
            REQUIRE(change.delta1 == prop.delta1);
            REQUIRE(change.delta2 == prop.delta2);

            double oldMass = domain.totalMass();

            domain.acceptProposal(prop);
        
            REQUIRE(domain.totalMass() == oldMass + prop.delta1 + prop.delta2);
            REQUIRE(domain.numAtoms() <= maxAtoms);
        }
    }
}