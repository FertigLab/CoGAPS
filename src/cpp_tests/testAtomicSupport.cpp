#include "catch.h"
#include "../AtomicSupport.h"

TEST_CASE("Test AtomicSupport.h")
{
    unsigned nrow = 100, ncol = 500;

    SECTION("Make and Convert proposals")
    {
        AtomicSupport domain('A', nrow, ncol, 1.0, 0.5);
        REQUIRE(domain.alpha() == 1.0);
        REQUIRE(domain.lambda() == 0.5);
        
        gaps::random::setSeed(1);
        AtomicProposal prop = domain.makeProposal();
        
        REQUIRE(prop.label == 'A');
        REQUIRE(prop.type == 'B');
        REQUIRE(prop.nChanges == 1);
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

            uint64_t nOld = domain.numAtoms();

            domain.acceptProposal(prop);

            if (prop.type == 'B')
            {
                REQUIRE(domain.numAtoms() == nOld + 1);
            }
            else if (prop.type == 'D')
            {
                REQUIRE(domain.numAtoms() == nOld - 1);
            }
            /*else
            {
                REQUIRE(domain.numAtoms() == nOld);
            }*/
        
            REQUIRE(domain.totalMass() == oldMass + prop.delta1 + prop.delta2);
        }
    }

    SECTION("Proposal Distribution")
    {
        AtomicSupport domain('A', nrow, ncol, 0.01, 0.05);
        
        gaps::random::setSeed(1);

        unsigned nB = 0, nD = 0, nM = 0, nE = 0;
        for (unsigned i = 0; i < 5000; ++i)
        {
            AtomicProposal prop = domain.makeProposal();
            domain.acceptProposal(prop);

            switch(prop.type)
            {
                case 'B': nB++; break;
                case 'D': nD++; break;
                case 'M': nM++; break;
                case 'E': nE++; break;
            }
        }
        REQUIRE(domain.numAtoms() > 100);
        REQUIRE(nB > 500);
        //REQUIRE(nD > 500);
        REQUIRE(nM > 500);
        REQUIRE(nE > 500);
    }
}