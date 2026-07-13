#include <testthat.h>
#include "../testthat-tweak.h"
#include "../gibbs_sampler/SingleThreadedGibbsSampler.h"
#include "../gibbs_sampler/SparseNormalModel.h"
#include "../gibbs_sampler/DenseNormalModel.h"
#include "../gibbs_sampler/AlphaParameters.h"
#include "../math/Math.h"
#include "../math/MatrixMath.h"

#include <cmath>

#define TEST_APPROX(x) Approx(x).epsilon(0.001f)

// Test-only subclass that exposes the protected alphaParameters()/mMatrix of the
// DataModel, so the dense-vs-sparse consistency can be checked without changing
// production access levels.
template <class DataModel>
class ExposedSampler : public SingleThreadedGibbsSampler<DataModel>
{
public:
    template <class DataType>
    ExposedSampler(const DataType &data, bool transpose, bool subsetRows,
        float alpha, float maxGibbsMass, const GapsParameters &params,
        GapsRandomState *randState)
    : SingleThreadedGibbsSampler<DataModel>(data, transpose, subsetRows, alpha,
        maxGibbsMass, params, randState) {}

    AlphaParameters a1(unsigned r, unsigned c)
        { return this->alphaParameters(r, c); }
    AlphaParameters a2(unsigned r1, unsigned c1, unsigned r2, unsigned c2)
        { return this->alphaParameters(r1, c1, r2, c2); }
    AlphaParameters aC(unsigned r, unsigned c, float ch)
        { return this->alphaParametersWithChange(r, c, ch); }
    float matrixSum() const { return gaps::sum(this->mMatrix); }
};

static void requireAlphaEqual(AlphaParameters sa, AlphaParameters da)
{
    REQUIRE(sa.s >= 0.f);
    REQUIRE(da.s >= 0.f);
    if (sa.s <= gaps::epsilon || da.s <= gaps::epsilon)
    {
        REQUIRE(sa.s <= gaps::epsilon);
        REQUIRE(da.s <= gaps::epsilon);
    }
    REQUIRE(sa.s == TEST_APPROX(da.s));
    REQUIRE(sa.s_mu == TEST_APPROX(da.s_mu));
}

TEST_CASE("Test SparseGibbsSampler", "[sparsegibbs]")
{
    SECTION("Construct from data matrix")
    {
        Matrix data(25, 50);
        for (unsigned i = 0; i < data.nRow(); ++i)
        {
            for (unsigned j = 0; j < data.nCol(); ++j)
            {
                data(i,j) = i + j + 1.f;
            }
        }

        GapsRandomState randState(123);
        GapsParameters params(data);
        SingleThreadedGibbsSampler<SparseNormalModel> ASampler(data, true, false, params.alphaA,
            params.maxGibbsMassA, params, &randState);
        SingleThreadedGibbsSampler<SparseNormalModel> PSampler(data, false, false, params.alphaP,
            params.maxGibbsMassP, params, &randState);

        ASampler.sync(PSampler);
        PSampler.sync(ASampler);

        REQUIRE(ASampler.chiSq() == 100.f * data.nRow() * data.nCol());
        REQUIRE(PSampler.chiSq() == 100.f * data.nRow() * data.nCol());

    #ifdef GAPS_DEBUG
        REQUIRE(ASampler.internallyConsistent());
        REQUIRE(PSampler.internallyConsistent());
    #endif
    }

    SECTION("Update decreases chiSq")
    {
        Matrix data(25, 50);
        for (unsigned i = 0; i < data.nRow(); ++i)
            for (unsigned j = 0; j < data.nCol(); ++j)
                data(i,j) = i + j + 1.f;

        GapsRandomState randState(42);
        GapsParameters params(data);
        SingleThreadedGibbsSampler<SparseNormalModel> ASampler(data, true, false, params.alphaA,
            params.maxGibbsMassA, params, &randState);
        SingleThreadedGibbsSampler<SparseNormalModel> PSampler(data, false, false, params.alphaP,
            params.maxGibbsMassP, params, &randState);

        // chiSq() uses mOtherMatrix, which is only valid after sync()
        ASampler.sync(PSampler);
        PSampler.sync(ASampler);
        ASampler.extraInitialization();
        PSampler.extraInitialization();

        double AChiInit = ASampler.chiSq();
        double PChiInit = PSampler.chiSq();

        // interleave update+sync as the algorithm does: each sampler re-syncs to
        // the other's updated matrix before updating (sparse rebuilds its lookup
        // tables in sync(), so it must sync after the other matrix changes)
        ASampler.update(100);
        PSampler.sync(ASampler);
        PSampler.update(100);
        ASampler.sync(PSampler);
        ASampler.extraInitialization();
        PSampler.extraInitialization();

        // sampling improved the fit
        REQUIRE(ASampler.chiSq() < AChiInit);
        REQUIRE(PSampler.chiSq() < PChiInit);
    }

    SECTION("chiSq() before sync() does not crash and matches dense (regression)")
    {
        // Regression for issue #18: SparseNormalModel::chiSq() dereferenced the
        // NULL mOtherMatrix when called before sync(), segfaulting -- whereas
        // DenseNormalModel::chiSq() is safe in the same (un-synced) state. Both
        // must return the "no fit" chiSq (A*P treated as 0), matching each other.
        Matrix data(25, 50);
        for (unsigned i = 0; i < data.nRow(); ++i)
            for (unsigned j = 0; j < data.nCol(); ++j)
                data(i,j) = i + j + 1.f; // all non-zero, all >= 1 -> S = 0.1*D

        GapsRandomState randState(123);
        GapsParameters params(data);
        SingleThreadedGibbsSampler<SparseNormalModel> sparse(data, true, false,
            params.alphaA, params.maxGibbsMassA, params, &randState);
        SingleThreadedGibbsSampler<DenseNormalModel> dense(data, true, false,
            params.alphaA, params.maxGibbsMassA, params, &randState);

        // no sync() -> mOtherMatrix is NULL. Must not crash.
        float sparseChi = sparse.chiSq();
        float denseChi = dense.chiSq();

        REQUIRE(std::isfinite(sparseChi));
        // no fit: (D-0)^2/S^2 with S = 0.1*D is 100 per entry, over nRow*nCol entries
        REQUIRE(sparseChi == 100.f * data.nRow() * data.nCol());
        REQUIRE(sparseChi == denseChi);
    }

    SECTION("Test consistency between alpha parameters calculations")
    {
        // create the "data"
        Matrix data(100, 75);
        GapsRandomState randState(123);
        GapsRng rng(&randState);
        for (unsigned i = 0; i < data.nRow(); ++i)
        {
            for (unsigned j = 0; j < data.nCol(); ++j)
            {
                // continuous values include 0 < D < 1 (exercise the S floor) and
                // D > 1 (S = factor*D regime), with ~50% zeros
                data(i,j) = rng.uniform(0.f, 3.f) * (rng.uniform() < 0.5f ? 0.f : 1.f);
            }
        }

        GapsParameters params(data);

        // pair of sparse and pair of dense samplers over the same data
        ExposedSampler<SparseNormalModel> sparse_A(data, true, false, params.alphaA,
            params.maxGibbsMassA, params, &randState);
        ExposedSampler<SparseNormalModel> sparse_P(data, false, false, params.alphaP,
            params.maxGibbsMassP, params, &randState);
        sparse_A.sync(sparse_P);
        sparse_P.sync(sparse_A);

        ExposedSampler<DenseNormalModel> dense_A(data, true, false, params.alphaA,
            params.maxGibbsMassA, params, &randState);
        ExposedSampler<DenseNormalModel> dense_P(data, false, false, params.alphaP,
            params.maxGibbsMassP, params, &randState);
        dense_A.sync(dense_P);
        dense_P.sync(dense_A);

        // set the A matrix (genes x patterns) to the same thing on both samplers
        Matrix AMat(data.nRow(), params.nPatterns);
        for (unsigned i = 0; i < data.nRow(); ++i)
            for (unsigned k = 0; k < params.nPatterns; ++k)
                AMat(i,k) = rng.uniform(0.f, 10.f) * (rng.uniform() < 0.2f ? 0.f : 1.f);
        sparse_A.setMatrix(AMat);
        dense_A.setMatrix(AMat);
        REQUIRE(dense_A.matrixSum() == sparse_A.matrixSum());

        // and the P matrix (samples x patterns)
        Matrix PMat(data.nCol(), params.nPatterns);
        for (unsigned j = 0; j < data.nCol(); ++j)
            for (unsigned k = 0; k < params.nPatterns; ++k)
                PMat(j,k) = rng.uniform(0.f, 10.f) * (rng.uniform() < 0.2f ? 0.f : 1.f);
        sparse_P.setMatrix(PMat);
        dense_P.setMatrix(PMat);
        REQUIRE(dense_P.matrixSum() == sparse_P.matrixSum());

        // sync them back up (dense needs extraInitialization to build the AP matrix;
        // sparse's is a nop, its sync regenerates the lookup tables)
        sparse_A.sync(sparse_P);
        sparse_P.sync(sparse_A);
        dense_A.sync(dense_P);
        dense_P.sync(dense_A);
        dense_A.extraInitialization();
        dense_P.extraInitialization();

        // chiSq must also match between sparse and dense. This exercises the S
        // floor in chiSq() specifically: the data has 0 < D < 1 entries where the
        // unfloored S = D disagrees with the floored S = max(0.1*D, 0.1) (issue #18
        // follow-up). Both must use the single invSSq() model.
        REQUIRE(sparse_A.chiSq() == TEST_APPROX(dense_A.chiSq()));
        REQUIRE(sparse_P.chiSq() == TEST_APPROX(dense_P.chiSq()));

        // 1D alphaParameters must match between sparse and dense
        for (unsigned i = 0; i < data.nRow(); ++i)
            for (unsigned k = 0; k < params.nPatterns; ++k)
                requireAlphaEqual(sparse_A.a1(i,k), dense_A.a1(i,k));
        for (unsigned j = 0; j < data.nCol(); ++j)
            for (unsigned k = 0; k < params.nPatterns; ++k)
                requireAlphaEqual(sparse_P.a1(j,k), dense_P.a1(j,k));

        // 2D alphaParameters (and symmetry) must match
        for (unsigned i = 0; i < data.nRow(); ++i)
            for (unsigned k1 = 0; k1 < params.nPatterns; ++k1)
                for (unsigned k2 = k1+1; k2 < params.nPatterns; ++k2)
                {
                    requireAlphaEqual(sparse_A.a2(i,k1,i,k2), dense_A.a2(i,k1,i,k2));
                    requireAlphaEqual(sparse_A.a2(i,k2,i,k1), dense_A.a2(i,k2,i,k1));
                }
        for (unsigned j = 0; j < data.nCol(); ++j)
            for (unsigned k1 = 0; k1 < params.nPatterns; ++k1)
                for (unsigned k2 = k1+1; k2 < params.nPatterns; ++k2)
                {
                    requireAlphaEqual(sparse_P.a2(j,k1,j,k2), dense_P.a2(j,k1,j,k2));
                    requireAlphaEqual(sparse_P.a2(j,k2,j,k1), dense_P.a2(j,k2,j,k1));
                }

        // alphaParametersWithChange must match
        for (unsigned i = 0; i < data.nRow(); ++i)
            for (unsigned k = 0; k < params.nPatterns; ++k)
            {
                float ch = rng.uniform(0.f, 25.f);
                requireAlphaEqual(sparse_A.aC(i,k,ch), dense_A.aC(i,k,ch));
            }
        for (unsigned j = 0; j < data.nCol(); ++j)
            for (unsigned k = 0; k < params.nPatterns; ++k)
            {
                float ch = rng.uniform(0.f, 25.f);
                requireAlphaEqual(sparse_P.aC(j,k,ch), dense_P.aC(j,k,ch));
            }
    }
}
