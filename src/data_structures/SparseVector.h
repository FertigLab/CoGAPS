// can access without iterator, can set elements with accesor
class Vector
{
public:

    explicit Vector(unsigned size);
    explicit Vector(const std::vector<float> &v);

    float& operator()(unsigned i, unsigned j); // set value
    const float* ptr() const; // access without iterator

    friend Archive& operator<<(Archive &ar, Vector &vec);
    friend Archive& operator>>(Archive &ar, Vector &vec);

private:

    aligned_vector mData;
};

// can only access through iterator, all data is const
class SparseVector
{
public:

    explicit SparseVector(unsigned size);
    explicit SparseVector(const std::vector<float> &v);

    friend Archive& operator<<(Archive &ar, Vector &vec);
    friend Archive& operator>>(Archive &ar, Vector &vec);

private:
    
    std::vector<uint64_t> mIndexBitFlags;
    std::vector<float> mData;
};

// stored as a dense vector (efficient setting of values) but maintains
// index bit flags of non-zeros so it can be used with SparseIterator
class HybridVector
{
public:

    explicit HybridVector(unsigned size);
    explicit HybridVector(const std::vector<float> &v);

    void change(unsigned i, unsigned j, float v);

    friend Archive& operator<<(Archive &ar, Vector &vec);
    friend Archive& operator>>(Archive &ar, Vector &vec);

private:

    std::vector<uint64_t> mIndexBitFlags;
    std::vector<float> mData;
};

class SparseIterator
{
public:

    SparseIterator(const HybridVector &A, const SparseVector &B);

    bool atEnd() const;
    void next();
    float firstValue();
    float secondValue();
    unsigned firstIndex();
    unsigned secondIndex();
};

