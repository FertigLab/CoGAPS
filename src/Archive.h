#ifndef __COGAPS_ARCHIVE_H__
#define __COGAPS_ARCHIVE_H__

#include <fstream>

#include <boost/random/mersenne_twister.hpp>

// flags for opening an archive
#define ARCHIVE_READ  std::ios::in
#define ARCHIVE_WRITE (std::ios::out | std::ios::trunc)

// magic number written to beginning of archive files
// needs to be updated everytime to method of checkpointing changes
#define ARCHIVE_MAGIC_NUM 0xCE45D32A

class Archive
{
private:

    std::fstream mStream;

public:

    Archive(const std::string &path, std::ios_base::openmode flags)
        : mStream(path.c_str(), std::ios::binary | flags)
    {}

    void close() {mStream.close();}

    template<typename T>
    friend Archive& writeToArchive(Archive &ar, T val)
    {
        ar.mStream.write(reinterpret_cast<char*>(&val), sizeof(T)); // NOLINT
        return ar;
    }

    template<typename T>
    friend Archive& readFromArchive(Archive &ar, T &val)
    {
        ar.mStream.read(reinterpret_cast<char*>(&val), sizeof(T)); // NOLINT
        return ar;
    }    

    // explicitly define which types can be automatically written/read
    // don't have C++11 and don't want to add another dependency on boost,
    // so no template tricks

    friend Archive& operator<<(Archive &ar, char val)     { return writeToArchive(ar, val); }
    friend Archive& operator<<(Archive &ar, bool val)     { return writeToArchive(ar, val); }
    friend Archive& operator<<(Archive &ar, int val)      { return writeToArchive(ar, val); }
    friend Archive& operator<<(Archive &ar, unsigned val) { return writeToArchive(ar, val); }
    friend Archive& operator<<(Archive &ar, uint64_t val) { return writeToArchive(ar, val); }
    friend Archive& operator<<(Archive &ar, int64_t val)  { return writeToArchive(ar, val); }
    friend Archive& operator<<(Archive &ar, float val)    { return writeToArchive(ar, val); }
    friend Archive& operator<<(Archive &ar, double val)   { return writeToArchive(ar, val); }
    friend Archive& operator<<(Archive &ar, boost::random::mt11213b val) { return writeToArchive(ar, val); }

    friend Archive& operator>>(Archive &ar, char &val)     { return readFromArchive(ar, val); }
    friend Archive& operator>>(Archive &ar, bool &val)     { return readFromArchive(ar, val); }
    friend Archive& operator>>(Archive &ar, int &val)      { return readFromArchive(ar, val); }
    friend Archive& operator>>(Archive &ar, unsigned &val) { return readFromArchive(ar, val); }
    friend Archive& operator>>(Archive &ar, uint64_t &val) { return readFromArchive(ar, val); }
    friend Archive& operator>>(Archive &ar, int64_t &val)  { return readFromArchive(ar, val); }
    friend Archive& operator>>(Archive &ar, float &val)    { return readFromArchive(ar, val); }
    friend Archive& operator>>(Archive &ar, double &val)   { return readFromArchive(ar, val); }
    friend Archive& operator>>(Archive &ar, boost::random::mt11213b &val) { return readFromArchive(ar, val); }
};

#endif // __COGAPS_ARCHIVE_H__
