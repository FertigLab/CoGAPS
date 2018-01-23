#ifndef __COGAPS_ARCHIVE_H__
#define __COGAPS_ARCHIVE_H__

#include <fstream>

#define ARCHIVE_READ  std::ios::in
#define ARCHIVE_WRITE std::ios::out | std::ios::trunc

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
    friend Archive& operator<<(Archive &ar, T val);

    template<typename T>
    friend Archive& operator>>(Archive &ar, T &val);
};

template<typename T>
Archive& operator<<(Archive &ar, T val)
{
    ar.mStream.write(reinterpret_cast<char*>(&val), sizeof(T));
    return ar;
}

template<typename T>
Archive& operator>>(Archive &ar, T &val)
{
    ar.mStream.read(reinterpret_cast<char*>(&val), sizeof(T));
    return ar;
}

#endif