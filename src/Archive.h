#ifndef __COGAPS_ARCHIVE_H__
#define __COGAPS_ARCHIVE_H__

#include <fstream>

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
    {
        /*if (flags == ARCHIVE_WRITE)
        {
            *this << ARCHIVE_MAGIC_NUM;
        }
        else if (flags == ARCHIVE_READ)
        {
            uint32_t magicNum = 0;
            *this >> magicNum;
            if (magicNum != ARCHIVE_MAGIC_NUM)
            {
                Rcpp::Rcout << "warning: invalid checkpoint file" << std::endl;
            }
        }*/
    }

    void close() {mStream.close();}

    template<typename T>
    friend Archive& operator<<(Archive &ar, T val)
    {
        ar.mStream.write(reinterpret_cast<char*>(&val), sizeof(T)); // NOLINT
        return ar;
    }

    template<typename T>
    friend Archive& operator>>(Archive &ar, T &val)
    {
        ar.mStream.read(reinterpret_cast<char*>(&val), sizeof(T)); // NOLINT
        return ar;
    }
};

#endif
