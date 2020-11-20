//
//  utils.hpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#ifndef utils_hpp
#define utils_hpp

#include <string>
#include <vector>
#include <mutex>
#include <set>
#include <limits>
#include <functional>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <iostream>
#include <fstream>

#include <unistd.h>

#include <spdlog/spdlog.h>
#include <mio/mmap.hpp>
#include <omp.h>

namespace vcfbwt
{

typedef std::uint64_t long_type;
typedef std::uint32_t short_type;
static long_type  long_prime   = 27162335252586509;
static short_type short_prime  =        1999999973;

#ifdef PFP_LONG_TYPE
//#pragma message("Using 64 bit parse")
typedef long_type  size_type;
#else
//#pragma message("Using 32 bit parse")
typedef short_type size_type;
#endif
typedef long_type  hash_type;

constexpr std::size_t KILOBYTE     = 1024;
constexpr std::size_t MEGABYTE     = KILOBYTE * KILOBYTE;
constexpr std::size_t GIGABYTE     = KILOBYTE * MEGABYTE;

constexpr double KILOBYTE_DOUBLE = 1024.0;
constexpr double MILLION_DOUBLE  = 1000000.0;
constexpr double MEGABYTE_DOUBLE = KILOBYTE_DOUBLE * KILOBYTE_DOUBLE;
constexpr double GIGABYTE_DOUBLE = KILOBYTE_DOUBLE * MEGABYTE_DOUBLE;

constexpr std::size_t MILLION      = 1000000;
constexpr std::size_t BILLION      = 1000 * MILLION;


//------------------------------------------------------------------------------

inline double
inMegabytes(std::size_t bytes)
{
    return bytes / MEGABYTE_DOUBLE;
}

inline double
inGigabytes(std::size_t bytes)
{
    return bytes / GIGABYTE_DOUBLE;
}

inline std::size_t pid()
{
#ifdef MSVC_COMPILER
    return _getpid();
#else
    return getpid();
#endif
}

//------------------------------------------------------------------------------

inline void
set_prime(std::size_t p)
{
    if (sizeof(hash_type) > 4) { long_prime = p; }
    else if (p < std::numeric_limits<hash_type>::max()) { short_prime = p; }
    else { spdlog::error("Prime is too big for {} bytes", sizeof(hash_type)); std::exit(EXIT_FAILURE); }
}

/*
  Thomas Wang's integer hash function. In many implementations, std::hash
  is identity function for integers, which leads to performance issues.
*/

inline std::size_t
wang_hash_64(std::size_t key)
{
    key = (~key) + (key << 21); // key = (key << 21) - key - 1;
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8); // key * 265
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4); // key * 21
    key = key ^ (key >> 28);
    key = key + (key << 31);
    return key;
}

inline hash_type
string_hash(const char* s, std::size_t size)
{
    hash_type hash = 0;
    hash_type prime;
    if (sizeof(hash_type) > 4) { prime = long_prime; }
    else { prime = short_prime; }
    for(std::size_t k = 0; k < size; k++)
    {
        int c = (unsigned char) s[k];
        assert(c >= 0 && c < 256);
        hash = (256 * hash + c) % prime;    //  add char k
    }
    return hash;
}

/*
  Window in a string and its KR fingerprint
*/
struct KR_window {
    int wsize;
    int *window;
    int asize;
    hash_type prime;
    hash_type hash;
    size_type tot_char;
    size_type asize_pot;   // asize^(wsize-1) mod prime
    
    KR_window(int w): wsize(w) {
        asize = 256;
        asize_pot = 1;
        if (sizeof(hash_type) > 4) { prime = long_prime; }
        else { prime = short_prime; }
        
        for(int i=1;i<wsize;i++)
            asize_pot = (asize_pot*asize)% prime; // ugly linear-time power algorithm
        // alloc and clear window
        window = new int[wsize];
        reset();
    }
    
    // init window, hash, and tot_char
    void reset() {
        for(int i=0;i<wsize;i++) window[i]=0;
        // init hash value and related values
        hash=tot_char=0;
    }
    
    hash_type addchar(int c) {
        std::size_t k = tot_char++ % wsize;
        // complex expression to avoid negative numbers
        hash += (prime - (window[k]*asize_pot) % prime); // remove window[k] contribution
        hash = (asize*hash + c) % prime;      //  add char i
        window[k]=c;
        // cerr << get_window() << " ~~ " << window << " --> " << hash << endl;
        return hash;
    }
    
    ~KR_window() {
        delete[] window;
    }
    
};


//------------------------------------------------------------------------------

/*
  Temporary file names have the pattern "prefix_hostname_pid_counter", where
  - prefix is given as an argument to getName();
  - hostname is the name of the host;
  - pid is the process id; and
  - counter is a running counter starting from 0.
  The generated names are stored until the file is deleted with remove(). All
  remaining temporary files are deleted when the program exits (normally or
  with std::exit()).
  TempFile is thread-safe.
*/

namespace TempFile
{
    extern const std::string DEFAULT_TEMP_DIR;
    extern std::string temp_dir;
    
    void setDirectory(const std::string& directory);
    std::string getName(const std::string& name_part);
    void remove(std::string& filename);  // Also clears the filename.
}

void truncate_file(std::string file_name, std::size_t new_size_in_bytes);

bool is_gzipped(std::ifstream& in);

//------------------------------------------------------------------------------

} // end namespace vcfbwt

#endif //utils_hpp
