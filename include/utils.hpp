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
#include <list>
#include <unordered_map>

#include <unistd.h>
#include <math.h>

#include <spdlog/spdlog.h>
#include <mio/mmap.hpp>
#include <omp.h>

#include <malloc_count/malloc_count.h>
#include <murmur/MurmurHash3.h>

namespace vcfbwt
{

typedef std::uint64_t long_type;
typedef std::uint32_t short_type;
static long_type  long_prime   = 27162335252586509;
static short_type short_prime  =        1999999973;

#ifdef PFP_LONG_TYPE
//pragma message("Using 64 bit parse")
typedef long_type  size_type;
#else
//#pragma message("Using 32 bit parse")
typedef short_type size_type;
#endif
typedef long_type  hash_type;

constexpr std::size_t KILOBYTE      = 1024;
constexpr std::size_t MEGABYTE      = KILOBYTE * KILOBYTE;
constexpr std::size_t GIGABYTE      = KILOBYTE * MEGABYTE;

constexpr double KILOBYTE_DOUBLE    = 1024.0;
constexpr double MILLION_DOUBLE     = 1000000.0;
constexpr double MEGABYTE_DOUBLE    = KILOBYTE_DOUBLE * KILOBYTE_DOUBLE;
constexpr double GIGABYTE_DOUBLE    = KILOBYTE_DOUBLE * MEGABYTE_DOUBLE;

constexpr std::size_t MILLION       = 1000000;
constexpr std::size_t BILLION       = 1000 * MILLION;

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

inline hash_type
string_hash(const char* s, std::size_t size)
{
    long_type hash[2] = {0};
    if (size >= std::numeric_limits<int>::max()) { spdlog::error("String too big"); exit(EXIT_FAILURE); }
    MurmurHash3_x64_128(s, size, short_prime, hash);
    return hash[0];
}

class KarpRabinHash
{
public:
    KarpRabinHash(size_type n);
    
    void initialize(const std::string& window);
    void update(char char_out, char char_in);
    void reset();
    
    const hash_type& get_hash() const { return this->hash_value; }
    void set_constant(const hash_type& c) { this->constant = c; }
    void set_prime(const hash_type& p) { this->prime = p; }
    
    
private:
    hash_type constant;
    hash_type prime;
    hash_type constant_to_n_minus_one_mod;
    size_type window_length;
    hash_type hash_value = 0;
    
public:
    constexpr static hash_type kr_prime    = 1999999973;
    constexpr static hash_type kr_constant = 256;
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

// Disk space statistics
namespace DiskWrites
{
    void update(std::size_t num_of_bytes);
};

//------------------------------------------------------------------------------


} // end namespace vcfbwt

#endif //utils_hpp
