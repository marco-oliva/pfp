//
//  utils.cpp
//


#include <utils.hpp>

//------------------------------------------------------------------------------

vcfbwt::hash_type
modular_pow(vcfbwt::hash_type base, vcfbwt::hash_type exponent, vcfbwt::hash_type modulus)
{
    if (modulus == 1) { return 0; }
    if (exponent == 0) { return 1; }
    unsigned __int128 c = 1;
    for (vcfbwt::hash_type e_prime = 0; e_prime <= exponent - 1; e_prime++)
    {
        c = (c * base) % modulus;
    }
    return c;
}


vcfbwt::KarpRabinHash::KarpRabinHash(std::size_t n, bool debug) : window_length(n), debug_(debug)
{
    this->constant = kr_constant;
    this->prime = kr_prime;
}

void
vcfbwt::KarpRabinHash::reset() { this->hash_value = 0; }

void
vcfbwt::KarpRabinHash::initialize(const char_type* data, std::size_t length)
{
    if (debug_) { debug_content_ = std::string((char*) data, length); }
    
    constant_to_n_minus_one_mod = modular_pow(constant, window_length - 1, prime);
    
    for (hash_type i = 0; i < this->window_length; i++)
    {
        char_type c = data[length - 1 - i];
        hash_value += c * modular_pow(constant, i, prime);
        hash_value = hash_value % prime;
    }
}

void
vcfbwt::KarpRabinHash::update(vcfbwt::char_type char_out, vcfbwt::char_type char_in)
{
    if (debug_) { assert(debug_content_[0] == char_out); debug_content_.erase(0, 1); debug_content_.append(1, char_in); }
    
    hash_value = hash_value + prime; // negative avoider
    hash_value = hash_value - ((constant_to_n_minus_one_mod * char_out) % prime);
    hash_value = ((constant * hash_value) % prime) + char_in;
    hash_value = hash_value % prime;
}

vcfbwt::hash_type
vcfbwt::KarpRabinHash::string_hash(const char_type* data, std::size_t length)
{
    hash_type result = 0;

    for (std::size_t i = 0; i < length; i++)
    {
        char_type c = data[length - 1 - i];
        result += c * modular_pow(kr_constant, i, kr_prime);
        result = result % kr_prime;
    }

    return result;
}

//------------------------------------------------------------------------------

vcfbwt::KarpRabinHash4::KarpRabinHash4(std::size_t n, bool debug) : window_length(n), debug_(debug)
{
    this->window_length_32 = window_length / 4;
    this->constant = kr_constant;
    this->prime = kr_prime;
}

void
vcfbwt::KarpRabinHash4::reset() { this->hash_value = 0; }

void
vcfbwt::KarpRabinHash4::initialize(const char_type* data, std::size_t length)
{
    if (debug_) { debug_content_ = std::string((char*) data, length); }

    constant_to_n_minus_one_mod = modular_pow(constant, window_length_32 - 1, prime);

    uint32_t* string32 = (uint32_t*) data;
    for (std::size_t i = 0; i < this->window_length_32; i++) // window always multiple of 4
    {
        hash_value += string32[i] * modular_pow(constant, i, prime);
        hash_value = hash_value % prime;
    }
}

void
vcfbwt::KarpRabinHash4::update(const vcfbwt::char_type* chars_out, const vcfbwt::char_type* chars_in)
{
    if (debug_)
    {
        for (std::size_t  i = 0; i < 4; i++) { assert(debug_content_[i] == chars_out[i]); }
        debug_content_.erase(0, 4);
        debug_content_.append((char*) chars_in, 4);
    }

    uint32_t *s32c_in = (uint32_t*) chars_in, *s32c_out = (uint32_t*) chars_out;
    hash_value = hash_value + prime; // negative avoider
    hash_value = hash_value - ((constant_to_n_minus_one_mod * (*s32c_out)) % prime);
    hash_value = ((constant * hash_value) % prime) + (*s32c_in);
    hash_value = hash_value % prime;
}

vcfbwt::hash_type
vcfbwt::KarpRabinHash4::string_hash(const char_type* data, std::size_t length)
{
    hash_type result = 0;
    uint32_t* string32 = (uint32_t*) data;
    for (std::size_t i = 0; i < length / 4; i++) // window always multiple of 4
    {
        result += string32[i] * modular_pow(kr_constant, i, kr_prime);
        result = result % kr_prime;
    }

    return result;
}

//------------------------------------------------------------------------------

vcfbwt::Mersenne_KarpRabinHash::Mersenne_KarpRabinHash(std::size_t n, bool debug) : window_length(n), debug_(debug)
{
    this->constant_to_n_minus_one_mod = modular_pow(kr_base, window_length - 1, kr_prime);
}

void
vcfbwt::Mersenne_KarpRabinHash::reset() { this->hash_value = 0; }

void
vcfbwt::Mersenne_KarpRabinHash::initialize(const char_type* data, std::size_t length)
{
    if (debug_) { debug_content_ = std::string((char*) data, length); }

    hash_type h = 0;
    for (std::size_t i = 0; i < this->window_length; i++) // window always multiple of 4
    {
        unsigned __int128 hc = (unsigned __int128)h*kr_base + (char_type) data[i]; //x64 compilers generate standard mul instruction
        hash_type lo = (hash_type)hc, hi = (hash_type)(hc >> 64);
        h = mersenne_modulo(lo, hi, kr_prime, kr_p_pow);
    }
    this->hash_value = h;
}

void
vcfbwt::Mersenne_KarpRabinHash::update(const vcfbwt::char_type char_out, const vcfbwt::char_type char_in)
{
    if (debug_) { assert(debug_content_[0] == char_out); debug_content_.erase(0, 1); debug_content_.append(1, char_in); }

    unsigned __int128 hc = (unsigned __int128) hash_value + kr_prime; // negative avoiders

    unsigned __int128 to_remove = (unsigned __int128) constant_to_n_minus_one_mod * char_out;
    hash_type lo = (hash_type)to_remove, hi = (hash_type)(to_remove >> 64);
    to_remove = mersenne_modulo(lo, hi, kr_prime, kr_p_pow);

    hc = hc - to_remove; // take char_out out
    hc = (hc * kr_base) + char_in;

    lo = (hash_type)hc, hi = (hash_type)(hc >> 64);
    hash_value = mersenne_modulo(lo, hi, kr_prime, kr_p_pow);
}

vcfbwt::hash_type
vcfbwt::Mersenne_KarpRabinHash::string_hash(const char_type* data, std::size_t length)
{
    hash_type h = 0; //not 0!
    for (std::size_t i = 0; i < length; i++) // window always multiple of 4
    {
        unsigned __int128 hc = (unsigned __int128)h*kr_base + (char_type) data[i]; //x64 compilers generate standard mul instruction
        hash_type lo = (hash_type)hc, hi = (hash_type)(hc >> 64);
        h = mersenne_modulo(lo, hi, kr_prime, kr_p_pow);
    }
    return h;
}

//------------------------------------------------------------------------------

vcfbwt::Mersenne_KarpRabinHash4::Mersenne_KarpRabinHash4(std::size_t n, bool debug) : window_length(n), debug_(debug)
{
    assert(this->window_length % 4 == 0 and this->window_length >= 4);
    this->window_length_32 = window_length / 4;

    this->constant_to_n_minus_one_mod = modular_pow(kr_base, window_length_32 - 1, kr_prime);
}

void
vcfbwt::Mersenne_KarpRabinHash4::reset() { this->hash_value = 0; }

void
vcfbwt::Mersenne_KarpRabinHash4::initialize(const char_type* data, std::size_t length)
{
    if (debug_) { debug_content_ = std::string((char*) data, length); }

    uint32_t* s32c_ins = (uint32_t*) data;
    hash_type h = 0;
    for (std::size_t i = 0; i < this->window_length_32; i++) // window always multiple of 4
    {
        unsigned __int128 hc = (unsigned __int128)h*kr_base + s32c_ins[i]; //x64 compilers generate standard mul instruction
        hash_type lo = (hash_type)hc, hi = (hash_type)(hc >> 64);
        h = mersenne_modulo(lo, hi, kr_prime, kr_p_pow);
    }
    this->hash_value = h;
}

void
vcfbwt::Mersenne_KarpRabinHash4::update(const vcfbwt::char_type* chars_out, const vcfbwt::char_type* chars_in)
{
    if (debug_)
    {
        for (std::size_t  i = 0; i < 4; i++) { assert(debug_content_[i] == chars_out[i]); }
        debug_content_.erase(0, 4);
        debug_content_.append((char*) chars_in, 4); // this goes byte by byte
    }

    uint32_t* s32c_in  = (uint32_t*)  chars_in;
    uint32_t* s32c_out = (uint32_t*) chars_out;

    unsigned __int128 hc = (unsigned __int128) hash_value + kr_prime; // negative avoiders

    unsigned __int128 to_remove = (unsigned __int128) constant_to_n_minus_one_mod * (*s32c_out);
    hash_type lo = (hash_type)to_remove, hi = (hash_type)(to_remove >> 64);
    to_remove = mersenne_modulo(lo, hi, kr_prime, kr_p_pow);

    hc = hc - to_remove; // take char_out out
    hc = (hc * kr_base) + (*s32c_in);

    lo = (hash_type)hc, hi = (hash_type)(hc >> 64);
    hash_value = mersenne_modulo(lo, hi, kr_prime, kr_p_pow);
}

vcfbwt::hash_type
vcfbwt::Mersenne_KarpRabinHash4::string_hash(const char_type* data, std::size_t length)
{
    uint32_t* string32 = (uint32_t*) data;
    hash_type h = 0;
    for (std::size_t i = 0; i < length / 4; i++) // window always multiple of 4
    {
        unsigned __int128 hc = (unsigned __int128)h*kr_base + string32[i]; //x64 compilers generate standard mul instruction
        hash_type lo = (hash_type)hc, hi = (hash_type)(hc >> 64);
        h = mersenne_modulo(lo, hi, kr_prime, kr_p_pow);
    }
    return h;
}

//------------------------------------------------------------------------------

const std::string vcfbwt::TempFile::DEFAULT_TEMP_DIR = ".";
std::string vcfbwt::TempFile::temp_dir = vcfbwt::TempFile::DEFAULT_TEMP_DIR;

// By storing the filenames in a static object, we can delete the remaining
// temporary files when std::exit() is called.
struct Handler
{
    std::mutex tempfile_lock;
    std::size_t counter;
    std::set<std::string> filenames;
    
    Handler() : counter(0) {}
    
    ~Handler()
    {
        std::lock_guard<std::mutex> lock_guard(this->tempfile_lock);
        for(auto& filename : this->filenames)
        {
            std::remove(filename.c_str());
        }
        this->filenames.clear();
    }
} handler;

void
vcfbwt::TempFile::setDirectory(const std::string& directory)
{
    std::lock_guard<std::mutex> lock_guard(handler.tempfile_lock);
    if(directory.empty()) { temp_dir = DEFAULT_TEMP_DIR; }
    else if(directory[directory.length() - 1] != '/') { temp_dir = directory; }
    else { temp_dir = directory.substr(0, directory.length() - 1); }
}

std::string
vcfbwt::TempFile::getName(const std::string& name_part)
{
    char hostname[32];
    gethostname(hostname, 32); hostname[31] = 0;
    
    std::string filename;
    {
        std::lock_guard<std::mutex> lock_guard(handler.tempfile_lock);
        filename = temp_dir + '/' + name_part + '_'
        + std::string(hostname) + '_'
        + std::to_string(pid()) + '_'
        + std::to_string(handler.counter);
        handler.filenames.insert(filename);
        handler.counter++;
    }
    
    return filename;
}

void
vcfbwt::TempFile::remove(std::string& filename)
{
    if(!(filename.empty()))
    {
        std::remove(filename.c_str());
        {
            std::lock_guard<std::mutex> lock_guard(handler.tempfile_lock);
            handler.filenames.erase(filename);
        }
        filename.clear();
    }
}

//------------------------------------------------------------------------------

#include <unistd.h>
void
vcfbwt::truncate_file(std::string file_name, std::size_t new_size_in_bytes)
{
    int res = truncate(file_name.c_str(), new_size_in_bytes);
    if (res != 0) { spdlog::error("Error while truncating {} to {} bytes", file_name, new_size_in_bytes); }
}

bool
vcfbwt::is_gzipped(std::ifstream& in)
{
    in.seekg(0);
    char f, s; in.get(f); in.get(s);
    in.seekg(0);
    
    return ((f == '\x1F') and (s == '\x8B'));
}

//------------------------------------------------------------------------------

struct WritesCounter
{
    std::mutex write_stats_lock;
    std::size_t bytes_wrote;
    
    WritesCounter() : bytes_wrote(0) {};
    ~WritesCounter()
    {
        // Can't use spdlog here, could be destroyed after spdlog's objects
        std::cout << "[Disk Write (bytes): " << std::to_string(bytes_wrote) << "]" << std::endl;
    }
    
} writes_counter;

void vcfbwt::DiskWrites::update(std::size_t num_of_bytes)
{
    std::lock_guard<std::mutex> lock_guard(writes_counter.write_stats_lock);
    writes_counter.bytes_wrote += num_of_bytes;
}

//------------------------------------------------------------------------------
