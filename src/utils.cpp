//
//  utils.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//


#include <utils.hpp>

//------------------------------------------------------------------------------

vcfbwt::hash_type
modular_pow(vcfbwt::hash_type base, vcfbwt::hash_type exponent, vcfbwt::hash_type modulus)
{
    if (modulus == 1) { return 0; }
    if (exponent == 0) { return 1; }
    vcfbwt::hash_type c = 1;
    for (vcfbwt::hash_type e_prime = 0; e_prime <= exponent - 1; e_prime++)
    {
        c = (c * base) % modulus;
    }
    return c;
}


vcfbwt::KarpRabinHash::KarpRabinHash(size_type n) : window_length(n)
{
    this->constant = kr_constant;
    this->prime = kr_prime;
    constant_to_n_minus_one_mod = modular_pow(constant, window_length - 1, prime);
}

void
vcfbwt::KarpRabinHash::reset() { this->hash_value = 0; }

void
vcfbwt::KarpRabinHash::initialize(const std::string& window)
{
    constant_to_n_minus_one_mod = modular_pow(constant, window_length - 1, prime);
    
    assert(window.size() == this->window_length);
    for (hash_type i = 0; i < this->window_length; i++)
    {
        hash_value += window[window.size() - 1 - i] * modular_pow(constant, i, prime);
        hash_value = hash_value % prime;
    }
}

void
vcfbwt::KarpRabinHash::update(char char_out, char char_in)
{
    hash_value = hash_value + prime; // negative avoider
    hash_value = hash_value - ((constant_to_n_minus_one_mod * char_out) % prime);
    hash_value = ((constant * hash_value) % prime) + char_in;
    hash_value = hash_value % prime;
}

vcfbwt::hash_type
vcfbwt::KarpRabinHash::string_hash(const std::string_view& s)
{
    hash_type result = 0;
    hash_type constant_to_n_minus_one_mod = modular_pow(kr_constant, s.size() - 1, kr_prime);
    
    for (hash_type i = 0; i < s.size(); i++)
    {
        result += s[s.size() - 1 - i] * modular_pow(kr_constant, i, kr_prime);
        result = result % kr_prime;
    }
    
    return result;
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
        std::lock_guard<std::mutex>(this->tempfile_lock);
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
    std::lock_guard<std::mutex>(handler.tempfile_lock);
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
        std::lock_guard<std::mutex>(handler.tempfile_lock);
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
            std::lock_guard<std::mutex>(handler.tempfile_lock);
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
    std::lock_guard<std::mutex> lock(writes_counter.write_stats_lock);
    writes_counter.bytes_wrote += num_of_bytes;
}

bool
vcfbwt::ref_smaller(
std::pair<std::reference_wrapper<std::string>, vcfbwt::hash_type> a,
std::pair<std::reference_wrapper<std::string>, vcfbwt::hash_type> b)
{
    return (a.first.get() < b.first.get());
}


//------------------------------------------------------------------------------
