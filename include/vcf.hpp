//
//  vcf.hpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#ifndef vcf_hpp
#define vcf_hpp

#include <string_view>
#include <vector>
#include <unordered_map>
#include <functional>

#include <spdlog/spdlog.h>

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <stdlib.h>

#include <zstr/zstr.hpp>

#include <hts.h>
#include <vcf.h>
#include <tbx.h>
#include <vcfutils.h>

#include <utils.hpp>


namespace vcfbwt
{

//------------------------------------------------------------------------------

struct Variation
{
    std::size_t pos = 0; // already adjusted when read from the vcf by htslib. No need of -1
    std::size_t ref_len = 0;
    std::string alt;
    
    double freq = 0.0;
    bool used = false;
};

class Sample
{
private:
    
    const std::string& reference_;
    const std::vector<Variation>& variations_list;
    
    std::string sample_id;
    
    friend class iterator;
    
public:
    
    std::vector<std::size_t> variations;

    const std::string& id() const { return this->sample_id; }
    
    Sample(const std::string& id, const std::string& ref, const std::vector<Variation>& variations)
    : sample_id(id), reference_(ref), variations_list(variations) {}
    
    const Variation& get_variation(size_type i) const { return this->variations_list[variations[i]]; }
    
    const std::string& get_reference() const { return this->reference_; }
    
    class iterator
    {
    private:
        
        const Sample& sample_;
        std::size_t ref_it_;
        std::size_t sam_it_;
        std::size_t var_it_;
        std::size_t curr_var_it_;
        std::size_t sample_length_;
        const char* curr_char_;
        
    public:
        
        iterator(const Sample& sample);
        
        bool end();
        void operator++();
        void go_to(std::size_t i);
        const char operator*();
        
        bool in_a_variation();
        
        std::size_t get_var_it() const { return var_it_; }
        std::size_t get_sam_it() const { return sam_it_; }
        std::size_t get_ref_it() const { return ref_it_ - 1; } // -1 because the iterator is pointing the next one
        std::size_t next_variation() const;
        std::size_t next_variation_distance() const;
        std::size_t prev_variation() const;
        
        std::size_t length() const { return sample_length_; }
    };
};

//------------------------------------------------------------------------------

class VCF
{

private:
    
    std::string reference;
    
    std::vector<Variation> variations;
    std::vector<Sample> samples;
    std::unordered_map<std::string, std::size_t> samples_id;
    
    std::size_t max_samples = 0;
    
    std::vector<std::size_t> ref_sum_lengths;
    
    void init_vcf(const std::string& vcf_path, std::vector<Variation>& l_variations,
                  std::vector<Sample>& l_samples, std::unordered_map<std::string, std::size_t>& l_samples_id,
                  std::size_t i = 0);
    void init_vcf(const std::string& vcf_path, std::size_t i =0);
    void init_ref(const std::string& ref_path, bool last = true);
    
    void init_multi_vcf(const std::vector<std::string>& vcfs_path);
    void init_multi_ref(const std::vector<std::string>& refs_path);
    

public:
    
    static const std::string vcf_freq;
    
    VCF(const std::string& ref_path, const std::string& vcf_path, std::size_t ms = 0) : max_samples(ms) { init_ref(ref_path); init_vcf(vcf_path); }
    VCF(const std::vector<std::string>& refs_path, const std::vector<std::string>& vcfs_path, std::size_t ms = 0) : max_samples(ms) { init_multi_ref(refs_path); init_multi_vcf(vcfs_path); }
    ~VCF() = default;
    
    std::size_t size() const { if (max_samples != 0) { return max_samples; } else { return samples.size(); } }
    Sample& operator[](std::size_t i) { assert(i < size()); return samples.at(i); }
    const std::vector<Variation>& get_variations() const { return this->variations; }
    const std::string& get_reference() const { return this->reference; }
    void set_max_samples(std::size_t max) { this->max_samples = max; }
    
};

//------------------------------------------------------------------------------

} // end namespace vcfbwt


#endif //vcf_hpp
