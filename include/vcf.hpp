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

#include <zstr.hpp>

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
    std::vector<std::string> alt;
    
    double freq = 0.0;
    bool used = false;
    
    enum variation_type { H1, H2, H12 };
};

class Sample
{
private:
    
    const std::string& reference_;
    const std::vector<Variation>& variations_list;
    
    std::string sample_id;
    bool is_last_sample = false;
    int last_variation_type = 0;

    friend class iterator;
    
public:
    void set_last(const int type){  this->is_last_sample = true; this->last_variation_type = type; }
    bool last(const int type) const { return this->is_last_sample and ( type == this->last_variation_type ); }

    // variation, variation type
    std::vector<std::size_t> variations;
    std::vector<std::vector<int>> genotypes;

    const std::string& id() const { return this->sample_id; }
    
    Sample(const std::string& id, const std::string& ref, const std::vector<Variation>& variations)
    : sample_id(id), reference_(ref), variations_list(variations) {}
    
    const Variation& get_variation(std::size_t i) const { return this->variations_list[variations[i]]; }
    
    const std::string& get_reference() const { return this->reference_; }
    
    class iterator
    {
    private:
        
        const Sample& sample_;
        std::size_t genotype;
        std::size_t ref_it_;
        std::size_t sam_it_;
        std::size_t var_it_;
        std::size_t prev_variation_it;
        std::size_t curr_var_it_;
        std::size_t sample_length_;
        const char* curr_char_;
        
    public:
        
        iterator(const Sample& sample, std::size_t genotype = 0);
        
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
        std::size_t prev_variation_end() const;
        
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
    std::set<std::string> input_samples;
    std::vector<std::size_t> populated_samples;

    std::vector<std::size_t> ref_sum_lengths;

    void init_samples(const std::string& samples_path);
    
    void init_vcf(const std::string& vcf_path, std::vector<Variation>& l_variations,
                  std::vector<Sample>& l_samples, std::unordered_map<std::string, std::size_t>& l_samples_id,
                  std::size_t i = 0);
    void init_vcf(const std::string& vcf_path, std::size_t i = 0);
    void init_ref(const std::string& ref_path, bool last = true);
    
    void init_multi_vcf(const std::vector<std::string>& vcfs_path);
    void init_multi_ref(const std::vector<std::string>& refs_path);


public:
    
    static const std::string vcf_freq;

    VCF(const std::string &ref_path, const std::string &vcf_path, const std::string &samples_path, std::size_t ms = 0, const int last_genotype = 0) : max_samples(ms)
    {
        if (samples_path != "") { init_samples(samples_path); }
        init_ref(ref_path); init_vcf(vcf_path);
        for (std::size_t i = 0; i < samples.size(); i++)
        { if (not samples.at(i).variations.empty()) { this->populated_samples.push_back(i); } }

        this->samples.at(populated_samples.back()).set_last(last_genotype);
    }

    VCF(const std::vector<std::string> &refs_path, const std::vector<std::string> &vcfs_path, const std::string &samples_path, std::size_t ms = 0, const int last_genotype = 0) : max_samples(ms)
    {
        if (samples_path != "") { init_samples(samples_path); }
        init_multi_ref(refs_path); init_multi_vcf(vcfs_path);
        for (std::size_t i = 0; i < samples.size(); i++)
        { if (not samples.at(i).variations.empty()) { this->populated_samples.push_back(i); } }

        this->samples.at(populated_samples.back()).set_last(last_genotype);
    }
    
    ~VCF() = default;
    
    std::size_t size() const { return this->populated_samples.size(); }
    Sample& operator[](std::size_t i) { assert(i < size()); return samples.at(populated_samples.at(i)); }
    const std::vector<Variation>& get_variations() const { return this->variations; }
    const std::string& get_reference() const { return this->reference; }
    void set_max_samples(std::size_t max) { this->max_samples = max; }
};

//------------------------------------------------------------------------------

} // end namespace vcfbwt


#endif //vcf_hpp
