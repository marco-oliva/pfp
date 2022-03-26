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
    std::vector<std::string> alt;
    std::vector<int> types; // Type of variation
    // std::vector<int> len; // Number of bases affected. Negative for deletions.
    
    double freq = 0.0;
    bool used = false;
    
    enum variation_type { H1, H2, H12 };
};


class Contig
{
private:
    
    const std::string& reference_;
    const std::vector<Variation>& variations_list;
    
    std::size_t offset_ = 0;
    std::size_t ref_index_ = 0;
    
    std::string contig_id;
    bool is_last_contig = false;
    int last_variation_type = 0;
    int max_ploidy = 0; // Number of haplotypes in the contig
    
    friend class iterator;
    
public:
    
    void set_last(int type) { this->is_last_contig = true; last_variation_type = type; }
    bool last(const int type) const { return this->is_last_contig and (type == last_variation_type); }
    size_t offset() const { return this->offset_; }
    size_t get_reference_index() const { return this->ref_index_; }
    void set_ploidy(int ploidy) { this->max_ploidy = ploidy; }
    void update_ploidy(int ploidy) { this->max_ploidy = std::max(this->max_ploidy, ploidy); }
    int get_ploidy() const { return this->max_ploidy; }
    
    // variation, variation type
    std::vector<std::size_t> variations;
    std::vector<std::vector<int>> genotypes;

    const std::string& id() const { return this->contig_id; }
    
    Contig(const std::string& id, const std::string& ref, const std::vector<Variation>& variations, const size_t offset, const size_t ref_idx)
    : contig_id(id), reference_(ref), variations_list(variations), offset_(offset), ref_index_(ref_idx) {}
    
    Contig(Contig&& other)
    :   reference_(other.reference_), 
        variations_list(other.variations_list),
        contig_id(std::move(other.contig_id)),
        offset_(other.offset_),
        ref_index_(other.ref_index_),
        is_last_contig(other.is_last_contig),
        last_variation_type(other.last_variation_type),
        max_ploidy(other.max_ploidy),
        variations(std::move(other.variations)),
        genotypes(std::move(other.genotypes))
    {

    }
    
    Contig(const Contig& other)
    :   reference_(other.reference_), 
        variations_list(other.variations_list),
        contig_id(other.contig_id),
        offset_(other.offset_),
        ref_index_(other.ref_index_),
        is_last_contig(other.is_last_contig),
        last_variation_type(other.last_variation_type),
        max_ploidy(other.max_ploidy),
        variations(other.variations),
        genotypes(other.genotypes)
    {

    }

    const Variation& get_variation(size_type i) const { return this->variations_list[variations[i]]; }
    
    const std::string& get_reference() const { return this->reference_; }
    
    class iterator
    {
    private:
        
        const Contig& contig_;
        std::size_t genotype;
        std::size_t ref_it_;
        std::size_t sam_it_;
        std::size_t var_it_;
        std::size_t prev_variation_it;
        std::size_t curr_var_it_;
        std::size_t contig_length_;
        const char* curr_char_;
        
    public:
        
        iterator(const Contig& contig, std::size_t genotype = 0);
        
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
        
        std::size_t length() const { return contig_length_; }
    };
};


class Sample
{  
protected:

    std::string sample_id;

    bool is_last_sample = false;
    int last_variation_type = 0;

    friend class iterator;

public:

    Sample(const std::string& id)
    : sample_id(id){}

    bool set_last(int type) 
    { 
        this->is_last_sample = true;
        this->last_variation_type = type;
        size_t i = 0;
        const size_t n = contigs.size();

        while (i < n and this->contigs[n - i - 1].get_ploidy() <= type){++i;}

        if(i < n) this->contigs[n - i - 1].set_last(type);
        else return false;

        return true;
    } //TODO: Check if we need to set also the contig flag
    bool last(const int type) const { return this->is_last_sample and (type == last_variation_type); }

    const std::string& id() const { return this->sample_id; }
    
    // The actual list of contigs of the sample.
    std::vector<Contig> contigs;

    class iterator
    {
    private:
        
        const Sample& sample_;
        std::size_t genotype;
        vcfbwt::Contig::iterator* c_it_ = nullptr;
        std::vector<Contig>::const_iterator v_it_;

    public:
        
        iterator(const Sample& sample, std::size_t genotype = 0);
        ~iterator() {if(c_it_ != nullptr) delete c_it_;}
        bool end();
        void operator++();
        const char operator*();
    };
};

//------------------------------------------------------------------------------

class VCF
{

private:
    
    std::string reference; // TODO: Replace it with as many parses as references.
    std::vector<std::string> references;
    std::vector<std::string> references_name;
    std::unordered_map<std::string, std::size_t> references_id;
    
    std::vector<std::vector<Variation>> variations;
    std::vector<Contig> contigs;
    std::vector<Sample> samples;
    std::unordered_map<std::string, std::size_t> samples_id;
    
    std::size_t max_samples = 0;
    std::set<std::string> input_samples;
    std::vector<std::size_t> populated_samples;

    std::vector<std::size_t> ref_sum_lengths;

    void init_samples(const std::string& samples_path);
    
    void init_vcf(const std::string& vcf_path,
                  std::vector<Sample>& l_samples, std::unordered_map<std::string, 
                  std::size_t>& l_samples_id);
    void init_vcf(const std::string& vcf_path, std::size_t i = 0);
    void init_ref(const std::string &ref_path, const size_t w, bool last = true);

    void init_multi_vcf(const std::vector<std::string>& vcfs_path);
    void init_multi_ref(const std::vector<std::string> &refs_path, const size_t w);

    void init_contigs();

public:
    
    static const std::string vcf_freq;

    VCF(const std::string &ref_path, const std::string &vcf_path, const std::string &samples_path, const size_t w, const int last_genotype = 0, std::size_t ms = 0) : max_samples(ms)
    {
        if (samples_path != "") { init_samples(samples_path); }
        init_ref(ref_path, w); init_contigs(); init_vcf(vcf_path);
        for (std::size_t i = 0; i < samples.size(); i++)
        { if (not samples.at(i).contigs.empty()) { this->populated_samples.push_back(i); } }

        this->samples.at(populated_samples.back()).set_last(last_genotype);
    }

    VCF(const std::vector<std::string> &refs_path, const std::vector<std::string> &vcfs_path, const std::string &samples_path, const size_t w, const int last_genotype = 0, std::size_t ms = 0) : max_samples(ms)
    {
        if (samples_path != "") { init_samples(samples_path); }
        init_multi_ref(refs_path, w); init_contigs(); init_multi_vcf(vcfs_path);
        for (std::size_t i = 0; i < samples.size(); i++)
        { if (not samples.at(i).contigs.empty()) { this->populated_samples.push_back(i); } }

        this->samples.at(populated_samples.back()).set_last(last_genotype);
    }
    
    ~VCF() = default;
    
    std::size_t size() const { return this->populated_samples.size(); }
    Sample& operator[](std::size_t i) { assert(i < size()); return samples.at(populated_samples.at(i)); }
    const std::vector<std::vector<Variation>>& get_variations() const { return this->variations; }
    const std::string& get_reference() const { return this->reference; }
    const std::vector<std::string>& get_references() const { return this->references; }
    const std::vector<std::string>& get_references_name() const { return this->references_name; }
    void set_max_samples(std::size_t max) { this->max_samples = max; }
};

//------------------------------------------------------------------------------

} // end namespace vcfbwt


#endif //vcf_hpp
