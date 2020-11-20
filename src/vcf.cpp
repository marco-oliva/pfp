//
//  vcf.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#include <vcf.hpp>

const std::string vcfbwt::VCF::vcf_freq = "AF";


//void
//vcfbwt::Sample::init_sum()
//{
//    if (variations.size() == 0) { return; }
//
//    variations_length_sum.reserve(variations.size());
//    variations_length_sum.push_back(variations[0]);
//    for (std::size_t i = 1; i < variations.size(); i++)
//    {
//        variations_length_sum.push_back(variations[i] + variations[i - 1]);
//    }
//}

//------------------------------------------------------------------------------

vcfbwt::Sample::iterator::iterator(const Sample& sample) :
sample_(sample), ref_it_(0), sam_it_(0), var_it_(0), curr_var_it_(0), curr_char_(NULL), sample_length_(sample.reference_.size())
{
    // Compute sample length, might take some time
    long long int indels = 0; // could be negative, so int
    for (auto var_id : sample.variations)
    {
        indels = indels + (sample_.variations_list[var_id].alt.size() - sample_.variations_list[var_id].ref_len);
    }
    sample_length_ = sample_length_ + indels;
    this->operator++();
}

bool
vcfbwt::Sample::iterator::end() { return ref_it_ > this->sample_.reference_.size(); }

const char
vcfbwt::Sample::iterator::operator*() { return *curr_char_; }

bool
vcfbwt::Sample::iterator::in_a_variation()
{
    const Variation& curr_variation = sample_.variations_list[sample_.variations[var_it_]];
    return (ref_it_ == curr_variation.pos);
}

std::size_t
vcfbwt::Sample::iterator::next_variation() const
{
    if (var_it_ < sample_.variations.size()) { return sample_.get_variation(var_it_).pos; }
    else { return sample_.reference_.size() - 1; }
}

std::size_t
vcfbwt::Sample::iterator::prev_variation() const
{
    if (var_it_ == 0) { spdlog::error("vcfbwt::Sample::iterator::prev_variation() var_it == 0"); std::exit(EXIT_FAILURE); }
    return sample_.get_variation(var_it_ - 1).pos;
}

std::size_t
vcfbwt::Sample::iterator::next_variation_distance() const
{
    return this->next_variation() - (this->ref_it_ - 1);
}

void
vcfbwt::Sample::iterator::operator++()
{
    // Ci sono ancora variazioni da processare
    if (var_it_ < sample_.variations.size())
    {
        const Variation& curr_variation = sample_.variations_list[sample_.variations[var_it_]];
        
        if (ref_it_ < curr_variation.pos)
        {
            curr_char_ = &(sample_.reference_[ref_it_]); ref_it_++; sam_it_++;
            return;
        }
        
        // Più nucleotidi nella variaizione
        if (curr_variation.alt.size() > 1)
        {
            if (curr_var_it_ < curr_variation.alt.size() - 1)
            {
                curr_char_ = & curr_variation.alt[curr_var_it_];
                curr_var_it_++; sam_it_++;
            }
            else
            {
                curr_char_ = & curr_variation.alt.back();
                var_it_++; curr_var_it_ = 0; ref_it_ += curr_variation.ref_len; sam_it_++;
            }
            
        }
        else
        {
            curr_char_ = & curr_variation.alt.front();
            ref_it_ += curr_variation.ref_len; sam_it_++; var_it_++;
        }
    }
    // Non ci sono più variazioni, itera sulla reference
    else { curr_char_ = &(sample_.reference_[ref_it_]); ref_it_++; sam_it_++; }
}

void
vcfbwt::Sample::iterator::go_to(std::size_t i)
{
    if (i >= this->sample_.reference_.size())
    {
        spdlog::error("vcfbwt::Sample::iterator::go_to(std::size_t i) Error: i >= reference.size()");
        std::exit(EXIT_FAILURE);
    }
    if (i < ref_it_)
    {
        spdlog::error("vcfbwt::Sample::iterator::go_to(std::size_t i) Error: i < ref_it_ (going back not allowed)");
        std::exit(EXIT_FAILURE);
    }
    
    // Not an efficient implementation, TODO: make it more efficient
    while (ref_it_ < i)
    {
        this->operator++();
    }
}

//------------------------------------------------------------------------------

void
vcfbwt::VCF::init_ref(const std::string& ref_path)
{
    spdlog::info("Reading reference file: {}", ref_path);
    std::ifstream in_stream(ref_path);
    if (not is_gzipped(in_stream))
    {
        spdlog::error("Reference file expected gzip compressed");
        std::exit(EXIT_FAILURE);
    }
    
    zstr::istream is(in_stream);
    std::string line;
    
    while (getline(is, line)) { if ( not (line.empty() or line[0] == '>') ) { reference.append(line); } }
    
    ref_sum_lengths.push_back(reference.size());
    
    spdlog::info("Done reading {}", ref_path);
}

//------------------------------------------------------------------------------

void
vcfbwt::VCF::init_vcf(const std::string& vcf_path, std::size_t i)
{
    // open VCF file
    htsFile * inf = bcf_open(vcf_path.c_str(), "r");
    if (inf == NULL)
    {
        spdlog::error("Can't open vcf file: {}", vcf_path);
        std::exit(EXIT_FAILURE);
    }
    
    spdlog::info("Parsing vcf: {}", vcf_path);
    
    // read header
    bcf_hdr_t *hdr = bcf_hdr_read(inf);
    
    // get samples ids from header
    std::size_t n_samples = bcf_hdr_nsamples(hdr);
    this->samples.reserve(n_samples);
    
    std::size_t size_before = samples.size();
    for (std::size_t i = 0; i < n_samples; i++)
    {
        vcfbwt::Sample s(std::string(hdr->samples[i]), this->reference, this->variations);
        samples.push_back(s);
        
        samples_id.insert(std::make_pair(s.id(), i));
    }
    spdlog::info("{} samples in the vcf", samples.size() - size_before);
    
    // struct for storing each record
    bcf1_t *rec = bcf_init();
    if (rec == NULL)
    {
        spdlog::error("Error while parsing vcf file: {}", vcf_path);
        bcf_close(inf);
        bcf_hdr_destroy(hdr);
        std::exit(EXIT_FAILURE);
    }
    
    // start parsing
    while (bcf_read(inf, hdr, rec) == 0)
    {
        vcfbwt::Variation var;
    
        int allocation_size = 0;
        float* freq = NULL;
        bcf_get_info_float(hdr, rec, vcf_freq.c_str(), &freq, &allocation_size);
    
        var.ref_len = rec->rlen;
        std::size_t offset = i != 0 ? ref_sum_lengths[i-1] : 0; // when using multiple vcfs
        var.pos = rec->pos + offset;
        var.freq = *freq;
    
        bcf_unpack(rec, BCF_UN_ALL);
        int type = bcf_get_variant_types(rec);
    
        this->variations.push_back(var);
        
        int32_t *gt_arr = NULL, ngt_arr = 0;
        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if ( ngt > 0 )
        {
            int max_ploidy = ngt/n_samples;
            bool skip_this_variation = false;
            for (std::size_t i = 0; i < n_samples; i++)
            {
                if (skip_this_variation) { break; }
                int32_t *ptr = gt_arr + i*max_ploidy;
                for (std::size_t j = 0; j < max_ploidy; j++)
                {
                    // if true, the sample has smaller ploidy
                    if ( ptr[j]==bcf_int32_vector_end ) break;

                    // missing allele
                    if ( bcf_gt_is_missing(ptr[j]) ) continue;
                    
                    if (bcf_gt_allele(ptr[j]))
                    {
                        // the VCF 0-based allele index
                        int allele_index = bcf_gt_allele(ptr[j]);
                        
                        if (variations.back().alt.size() == 0)
                            variations.back().alt = rec->d.allele[allele_index];
    
                        // Skip symbolic allele
                        if (variations.back().alt[0] == '<')
                        {
                            //spdlog::warn("vcfbwt::VCF::init_vcf: Skipping symbolic allele at pos {}", variations.back().pos);
                            skip_this_variation = true;
                            continue;
                        }
                        
                        auto id = samples_id.find(std::string(hdr->samples[i]));
                        if (id != samples_id.end())
                        {
                            // Adding this variation to a sample only if:
                            // - has not been inserted right before this insertion
                            if ( not (
                            ((this->samples[id->second].variations.size() > 0) and
                            (this->samples[id->second].variations.back() == this->variations.size() - 1))
                            ))
                            {
                                this->samples[id->second].variations.push_back(this->variations.size() - 1);
                            }
                        }
                    }
                }
            }
        }
        free(gt_arr);
        free(freq);
    }
    
    // free allocated memory
    bcf_hdr_destroy(hdr);
    bcf_close(inf);
    bcf_destroy(rec);
    
    spdlog::info("Parsing vcf done");
    
    // print some statistics
    spdlog::info("Variations size [{}]: {}GB", variations.size(), inGigabytes(variations.size() * sizeof(Variation)));
    spdlog::info("Reference size: {} GB", inGigabytes(reference.size()));
    
    std::size_t tot_a_s = 0;
    for (auto& s : this->samples)
    {
        tot_a_s += s.variations.size();
    }
    spdlog::info("Samples size: {} GB", inGigabytes(tot_a_s * 8));
}

//------------------------------------------------------------------------------

void
vcfbwt::VCF::init_multi_ref(const std::vector<std::string>& refs_path)
{
    if (refs_path.size() == 0) { spdlog::error("No reference file provided"); std::exit(EXIT_FAILURE); }
    
    spdlog::info("Opening {} ref files, assuming input oreder refelcts the intended genome order", refs_path.size());
    for (auto& path : refs_path) { init_ref(path); }
}

//------------------------------------------------------------------------------

void
vcfbwt::VCF::init_multi_vcf(const std::vector<std::string>& vcfs_path)
{
    if (vcfs_path.size() == 0) { spdlog::error("No vcf file provided"); std::exit(EXIT_FAILURE); }
    
    spdlog::info("Opening {} vcf files, assuming input oreder refelcts the intended genome order", vcfs_path.size());
    for (std::size_t i = 0; i < vcfs_path.size(); i++) { init_vcf(vcfs_path[i], i); }
}

//------------------------------------------------------------------------------


