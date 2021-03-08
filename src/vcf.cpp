//
//  vcf.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#include <vcf.hpp>
#include <pfp_algo.hpp>

const std::string vcfbwt::VCF::vcf_freq = "AF";


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
vcfbwt::VCF::init_ref(const std::string& ref_path, bool last)
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
    if (not last) { reference.push_back(pfp::SPECIAL_TYPES::DOLLAR_PRIME); }

    ref_sum_lengths.push_back(reference.size());

    spdlog::info("Done reading {}", ref_path);
}

//------------------------------------------------------------------------------

std::vector<int> decompose_genotype_fast(const string& genotype) {
    // Rather than doing lots of splits, we just squash atoi in with a single
    // pass over the characters in place.
    
    // We'll fill this in
    std::vector<int> to_return;
    
    // Guess at size
    to_return.reserve(2);
    
    // We use this for our itoa
    int number = 0;
    for (size_t i = 0; i < genotype.size(); i++) {
        switch(genotype[i]) {
        case '.':
            // We have a missing allele.
            number = vcflib::NULL_ALLELE;
            break;
        case '|':
        case '/':
            // We've terminated a field
            to_return.push_back(number);
            number = 0;
            break;
        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
            // We have a digit, so add it to the growing number
            number *= 10;
            number += (genotype[i] - '0');
            break;
        default:
            throw std::runtime_error("Invalid genotype character in " + genotype);
            break;
        }
    }
    if(!genotype.empty()) {
        // Finish the last field
        to_return.push_back(number);
    }
    return to_return;
}

void
vcfbwt::VCF::init_vcf(const std::string& vcf_path, std::vector<Variation>& l_variations,
                      std::vector<Sample>& l_samples, std::unordered_map<std::string, std::size_t>& l_samples_id,
                      std::size_t i)
{
    // offset
    std::size_t offset = i != 0 ? ref_sum_lengths[i-1] : 0; // when using multiple vcfs
    
    // open VCF file
    vcflib::VariantCallFile vcf_file;
    std::string vcf_path_non_cost = vcf_path;
    vcf_file.open(vcf_path_non_cost);
    
    if (not vcf_file.is_open())
    {
        spdlog::error("Can't open vcf file: {}", vcf_path);
        std::exit(EXIT_FAILURE);
    }
    
    spdlog::info("Parsing vcf: {}", vcf_path);
    
    std::size_t size_before = l_samples.size();
    for (auto& sample_name : vcf_file.sampleNames)
    {
        vcfbwt::Sample s(sample_name, this->reference, l_variations);
        if (l_samples_id.find(s.id()) == l_samples_id.end())
        {
            l_samples.push_back(s);
            l_samples_id.insert(std::make_pair(s.id(), l_samples.size() - 1));
        }
    }
    
    spdlog::debug("{} new l_samples in the vcf, tot: {}", l_samples.size() - size_before, l_samples.size());
    
    // start parsing
    vcflib::Variant variant(vcf_file);
    while (vcf_file.getNextVariant(variant))
    {
        // Skip symbolic allele
        if (variant.isSymbolicSV()) { continue; }
        
        vcfbwt::Variation var;
        var.ref_len = variant.ref.size();
        var.pos = variant.position + offset - 1;
        var.freq = 0;
        var.alt = variant.alt[0];
    
        for (auto& kv : variant.samples)
        {
            // Go through all the parsed FORMAT fields for each sample
            auto& sample_name = kv.first;
            
            if (l_samples_id.find(sample_name) == l_samples_id.end()) { spdlog::error("ERROR"); std::exit(EXIT_FAILURE); }
            std::size_t sample_idx = l_samples_id.at(sample_name);
            
            if (sample_idx > this->max_samples) { continue; }
            
            // Pull out the GT value. Explode if there isn't one (though
            // something like "." is acceptable)
            auto& gt_string = kv.second.at("GT").at(0);
            
            // Decompose it and fill in the genotype slot for this sample.
            std::vector<int> out = decompose_genotype_fast(gt_string);
            
            if (out[0] > 0)
            {
                // Add vartiation only if in one of the desired samples
                if (l_variations.size() == 0 or l_variations.back().pos != var.pos)
                {
                    l_variations.push_back(var);
                }
                
                l_variations.back().freq += 1;
                l_variations.back().used = true;
                // Add variation to sample
                l_samples[sample_idx].variations.push_back(l_variations.size() - 1);
            }
        }
    }
    
    // Compute normalized variations frequency
    std::size_t number_of_samples = 0;
    for (auto& s : l_samples) { if (s.variations.size() > 0) { number_of_samples += 1; } }
    for (auto& v : l_variations) { v.freq = v.freq / double(number_of_samples); }
    
    // print some statistics
    spdlog::info("Variations size [{}]: {}GB", l_variations.size(), inGigabytes(l_variations.size() * sizeof(Variation)));
    spdlog::info("Reference size: {} GB", inGigabytes(reference.size()));
    
    std::size_t tot_a_s = 0;
    for (auto& s : l_samples)
    {
        tot_a_s += s.variations.size();
    }
    spdlog::info("Samples size: {} GB", inGigabytes(tot_a_s * 8));
}

//------------------------------------------------------------------------------

void
vcfbwt::VCF::init_vcf(const std::string &vcf_path, std::size_t i)
{
    init_vcf(vcf_path, variations, samples, samples_id, i);
}

//------------------------------------------------------------------------------


void
vcfbwt::VCF::init_multi_ref(const std::vector<std::string>& refs_path)
{
    if (refs_path.size() == 0) { spdlog::error("No reference file provided"); std::exit(EXIT_FAILURE); }
    
    spdlog::info("Opening {} ref files, assuming input order reflects the intended genome order", refs_path.size());
    for (std::size_t i = 0; i < refs_path.size(); i++) { init_ref(refs_path[i], i == (refs_path.size() - 1)); }
}

//------------------------------------------------------------------------------

void
vcfbwt::VCF::init_multi_vcf(const std::vector<std::string>& vcfs_path)
{
    if (vcfs_path.empty()) { spdlog::error("No vcf file provided"); std::exit(EXIT_FAILURE); }
    
    spdlog::info("Opening {} vcf files, assuming input order reflects the intended genome order", vcfs_path.size());

    std::vector<std::vector<Sample>> tmp_samples_array;
    std::vector<std::vector<Variation>> tmp_variations_array;
    std::vector<std::unordered_map<std::string, std::size_t>> tmp_samples_id;

    tmp_samples_array.resize(vcfs_path.size());
    tmp_variations_array.resize(vcfs_path.size());
    tmp_samples_id.resize(vcfs_path.size());

    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < vcfs_path.size(); i++)
    {
        init_vcf(vcfs_path[i],
                 tmp_variations_array[i],
                 tmp_samples_array[i],
                 tmp_samples_id[i],
                 i);
    }

    // Merge tmp structures into global structures
    spdlog::info("Merging variations");
    for (std::size_t i = 0; i < vcfs_path.size(); i++)
    {
        std::size_t prev_variations_arr_size = variations.size();
        this->variations.insert(this->variations.end(), tmp_variations_array[i].begin(), tmp_variations_array[i].end());
        tmp_variations_array[i].clear();

        for (auto& sample : tmp_samples_array[i])
        {
            if (this->samples_id.find(sample.id()) ==this-> samples_id.end())
            {
                Sample s(sample.id(), this->reference, this->variations);
                this->samples.push_back(s);
                this->samples_id.insert(std::make_pair(sample.id(), this->samples.size() - 1));
            }

            for (auto& sample_variation : sample.variations)
            {
                this->samples[samples_id[sample.id()]].variations.push_back(sample_variation + prev_variations_arr_size);
            }
        }
        tmp_samples_array[i].clear();
        tmp_samples_id[i].clear();
    }

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


