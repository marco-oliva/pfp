//
//  vcf.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#include <vcf.hpp>
#include <pfp_algo.hpp>

const std::string vcfbwt::VCF::vcf_freq = "AF";


//------------------------------------------------------------------------------

vcfbwt::Sample::iterator::iterator(const Sample& sample, std::size_t genotype) :
sample_(sample), genotype(genotype),
ref_it_(0), sam_it_(0), var_it_(0), curr_var_it_(0), prev_variation_it(0),
curr_char_(NULL), sample_length_(sample.reference_.size())
{
    // Compute sample length, might take some time
    long long int indels = 0; // could be negative, so int
    for (std::size_t i = 0; i < this->sample_.variations.size(); i++)
    {
        std::size_t var_id = this->sample_.variations[i];
        int var_genotype = this->sample_.genotypes[i][this->genotype];
        if (var_genotype != 0)
        {
            indels += sample_.variations_list[var_id].alt[var_genotype].size() -
                      sample_.variations_list[var_id].ref_len;
        }
        
    }
    sample_length_ = sample_length_ + indels;
    while (var_it_ < sample_.variations.size() and sample_.genotypes[var_it_][genotype] == 0)
    { var_it_++; }
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
    return sample_.get_variation(prev_variation_it).pos;
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
        int var_genotype = sample_.genotypes[var_it_][genotype];
        // Handling same position insertions see bcftools consensus:
        // https://github.com/samtools/bcftools/blob/df43fd4781298e961efc951ba33fc4cdcc165a19/consensus.c#L723
        int gap = ref_it_ - curr_variation.pos;
        if (gap > curr_var_it_)
        {
            // Check length of unchanged bases
            int start = 0;
            int len = curr_variation.alt[var_genotype].size();
            assert(len >= gap);
            while (start < std::min(gap, len) && curr_variation.alt[var_genotype][start] == curr_variation.alt[0][start])
                ++start;
            if (start < gap)
            {
                spdlog::error("vcfbwt::Contig::iterator::operator++() Error: variant overriding previous variant.");
                std::exit(EXIT_FAILURE);
            }
            // We need to preserve preceeding variants
            curr_var_it_ = gap;
        }

        bool get_next_variant = true;

        if (curr_var_it_ < curr_variation.alt[var_genotype].size() - 1)
        {
            curr_char_ = &curr_variation.alt[var_genotype][curr_var_it_];
            curr_var_it_++;
            get_next_variant = false;
        }
        else if (curr_var_it_ < curr_variation.alt[var_genotype].size())
            curr_char_ = &curr_variation.alt[var_genotype].back();
        else
            curr_char_ = &(sample_.reference_[ref_it_++ + curr_variation.ref_len - gap]);

        sam_it_++;
        if (get_next_variant)
        {
            prev_variation_it = var_it_;
            var_it_++;
            while (var_it_ < sample_.variations.size() and sample_.genotypes[var_it_][genotype] == 0)
            {
                var_it_++;
            }
            curr_var_it_ = 0;
            ref_it_ += curr_variation.ref_len - gap; // Adding -gap to balance the sipping
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
vcfbwt::VCF::init_samples(const std::string& samples_path)
{
    spdlog::info("Reading samples file: {}", samples_path);
    std::ifstream in_stream(samples_path);
    if (not in_stream.is_open()) { spdlog::error("Error while opening {}", samples_path); }

    for (std::string line; getline( in_stream, line);)
    {
        this->input_samples.insert(line);
    }
    spdlog::info("Done reading {}, found {} samples", samples_path, this->input_samples.size());

    for (auto& sample_id : input_samples) { spdlog::info(sample_id); }
}

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

void
vcfbwt::VCF::init_vcf(const std::string& vcf_path, std::vector<Variation>& l_variations,
                      std::vector<Sample>& l_samples, std::unordered_map<std::string, std::size_t>& l_samples_id,
                      std::size_t i)
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
    
    // get l_samples ids from header
    std::size_t n_samples = bcf_hdr_nsamples(hdr);
    if (this->max_samples == 0) { set_max_samples(n_samples); }

    std::size_t size_before = l_samples.size();
    for (std::size_t i = 0; i < std::min(n_samples, this->max_samples); i++)
    {
        vcfbwt::Sample s(std::string(hdr->samples[i]), this->reference, l_variations);
        if (l_samples_id.find(s.id()) == l_samples_id.end())
        {
            l_samples.push_back(s);
            l_samples_id.insert(std::make_pair(s.id(), i));
        }
    }
    spdlog::debug("{} new l_samples in the vcf, tot: {}", l_samples.size() - size_before, l_samples.size());
    
    // struct for storing each record
    bcf1_t *rec = bcf_init();
    if (rec == NULL)
    {
        spdlog::error("Error while parsing vcf file: {}", vcf_path);
        bcf_close(inf);
        bcf_hdr_destroy(hdr);
        std::exit(EXIT_FAILURE);
    }

    std::vector<std::vector<int>> tppos(1, std::vector<int>(n_samples,0));
    std::vector<std::vector<bool>> prev_is_ins(1, std::vector<bool>(n_samples,false));
    
    // start parsing
    while (bcf_read(inf, hdr, rec) == 0)
    {
        vcfbwt::Variation var;
        
        var.ref_len = rec->rlen;
        std::size_t offset = i != 0 ? ref_sum_lengths[i-1] : 0; // when using multiple vcfs
        var.pos = rec->pos + offset;
        var.freq = 0;
        
        // get all alternate alleles
        bcf_unpack(rec, BCF_UN_ALL);
        int type = bcf_get_variant_types(rec);
        for (int allele_idx = 0; allele_idx < rec->n_allele; allele_idx++)
        {
            var.alt.push_back(rec->d.allele[allele_idx]);
        }
        
        int32_t *gt_arr = NULL, ngt_arr = 0;
        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if ( ngt > 0 )
        {
            int max_ploidy = ngt/n_samples;
            while (max_ploidy > tppos.size())
            {
                tppos.push_back(std::vector<int>(n_samples,0));
                prev_is_ins.push_back(std::vector<bool>(n_samples,false));
            }
            bool skip_this_variation = false;
            for (std::size_t i_s = 0; i_s < n_samples; i_s++)
            {
                if (skip_this_variation) { break; }
                int32_t *ptr = gt_arr + i_s * max_ploidy;
                std::vector<int> alleles_idx(max_ploidy, 0);
                bool alt_alleles_set = false;
                for (std::size_t j = 0; j < max_ploidy; j++)
                {
                    // if true, the sample has smaller ploidy
                    if ( ptr[j]==bcf_int32_vector_end ) { break; }

                    // missing allele
                    if ( bcf_gt_is_missing(ptr[j]) ) { continue; }
                    
                    if (bcf_gt_allele(ptr[j]))
                    {
                        // the VCF 0-based allele index
                        int allele_index = bcf_gt_allele(ptr[j]);
                        int var_type = bcf_get_variant_type(rec, allele_index);

                        // Determine if overlap. Logic copied from leviosam's: 
                        // https://github.com/alshai/levioSAM/blob/f72d84ad1141c84e4b315c0dc5d705d2c0d5b936/src/leviosam.hpp#L530
                        // copied from bcftools consensus`:
                        // https://github.com/samtools/bcftools/blob/df43fd4781298e961efc951ba33fc4cdcc165a19/consensus.c#L579

                        // For some variant types POS+REF refer to the base *before* the event; in such case set trim_beg
                        int trim_beg = 0;
                        int var_len  = rec->d.var[allele_index].n;
                        if ( var_type & VCF_INDEL ) trim_beg = 1;
                        else if ( (var_type & VCF_OTHER) && !strcasecmp(rec->d.allele[allele_index],"<DEL>") ) {
                            trim_beg = 1;
                            var_len  = 1 - var.ref_len;
                        }
                        else if ( (var_type & VCF_OTHER) && !strncasecmp(rec->d.allele[allele_index],"<INS",4) )
                            trim_beg = 1;

                        if (rec->pos <= tppos[j][i_s]) {
                            int overlap = 0;
                            if ( rec->pos < tppos[j][i_s] || !trim_beg || var_len==0 || prev_is_ins[j][i_s] ) overlap = 1;
                            if (overlap) {
                                spdlog::debug("vcfbwt::VCF::init_vcf: Skipping overlapping variantat sample {} in pos {}", i_s, var.pos);
                                continue;
                            }
                        }

                        // Skip symbolic allele
                        if (var.alt[allele_index][0] == '<')
                        {
                            spdlog::debug("vcfbwt::VCF::init_vcf: Skipping symbolic allele at pos {}", var.pos);
                            skip_this_variation = true;
                            continue;
                        }
                        // Update tppos and prev_is_ins
                        tppos[j][i_s] = rec->pos + rec->rlen - 1;
                        prev_is_ins[j][i_s] = (var.alt[0].size() < var.alt[allele_index].size());

                        alleles_idx[j] = allele_index;
                        alt_alleles_set = true;
                    }
                }
                
                if (alt_alleles_set)
                {
                    auto id = l_samples_id.find(std::string(hdr->samples[i_s]));
                    if ((id != l_samples_id.end() and id->second < max_samples) and
                    ((input_samples.empty()) or input_samples.contains(id->first))) // Process only wanted l_samples
                    {
                        // Update frequency, to be normalized by the number of samples when parsing ends
                        var.freq += 1;
                        var.used = true;
                        // Add variation to sample, size() because we have not added the variations to the list yet
                        l_samples[id->second].variations.emplace_back(l_variations.size());
                        l_samples[id->second].genotypes.emplace_back(alleles_idx);
                    }
                }
            }
        }
        if (var.used) { l_variations.push_back(var); }
        free(gt_arr);
    }
    
    // free allocated memory
    bcf_hdr_destroy(hdr);
    bcf_close(inf);
    bcf_destroy(rec);
    
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
    std::size_t tot_a_s = 0, tot_samples = 0;
    for (auto& s : this->samples)
    {
        tot_a_s += s.variations.size();
        if (s.variations.size() > 0) { tot_samples += 1; }
    }
    spdlog::info("Average variations per sample: {}", tot_a_s / tot_samples);
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
            if (this->samples_id.find(sample.id()) == this->samples_id.end())
            {
                Sample s(sample.id(), this->reference, this->variations);
                this->samples.push_back(s);
                this->samples_id.insert(std::make_pair(sample.id(), this->samples.size() - 1));
            }

            for (auto& sample_variation : sample.variations)
            {
                this->samples[samples_id[sample.id()]].variations.push_back(sample_variation + prev_variations_arr_size);
            }
    
            for (auto& genotype_info : sample.genotypes)
            {
                this->samples[samples_id[sample.id()]].genotypes.push_back(genotype_info);
            }
        }
        tmp_samples_array[i].clear();
        tmp_samples_id[i].clear();
    }

    // print some statistics
    spdlog::info("Variations size [{}]: {}GB", variations.size(), inGigabytes(variations.size() * sizeof(Variation)));
    spdlog::info("Reference size: {} GB", inGigabytes(reference.size()));
    
    std::size_t tot_a_s = 0, tot_samples = 0;
    for (auto& s : this->samples)
    {
        tot_a_s += s.variations.size();
        if (s.variations.size() > 0) { tot_samples += 1; }
    }
    spdlog::info("Samples size: {} GB", inGigabytes(tot_a_s * 8));
    spdlog::info("Average variations per sample: {}", tot_a_s / tot_samples);
}

//------------------------------------------------------------------------------

