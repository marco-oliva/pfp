//
//  pfp_algo.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#include <pfp_algo.hpp>

//------------------------------------------------------------------------------

bool
vcfbwt::pfp::Dictionary::contains(const std::string& phrase) const
{
    hash_type phrase_hash = string_hash(&(phrase[0]), phrase.size());
    const auto& ptr = hash_string_map.find(phrase_hash);
    
    return ((ptr != hash_string_map.end()) and (ptr->second.phrase == phrase));
}

vcfbwt::hash_type
vcfbwt::pfp::Dictionary::add(const std::string& phrase)
{
    this->sorted = false;
    
    hash_type phrase_hash = string_hash(&(phrase[0]), phrase.size());
    if (hash_string_map.find(phrase_hash) != hash_string_map.end())
    {
        spdlog::error("Hash collision! Hash already in the dictionary");
        std::exit(EXIT_FAILURE);
    }
    
    DictionaryEntry entry(phrase);
    hash_string_map.insert(std::make_pair(phrase_hash, entry));
    
    if (this->size() == std::numeric_limits<size_type>::max())
    { spdlog::error("Dictionary too big for {}", typeid(size_type).name()); std::exit(EXIT_FAILURE); }
    
    return phrase_hash;
}

void
vcfbwt::pfp::Dictionary::update_frequency(const std::string& phrase)
{
    if (this->contains(phrase))
    {
        this->hash_string_map[string_hash(&(phrase[0]), phrase.size())].occurrences++;
    }
}

vcfbwt::hash_type
vcfbwt::pfp::Dictionary::get(const std::string& phrase) const
{
    return string_hash(&(phrase[0]), phrase.size());
}

bool
ref_smaller(
std::pair<std::reference_wrapper<std::string>, vcfbwt::hash_type> a,
std::pair<std::reference_wrapper<std::string>, vcfbwt::hash_type> b)
{
    return (a.first.get() < b.first.get());
}

void
vcfbwt::pfp::Dictionary::sort()
{
    // sort the dictionary
    for (auto& entry : this->hash_string_map)
    {
        this->sorted_phrases.push_back(std::make_pair(std::ref(entry.second.phrase), entry.first));
    }
    std::sort(sorted_phrases.begin(), sorted_phrases.end(), ref_smaller);
    
    // insert in hashmap
    this->hash_to_ranks.clear();
    for (size_type i = 0; i < sorted_phrases.size(); i++)
    {
        hash_to_ranks.insert(std::make_pair(sorted_phrases[i].second, i + 1)); // 1 based
    }
    
    this->sorted = true;
}

vcfbwt::size_type
vcfbwt::pfp::Dictionary::hash_to_rank(hash_type hash)
{
    if (not this->sorted) { sort(); }
    auto it = hash_to_ranks.find(hash);
    if (it != hash_to_ranks.end()) { return it->second; }
    else { spdlog::error("Something went wrong"); std::exit(EXIT_FAILURE); }
}

//------------------------------------------------------------------------------

void
vcfbwt::pfp::ReferenceParse::init(const std::string& reference, const std::unordered_set<hash_type>& ts)
{
    std::string phrase;
    spdlog::info("Parsing reference");
    
    // Reference as first sample, just one dollar to be compatible with Giovanni's pscan.cpp
    phrase.append(1, DOLLAR);
    
    for (std::size_t ref_it = 0; ref_it < reference.size(); ref_it++)
    {
        char c = reference[ref_it];
        
        phrase.push_back(c);
        
        if (
        (phrase.size() > this->params.w) and
        (((string_hash(&(phrase[phrase.size() - this->params.w]), this->params.w) % this->params.p) == 0) or
        (ts.find(string_hash(&(phrase[phrase.size() - this->params.w]), this->params.w)) != ts.end())))
        {
            hash_type hash = 0;
            if (this->dictionary.contains(phrase))    { hash = this->dictionary.get(phrase); }
            else                                      { hash = this->dictionary.add(phrase); }
            
            this->parse.push_back(hash);
            this->trigger_strings_position.push_back(ref_it - this->params.w + 1);
    
            phrase.erase(phrase.begin(), phrase.end() - this->params.w); // Keep the last w chars
        }
    }
    
    // Last phrase
    if (phrase.size() >= this->params.w)
    {
        // Append w dollar prime at the end, reference as the first sample
        phrase.append(this->params.w, DOLLAR_PRIME);
    
        hash_type hash = 0;
        if (this->dictionary.contains(phrase))    { hash = this->dictionary.get(phrase); }
        else                                      { hash = this->dictionary.add(phrase); }
    
        this->parse.push_back(hash);
        this->trigger_strings_position.push_back(reference.size() - 1);
    }
    else { spdlog::error("The reference doesn't have w dollar prime at the end!"); std::exit(EXIT_FAILURE); }
}

//------------------------------------------------------------------------------

void
vcfbwt::pfp::Parser::init(const Params& params, const std::string& prefix, ReferenceParse& rp, std::size_t t)
{
    this->w = params.w; this->out_file_prefix = prefix; this->p = params.p; this->tags = t; this->parse_size = 0;
    if (not ((tags & MAIN) or (tags & WORKER))) { spdlog::error("A parser must be either the main parser or a worker"); std::exit(EXIT_FAILURE); }
    
    if (tags & MAIN) {this->out_file_name = out_file_prefix + ".parse"; }
    if ((tags & MAIN) and params.compress_dictionary) { tags = tags | COMPRESSED; }
    this->tmp_out_file_name = TempFile::getName("parse");
    this->out_file.open(tmp_out_file_name, std::ios::binary);
    this->reference_parse = &rp;
    
    this->params = params;
}

void
vcfbwt::pfp::Parser::operator()(const vcfbwt::Sample& sample, const std::unordered_set<hash_type>& ts)
{
    if (ts.empty()) {spdlog::error("Empty set of treigger strings!"); std::exit(EXIT_FAILURE); }
    
    Sample::iterator sample_iterator(sample);
    
    std::string phrase;
    
    // Every sample starts with w dollar prime
    phrase.append(this->w, DOLLAR_PRIME);
    
    // Shorthands
    std::vector<size_type>& tsp = reference_parse->trigger_strings_position;
    
    std::size_t start_window = 0, end_window = 0;
    while (not sample_iterator.end())
    {
        // Compute where we are on the reference
        std::size_t pos_on_reference = sample_iterator.get_ref_it();
        
        if ( not ((sample_iterator.get_var_it() > 0) and (sample_iterator.prev_variation() > (pos_on_reference - (4 * this->w)))))
        {
            // Set start postion to the position in the reference parse after the last computed phrase
            if (params.use_acceleration and ((phrase.size() == this->w) and ((pos_on_reference != 0) and (phrase[0] != DOLLAR_PRIME))))
            {
                start_window = end_window;
                while ((tsp[start_window] + this->w) <= pos_on_reference) { start_window++; }
    
                // Iterate over the parse up to the next variation
                while (tsp[end_window + 1] < (sample_iterator.next_variation() - (this->w + 1))) { end_window++; }
                
                // If the window is not empty
                if ((start_window < end_window) and (tsp[end_window] > pos_on_reference))
                {
                    spdlog::debug("------------------------------------------------------------");
                    spdlog::debug("from {}", sample.get_reference().substr(tsp[start_window - 1], this->w));
                    spdlog::debug("copied from {} to {}", tsp[start_window], tsp[end_window] + this->w);
                    spdlog::debug("next variatin: {}", sample_iterator.next_variation());
                    
                    
                    // copy from parse[start_window : end_window]
                    out_file.write((char*) &(this->reference_parse->parse[start_window]), sizeof(hash_type) * (end_window - start_window + 1));
                    parse_size += end_window - start_window + 1;
            
                    // move itarators and re initialize phrase
                    sample_iterator.go_to(tsp[end_window]);
                    phrase.clear();
                    for (std::size_t i = 0; i < this->w; i++) { ++sample_iterator; phrase.push_back(*sample_iterator);}
                    ++sample_iterator;
                    spdlog::debug("New phrase [{}]: {}", phrase.size(), phrase);
                    spdlog::debug("------------------------------------------------------------");
                }
            }
        }
        
        // Next phrase should contain a variation so parse as normal, also if we don't
        // want to use the acceleration we should always end up here
        phrase.push_back(*sample_iterator);
        ++sample_iterator;
    
        if (
        (phrase.size() > this->params.w) and
        ((((string_hash(&(phrase[phrase.size() - this->params.w]), this->w)) % this->params.p) == 0) or
        (ts.find(string_hash(&(phrase[phrase.size() - this->params.w]), this->w)) != ts.end())))
        {
            hash_type hash = 0;
            if (reference_parse->dictionary.contains(phrase))   { hash = reference_parse->dictionary.get(phrase); }
            else if (dictionary.contains(phrase))               { hash = dictionary.get(phrase); }
            else                                                { hash = dictionary.add(phrase); }
        
            out_file.write((char*) (&hash), sizeof(hash_type)); parse_size++;
    
            if (phrase[0] != DOLLAR_PRIME)
            {
                spdlog::debug("------------------------------------------------------------");
                spdlog::debug("Parsed phrase [{}] {}", phrase.size(), phrase);
                spdlog::debug("------------------------------------------------------------");
            }
            
            phrase.erase(phrase.begin(), phrase.end() - this->w); // Keep the last w chars
        }
    }
    
    // Last phrase
    if (phrase.size() >= this->w)
    {
        // Append w dollar prime at the end of each sample, DOLLAR if it's the last sample
        if (this->tags & LAST) { phrase.append(this->w, DOLLAR); }
        else { phrase.append(this->w, DOLLAR_PRIME); }
        
        hash_type hash = 0;
        if (reference_parse->dictionary.contains(phrase))   { hash = reference_parse->dictionary.get(phrase); }
        else if (dictionary.contains(phrase))               { hash = dictionary.get(phrase); }
        else                                                { hash = dictionary.add(phrase); }
        
        out_file.write((char*) (&hash), sizeof(hash_type));  parse_size++;
    }
    else { spdlog::error("A sample doesn't have w dollar prime at the end!"); std::exit(EXIT_FAILURE); }
    
}

void
vcfbwt::pfp::Parser::close()
{
    if (closed) return; closed = true;
    
    if ((tags & MAIN) or (tags & WORKER))
    {
        this->out_file.close();
    }
    
    // Output parse, substitute hash with rank
    if (tags & MAIN)
    {
        // close all the registered workers and merge their dictionaries
        spdlog::info("Main parser: closing all registered workers");
        for (auto worker : registered_workers) { worker.get().close(); }
    
        spdlog::info("Main parser: Merging reference dictionary into MAIN dictionary, WARNING (TODO) does not add up occurrences");
        this->dictionary.hash_string_map.insert(this->reference_parse->dictionary.hash_string_map.begin(), this->reference_parse->dictionary.hash_string_map.end());
    
        spdlog::info("Main parser: merging workers dictionaries, WARNING (TODO) does not add up occurrences");
        for (auto worker : registered_workers)
        {
            // there should not be any collision since the reference dictionary is checked during parsing
            // if there is a collision the hash should be the same anyway
            this->dictionary.hash_string_map.insert(worker.get().dictionary.hash_string_map.begin(), worker.get().dictionary.hash_string_map.end());
        }
        
        spdlog::info("Main parser: Replacing hash values with ranks in MAIN, WORKERS and reference");
    
        // MAIN mmap file and substitute
        if (this->parse_size != 0)
        {
            std::error_code error;
            mio::mmap_sink rw_mmap = mio::make_mmap_sink(tmp_out_file_name, 0, mio::map_entire_file, error);
            if (error) { spdlog::error(error.message()); std::exit(EXIT_FAILURE); }
    
            for (size_type i = 0; i < (rw_mmap.size() / sizeof(hash_type)); i++)
            {
                hash_type hash;
                std::memcpy(&hash, rw_mmap.data() + (i * sizeof(hash_type)), sizeof(hash_type));
                size_type rank = this->dictionary.hash_to_rank(hash);
                std::memcpy(rw_mmap.data() + (i * sizeof(size_type)), &rank, sizeof(size_type));
            }
            rw_mmap.unmap();
            truncate_file(tmp_out_file_name, this->parse_size * sizeof(size_type));
        }
        
        // repeat for reference
        if (not this->reference_parse->parse.empty())
        {
            for (size_type i = 0; i < this->reference_parse->parse.size(); i++)
            {
                hash_type rank = this->dictionary.hash_to_rank(this->reference_parse->parse[i]);
                this->reference_parse->parse[i] = rank;
            }
        }
        
        // repeat for every worker
        for (auto worker : registered_workers)
        {
            if (worker.get().parse_size != 0)
            {
                std::error_code error;
                mio::mmap_sink rw_mmap = mio::make_mmap_sink(worker.get().tmp_out_file_name, 0, mio::map_entire_file, error);
                if (error) { spdlog::error(error.message()); std::exit(EXIT_FAILURE); }
    
                for (size_type i = 0; i < (rw_mmap.size() / sizeof(hash_type)); i++)
                {
                    hash_type hash;
                    std::memcpy(&hash, rw_mmap.data() + (i * sizeof(hash_type)), sizeof(hash_type));
                    size_type rank = this->dictionary.hash_to_rank(hash);
                    std::memcpy(rw_mmap.data() + (i * sizeof(size_type)), &rank, sizeof(size_type));
                }
                rw_mmap.unmap();
                truncate_file(worker.get().tmp_out_file_name, worker.get().parse_size * sizeof(size_type));
            }
        }
    
        // Merging files toghether
        spdlog::info("Main parser: concatenating parsings from workers and reference, reference as first");
        std::ofstream merged(out_file_name, std::ios_base::binary);
        
        // Reference
        size_type out_parse_size = 0;
        out_parse_size += this->reference_parse->parse.size();
        for (auto& e : this->reference_parse->parse)
        {
            size_type out_e = e;
            merged.write((char*) &out_e, sizeof(size_type));
        }
    
        //Main
        if ((this->parse_size) != 0)
        {
            std::ifstream main_parse(this->tmp_out_file_name);
            merged << main_parse.rdbuf();
            out_parse_size += this->parse_size;
        }
        
        // Workers
        for (auto worker : registered_workers)
        {
            if (worker.get().parse_size != 0)
            {
                out_parse_size += worker.get().parse_size;
                std::ifstream worker_parse(worker.get().tmp_out_file_name);
                merged << worker_parse.rdbuf();
            }
        }
        merged.close();
        
        this->parse_size = out_parse_size;
        
        // Print dicitionary on disk
        if (tags & UNCOMPRESSED)
        {
            spdlog::info("Main parser: writing dictionary to disk NOT COMPRESSED");
            std::string dict_file_name = out_file_prefix + ".dict";
            std::ofstream dict(dict_file_name);
        
            for (size_type i = 0; i < this->dictionary.size(); i++)
            {
                dict.write(this->dictionary.sorted_entry_at(i).c_str(), this->dictionary.sorted_entry_at(i).size());
                dict.put(ENDOFWORD);
            }
        
            dict.put(ENDOFDICT);
            dict.close();
        }
        
        if (tags & COMPRESSED)
        {
            spdlog::info("Main parser: writing dictionary to disk COMPRESSED");
            std::ofstream dicz(out_file_prefix + ".dicz");
            std::ofstream lengths(out_file_prefix + ".dicz.len");
    
            for (size_type i = 0; i < this->dictionary.size(); i++)
            {
                std::size_t shift = 0;
                if (i != 0) { shift = this->w; }
                dicz.write(this->dictionary.sorted_entry_at(i).c_str() + shift,
                           this->dictionary.sorted_entry_at(i).size() - shift);
                int32_t len = this->dictionary.sorted_entry_at(i).size() - shift;
                lengths.write((char*) &len, sizeof(int32_t));
            }
            dicz.close();
            lengths.close();
        }
    
        // Outoput Occurrencies
        // NOP
    }
}

//------------------------------------------------------------------------------

void
vcfbwt::pfp::Parser::compute_trigger_strings(vcfbwt::VCF& vcf, const Params& params, std::unordered_set<hash_type>& trigger_string_set)
{
    // Put last w charachters of each sample in the VCF as a trigger strings
    std::string dollar_window; dollar_window.append(params.w, DOLLAR_PRIME);
    trigger_string_set.insert(string_hash(&(dollar_window[0]), dollar_window.size()));
    
    // Compute trigger strings from most common variations
    if (params.compute_seeded_trigger_strings)
    {
        const std::vector<Variation>& variations_ref = vcf.get_variations();
    
        for (std::size_t i = 0; i < variations_ref.size(); i++)
        {
            const Variation& variation = variations_ref[i];
        
            if (variation.freq > params.min_frequency)
            {
                std::size_t region_start_pos = variation.pos - 1;
                std::size_t region_end_pos = 0;
            
                // Keep going forward until we have a window without variations over a certain frequency threshold
                // and the region is at least of a certain length
                std::size_t j = i + 1;
                std::size_t last_over_frequency_variation_pos = variation.pos;
                region_end_pos = std::min(region_start_pos + params.min_seed_region_length(), variations_ref[j].pos);
            
                while (
                (region_start_pos - region_end_pos <= params.min_seed_region_length()) or
                (region_end_pos - last_over_frequency_variation_pos < params.w)
                )
                {
                    region_end_pos = variations_ref[j].pos + 1;
                    if (variations_ref[j].freq >= params.min_frequency)
                    {
                        last_over_frequency_variation_pos = variations_ref[j].pos;
                    }
                
                    ++j;
                }
            
                // Extract seed string from region limits
                std::string opening_trigger_string = vcf.get_reference().substr(region_start_pos - params.w, params.w);
                std::string closing_trigger_string = vcf.get_reference().substr(region_end_pos, params.w);
                
                trigger_string_set.insert(string_hash(&(opening_trigger_string[0]), opening_trigger_string.size()));
                trigger_string_set.insert(string_hash(&(closing_trigger_string[0]), closing_trigger_string.size()));
            
                // Don't process variations in seed region more than once
                i = j;
            }
        }
        
        spdlog::info("Generated {} trigger strings", trigger_string_set.size());
    }
    else
    {
        spdlog::info("Not generating seeded trigger strings");
    }
}

//------------------------------------------------------------------------------

void
vcfbwt::pfp::Parser::read_parse(std::string parse_file_name, std::vector<size_type>& parse)
{
    std::ifstream parse_file(parse_file_name, std::ios::binary);
    if (not parse_file.is_open()) { spdlog::error("ERROR!"); std::exit(EXIT_FAILURE); }
    
    while (not parse_file.eof()) { size_type i; parse_file.read((char*) &i, sizeof(size_type)); parse.push_back(i); }
    
    // remove last entry, apparently it reads the last twice, strano
    parse.pop_back();
}

void
vcfbwt::pfp::Parser::read_dictionary(std::string dic_file_name, std::vector<std::string>& dictionary_vector)
{
    std::ifstream dic_file(dic_file_name);
    if (not dic_file.is_open()) { spdlog::error("ERROR!"); std::exit(EXIT_FAILURE); }
    
    char c = '\5';
    while (c != ENDOFDICT)
    {
        std::string phrase = "";
        while ((c = dic_file.get()) != ENDOFWORD and c != ENDOFDICT)
        {
            phrase.push_back(c);
        }
        if (phrase.size() != 0) { dictionary_vector.push_back(phrase); }
    }
}

void
vcfbwt::pfp::Parser::merge(std::string left_prefix, std::string right_prefix, std::string out_prefix, const Params& params)
{
    // Merged Dictionary and Parse
    Dictionary dictionary;
    
    // Read Left and Right Dictionary
    spdlog::info("Loading dictionaries from disk");
    std::vector<std::string> left_dictionary;   read_dictionary(left_prefix + ".dict", left_dictionary);
    std::vector<std::string> right_dictionary; read_dictionary(right_prefix + ".dict", left_dictionary);
    
    // Read Left parse and substitute ranks with hash again storing in the merged dictionary
    spdlog::info("Adjusting ranks");
    std::error_code error;
    mio::mmap_sink rw_mmap = mio::make_mmap_sink(left_prefix + ".parse", 0, mio::map_entire_file, error);
    if (error) { spdlog::error(error.message()); std::exit(EXIT_FAILURE); }
    
    for (size_type i = 0; i < (rw_mmap.size() / sizeof(hash_type)); i++)
    {
        hash_type rank;
        std::memcpy(&rank, rw_mmap.data() + (i * sizeof(hash_type)), sizeof(hash_type));
        
        std::string& phrase = left_dictionary[rank];
        hash_type hash = 0;
        if (dictionary.contains(phrase))    { hash = dictionary.get(phrase); }
        else                                { hash = dictionary.add(phrase); }
        
        std::memcpy(rw_mmap.data() + (i * sizeof(hash_type)), &hash, sizeof(hash_type));
    }
    rw_mmap.unmap(); left_dictionary.resize(0);
    
    // Read right parse, changing first phrase;
    rw_mmap = mio::make_mmap_sink(left_prefix + ".parse", 0, mio::map_entire_file, error);
    if (error) { spdlog::error(error.message()); std::exit(EXIT_FAILURE); }
    
    for (size_type i = 0; i < (rw_mmap.size() / sizeof(hash_type)); i++)
    {
        hash_type rank;
        std::memcpy(&rank, rw_mmap.data() + (i * sizeof(hash_type)), sizeof(hash_type));
        
        std::string& phrase = right_dictionary[rank];
        
        if (rank == 1) { phrase[0] = DOLLAR_PRIME; phrase.insert( 0, params.w - 1, DOLLAR_PRIME); }
        
        hash_type hash = 0;
        if (dictionary.contains(phrase))    { hash = dictionary.get(phrase); }
        else                                { hash = dictionary.add(phrase); }
        
        std::memcpy(rw_mmap.data() + (i * sizeof(hash_type)), &hash, sizeof(hash_type));
    }
    rw_mmap.unmap(); right_dictionary.resize(0);
    
    // Sort the Dictionary
    spdlog::info("Sorting merged dictionary");
    std::vector<std::pair<std::string, hash_type>> entries;
    for (auto& entry : dictionary.hash_string_map)
    {
        entries.push_back(std::make_pair(entry.second.phrase, entry.first));
    }
    std::sort(entries.begin(), entries.end());
    
    // Insert in hashmap
    std::unordered_map<hash_type, hash_type> hash_to_ranks;
    for (hash_type i = 0; i < entries.size(); i++)
    {
        hash_to_ranks.insert(std::make_pair(entries[i].second, i + 1)); // 1 based
    }
    
    // Merge parsings
    spdlog::info("Creating output parse");
    std::ofstream out_parse(out_prefix + ".parse", std::ios::binary);
    
    // Left
    mio::mmap_source rl_mmap; rl_mmap.map(left_prefix + ".parse", error);
    if (error) { spdlog::error(error.message()); std::exit(EXIT_FAILURE); }
    for (size_type i = 0; i < (rl_mmap.size() / sizeof(hash_type)); i++)
    {
        hash_type hash;
        std::memcpy(&hash, rl_mmap.data() + (i * sizeof(hash_type)), sizeof(hash_type));
    
        auto it = hash_to_ranks.find(hash);
        if (it != hash_to_ranks.end()) { out_parse.write((char*) &(it->second), sizeof(hash_type)); }
        else { spdlog::error("Something wrong happend"); std::exit(EXIT_FAILURE); }
    }
    rl_mmap.unmap();
    
    // Right
    mio::mmap_source rr_mmap; rr_mmap.map(left_prefix + ".parse", error);
    if (error) { spdlog::error(error.message()); std::exit(EXIT_FAILURE); }
    for (size_type i = 0; i < (rr_mmap.size() / sizeof(hash_type)); i++)
    {
        hash_type hash;
        std::memcpy(&hash, rr_mmap.data() + (i * sizeof(hash_type)), sizeof(hash_type));
        
        auto it = hash_to_ranks.find(hash);
        if (it != hash_to_ranks.end()) { out_parse.write((char*) &(it->second), sizeof(hash_type)); }
        else { spdlog::error("Something wrong happend"); std::exit(EXIT_FAILURE); }
    }
    rr_mmap.unmap();
    out_parse.close();
    
    // Print dicitionary on disk
    spdlog::info("Writing dictionary to disk NOT COMPRESSED");
    std::string dict_file_name = out_prefix + ".dict";
    std::ofstream dict(dict_file_name);
    
    for (auto& entry : entries)
    {
        dict.write(entry.first.c_str(), entry.first.size());
        dict.put(ENDOFWORD);
    }
    
    dict.put(ENDOFDICT);
    dict.close();
}