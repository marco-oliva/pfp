//
//  pfp_algo.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#include <pfp_algo.hpp>

//------------------------------------------------------------------------------

bool
vcfbwt::pfp::Dictionary::contains(const std::string& phrase)
{
    // lock the dictionary
    std::lock_guard<std::mutex> guard(dictionary_mutex);

    hash_type phrase_hash = string_hash(&(phrase[0]), phrase.size());
    const auto& ptr = hash_string_map.find(phrase_hash);
    
    return ((ptr != hash_string_map.end()) and (ptr->second.phrase == phrase));
}

vcfbwt::hash_type
vcfbwt::pfp::Dictionary::add(const std::string& phrase)
{
    // lock the dictionary
    std::lock_guard<std::mutex> guard(dictionary_mutex);

    this->sorted = false;
    
    hash_type phrase_hash = string_hash(&(phrase[0]), phrase.size());
    if (hash_string_map.contains(phrase_hash))
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

vcfbwt::hash_type
vcfbwt::pfp::Dictionary::check_and_add(const std::string &phrase)
{
    // lock the dictionary
    std::lock_guard<std::mutex> guard(dictionary_mutex);

    // Check if present
    hash_type phrase_hash = string_hash(&(phrase[0]), phrase.size());
    const auto& ptr = hash_string_map.find(phrase_hash);

    if ((ptr != hash_string_map.end()) and (ptr->second.phrase != phrase))
    {
        spdlog::error("Hash collision! Hash already in the dictionary");
        std::exit(EXIT_FAILURE);
    }
    else if (ptr != hash_string_map.end()) { return phrase_hash; }

    this->sorted = false;

    DictionaryEntry entry(phrase);
    hash_string_map.insert(std::make_pair(phrase_hash, entry));

    if (this->size() == std::numeric_limits<size_type>::max())
    { spdlog::error("Dictionary too big for {}", typeid(size_type).name()); std::exit(EXIT_FAILURE); }

    return phrase_hash;
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
    // lock the dictionary
    std::lock_guard<std::mutex> guard(dictionary_mutex);

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

const std::string&
vcfbwt::pfp::Dictionary::sorted_entry_at(std::size_t i)
{
    std::lock_guard<std::mutex> guard(dictionary_mutex);
    if (not this->sorted) { sort(); } return sorted_phrases[i].first.get();
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
vcfbwt::pfp::ReferenceParse::init(const std::string& reference)
{
    if (params.ignore_ts_file.size() > 0)
    {
        spdlog::info("Reading trigger string to be ignored from: {}", params.ignore_ts_file);
        std::ifstream ts_file(params.ignore_ts_file);
        if (not ts_file.is_open()) { spdlog::error("ERROR!"); std::exit(EXIT_FAILURE); }
    
        char c = '\5';
        while (c != ENDOFDICT)
        {
            std::string ts = "";
            while ((c = ts_file.get()) != ENDOFWORD and c != ENDOFDICT)
            {
                ts.push_back(c);
            }
            if (ts.size() != 0) { this->to_ignore_ts_hash.insert(KarpRabinHash::string_hash(ts)); }
        }
        ts_file.close();
        spdlog::info("To be ingored trigger strings: {}", this->to_ignore_ts_hash.size());
    }
    
    std::string phrase;
    spdlog::info("Parsing reference");
    
    // Karp Robin Hash Function for sliding window
    KarpRabinHash kr_hash(this->params.w);
    
    // Reference as first sample, just one dollar to be compatible with Giovanni's pscan.cpp
    phrase.append(1, DOLLAR);
    
    for (std::size_t ref_it = 0; ref_it < reference.size(); ref_it++)
    {
        char c = reference[ref_it];
        
        phrase.push_back(c);
        if (phrase.size() == params.w) { kr_hash.initialize(phrase); }
        else if (phrase.size() > params.w) { kr_hash.update(phrase[phrase.size() - params.w - 1], phrase[phrase.size() - 1]); }
        
        if ((phrase.size() > this->params.w) and ((kr_hash.get_hash() % this->params.p) == 0))
        {
            std::string_view ts(&(phrase[phrase.size() - params.w]), params.w);
            hash_type ts_hash = KarpRabinHash::string_hash(ts);
            if (to_ignore_ts_hash.contains(ts_hash)) { continue; }
            
            hash_type hash = this->dictionary.check_and_add(phrase);
            
            this->parse.push_back(hash);
            this->trigger_strings_position.push_back(ref_it - this->params.w + 1);
    
            phrase.erase(phrase.begin(), phrase.end() - this->params.w); // Keep the last w chars
            
            kr_hash.reset(); kr_hash.initialize(phrase);
        }
    }
    
    // Last phrase
    if (phrase.size() >= this->params.w)
    {
        // Append w dollar prime at the end, reference as the first sample
        phrase.append(this->params.w, DOLLAR_PRIME);
    
        hash_type hash = this->dictionary.check_and_add(phrase);
    
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
    this->dictionary = &this->reference_parse->dictionary;
    
    this->params = params;
}

void
vcfbwt::pfp::Parser::operator()(const vcfbwt::Sample& sample)
{
    Sample::iterator sample_iterator(sample);
    this->samples_processed.push_back(sample.id());
    
    std::string phrase;
    
    // Karp Robin Hash Function for sliding window
    KarpRabinHash kr_hash(this->params.w);
    
    // Every sample starts with w dollar prime
    phrase.append(this->w, DOLLAR_PRIME);
    kr_hash.initialize(phrase);
    
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
                    this->parse_size += end_window - start_window + 1;
            
                    // move itarators and re initialize phrase
                    sample_iterator.go_to(tsp[end_window]);
                    phrase.clear();
                    for (std::size_t i = 0; i < this->w; i++) { ++sample_iterator; phrase.push_back(*sample_iterator);}
                    
                    kr_hash.reset(); kr_hash.initialize(phrase);
                    
                    ++sample_iterator;
                    spdlog::debug("New phrase [{}]: {}", phrase.size(), phrase);
                    spdlog::debug("------------------------------------------------------------");
                }
            }
        }
        
        // Next phrase should contain a variation so parse as normal, also if we don't
        // want to use the acceleration we should always end up here
        phrase.push_back(*sample_iterator);
        kr_hash.update(phrase[phrase.size() - params.w - 1], phrase[phrase.size() - 1]);
        ++sample_iterator;
    
        if ((phrase.size() > this->params.w) and ((kr_hash.get_hash() % this->params.p) == 0))
        {
            std::string_view ts(&(phrase[phrase.size() - params.w]), params.w);
            hash_type ts_hash = KarpRabinHash::string_hash(ts);
            if (this->reference_parse->to_ignore_ts_hash.contains(ts_hash)) { continue; }
            
            hash_type hash = this->dictionary->check_and_add(phrase);
        
            out_file.write((char*) (&hash), sizeof(hash_type)); this->parse_size += 1;
    
            if (phrase[0] != DOLLAR_PRIME)
            {
                spdlog::debug("------------------------------------------------------------");
                spdlog::debug("Parsed phrase [{}] {}", phrase.size(), phrase);
                spdlog::debug("------------------------------------------------------------");
            }
            
            phrase.erase(phrase.begin(), phrase.end() - this->w); // Keep the last w chars
    
            kr_hash.reset(); kr_hash.initialize(phrase);
        }
    }
    
    // Last phrase
    if (phrase.size() >= this->w)
    {
        // Append w dollar prime at the end of each sample, DOLLAR if it's the last sample
        if (this->tags & LAST) { phrase.append(this->w, DOLLAR); }
        else { phrase.append(this->w, DOLLAR_PRIME); }

        hash_type hash = this->dictionary->check_and_add(phrase);
        
        out_file.write((char*) (&hash), sizeof(hash_type));   this->parse_size += 1;
    }
    else { spdlog::error("A sample doesn't have w dollar prime at the end!"); std::exit(EXIT_FAILURE); }
}

void
vcfbwt::pfp::Parser::close()
{
    if (closed) return; closed = true;
    
    if ((tags & MAIN) or (tags & WORKER))
    {
        vcfbwt::DiskWrites::update(out_file.tellp()); // Disk Stats
        this->out_file.close();
    }
    
    // Output parse, substitute hash with rank
    if (tags & MAIN)
    {
        // close all the registered workers and merge their dictionaries
        spdlog::info("Main parser: closing all registered workers");
        for (auto worker : registered_workers) { worker.get().close(); }
        
        // Occurrences
        std::vector<size_type> occurrences(this->dictionary->size(), 0);
        
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
                size_type rank = this->dictionary->hash_to_rank(hash);
                std::memcpy(rw_mmap.data() + (i * sizeof(size_type)), &rank, sizeof(size_type));
                occurrences[rank - 1] += 1;
            }
            rw_mmap.unmap();
            truncate_file(tmp_out_file_name, this->parse_size * sizeof(size_type));
        }
        
        // repeat for reference
        if (not this->reference_parse->parse.empty())
        {
            for (size_type i = 0; i < this->reference_parse->parse.size(); i++)
            {
                hash_type rank = this->dictionary->hash_to_rank(this->reference_parse->parse[i]);
                this->reference_parse->parse[i] = rank;
                occurrences[rank - 1] += 1;
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
                    size_type rank = this->dictionary->hash_to_rank(hash);
                    std::memcpy(rw_mmap.data() + (i * sizeof(size_type)), &rank, sizeof(size_type));
                    occurrences[rank - 1] += 1;
                }
                rw_mmap.unmap();
                truncate_file(worker.get().tmp_out_file_name, worker.get().parse_size * sizeof(size_type));
            }
        }
    
        // Merging files together
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
    
        // Main
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
        vcfbwt::DiskWrites::update(merged.tellp()); // Disk Stats
        merged.close();
        
        this->parse_size = out_parse_size;
        
        // Print dicitionary on disk
        if (tags & UNCOMPRESSED)
        {
            spdlog::info("Main parser: writing dictionary to disk NOT COMPRESSED");
            std::string dict_file_name = out_file_prefix + ".dict";
            std::ofstream dict(dict_file_name);
        
            for (size_type i = 0; i < this->dictionary->size(); i++)
            {
                dict.write(this->dictionary->sorted_entry_at(i).c_str(), this->dictionary->sorted_entry_at(i).size());
                dict.put(ENDOFWORD);
            }
        
            dict.put(ENDOFDICT);
            
            vcfbwt::DiskWrites::update(dict.tellp()); // Disk Stats
            dict.close();
        }
        
        if (tags & COMPRESSED)
        {
            spdlog::info("Main parser: writing dictionary to disk COMPRESSED");
            std::ofstream dicz(out_file_prefix + ".dicz");
            std::ofstream lengths(out_file_prefix + ".dicz.len");
    
            for (size_type i = 0; i < this->dictionary->size(); i++)
            {
                std::size_t shift = 0;
                if (i != 0) { shift = this->w; }
                dicz.write(this->dictionary->sorted_entry_at(i).c_str() + shift,
                           this->dictionary->sorted_entry_at(i).size() - shift);
                int32_t len = this->dictionary->sorted_entry_at(i).size() - shift;
                lengths.write((char*) &len, sizeof(int32_t));
            }
    
            vcfbwt::DiskWrites::update(dicz.tellp()); // Disk Stats
            dicz.close();
    
            vcfbwt::DiskWrites::update(lengths.tellp()); // Disk Stats
            lengths.close();
        }
    
        // Outoput Occurrencies
        if(this->params.compute_occurrences)
        {
            spdlog::info("Main parser: writing occurrences to file");
            std::string occ_file_name = out_file_prefix + ".occ";
            std::ofstream occ(occ_file_name, std::ios::out | std::ios::binary);
            occ.write((char*)&occurrences[0], occurrences.size() * sizeof(size_type));
            
            vcfbwt::DiskWrites::update(occ.tellp()); // Disk Stats
            occ.close();
        }
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

//------------------------------------------------------------------------------

const std::string&
vcfbwt::pfp::AuPair::d_prime::at(std::size_t i) const
{
    if (i < d_prime_vector.size()) { return d_prime_vector.at(i); }
    else { return d_prime_map.at(i); }
}


int
vcfbwt::pfp::AuPair::cost_of_removing_trigger_string(const string_view& ts)
{
    auto& table_entry_pairs = T_table.at(ts);

    // compute the cost of removing this trigger string -----------
    int cost_of_removing_from_D = 0, cost_of_removing_from_P = 0, cost_of_removing_tot = 0;

    bool starts_and_ends_with_same_ts = false;

    // removing from P
    for (auto& pair : table_entry_pairs)
    {
        cost_of_removing_from_P -= (pair.second * sizeof(size_type));
    }

    // removing from D
    std::set<size_type> pair_seconds, pair_firsts;
    std::set<std::pair<hash_type, hash_type>> pairs;
    for (auto& pair : table_entry_pairs)
    {
        pair_firsts.insert(pair.first.first);
        pair_seconds.insert(pair.first.second);
        pairs.insert(pair.first);

        std::string_view f_ts(&(D_prime.at(pair.first.first )[0]), window_length);
        std::string_view l_ts(&(D_prime.at(pair.first.second)[D_prime.at(pair.first.second).size() - window_length]), window_length);
        if (f_ts == ts or l_ts == ts) { starts_and_ends_with_same_ts = true; }
    }

    for (auto& pair : pairs)
    {
        cost_of_removing_from_D += (this->D_prime.at(pair.first).size()
                                    + this->D_prime.at(pair.second).size() - window_length);
    }
    for (auto& p : pair_firsts)
    {
        cost_of_removing_from_D -= this->D_prime.at(p).size();
    }
    for (auto& p : pair_seconds)
    {
        cost_of_removing_from_D -= this->D_prime.at(p).size();
    }

    // flip sign of the cost so that can be compared with the threshold
    cost_of_removing_tot = cost_of_removing_from_D + cost_of_removing_from_P;
    cost_of_removing_tot = cost_of_removing_tot * -1;

    if (starts_and_ends_with_same_ts) { return 0; }
    else { return cost_of_removing_tot; }
}

void
vcfbwt::pfp::AuPair::init()
{
    spdlog::info("Reading Dictionary");
    vcfbwt::pfp::Parser::read_dictionary(in_prefix + ".dict", this->D_prime.d_prime_vector);
    this->curr_id = this->D_prime.d_prime_vector.size() + 1;
    
    spdlog::info("Reading Parse");
    mio::mmap_source in_parse(in_prefix + ".parse");
    for (std::size_t i = 0; i < in_parse.size() - sizeof(size_type); i += sizeof(size_type))
    {
        size_type curr_element = 0, next_element = 0;
        std::memcpy(&curr_element, &(in_parse[i]), sizeof(size_type));
        std::memcpy(&next_element, &(in_parse[i + sizeof(size_type)]), sizeof(size_type));
        
        const std::string& phrase_1 = D_prime.at(curr_element - 1);

        // update T
        std::string_view trigger_string(&(phrase_1[phrase_1.size() - window_length]), window_length);
        this->T_table[trigger_string][std::make_pair(curr_element - 1, next_element - 1)] += 1;

        if (i == (in_parse.size() - (2*(sizeof(size_type)))))
        {
            const std::string& phrase_2 = D_prime.at(next_element - 1);
            std::string_view trigger_string_last(&(phrase_2[phrase_2.size() - window_length]), window_length);
            auto& default_adding = this->T_table[trigger_string_last];
        }
    }

    spdlog::info("Initializing priority queue");
    this->priority_queue.init(T_table.size());
    int last_inserted = 0;
    for (const auto& table_entry : T_table)
    {
        this->trigger_string_pq_ids.insert(std::pair(table_entry.first, last_inserted));
        this->trigger_string_pq_ids_inv.insert(std::pair(last_inserted, table_entry.first));
        priority_queue.push(last_inserted, cost_of_removing_trigger_string(table_entry.first));

        last_inserted++;
    }
}

int
vcfbwt::pfp::AuPair::compress(std::set<std::string_view>& removed_trigger_strings, int threshold)
{
    spdlog::info("Start compressing with threshold {}", threshold);
    int bytes_removed = 0;
    
    spdlog::info("Initial TS count: {}", T_table.size());
    if (T_table.size() <= 1) { return 0; }

    // itearate over priority queue
    std::pair<int, int> max_cost_trigger_string = priority_queue.get_max();
    while (max_cost_trigger_string.first > threshold)
    {
        std::string_view& current_trigger_string = this->trigger_string_pq_ids_inv.at(max_cost_trigger_string.second);
        // remove trigger string if cost over threshold
        spdlog::debug("{}\tcost:\t{}\tbytes removed:\t{}", current_trigger_string, max_cost_trigger_string.first, bytes_removed);

        bytes_removed += max_cost_trigger_string.first;
        removed_trigger_strings.insert(current_trigger_string);

        std::set<std::string_view> to_update_cost;
        for (auto& pair_with_freq : T_table.at(current_trigger_string))
        {
            size_type first = pair_with_freq.first.first;
            size_type second = pair_with_freq.first.second;

            // New phrase
            std::string merged_phrase = D_prime.at(first) + D_prime.at(second).substr(window_length);
            size_type merged_phrase_id = curr_id++;
            
            D_prime.d_prime_map.insert(std::pair(merged_phrase_id, merged_phrase));

            // update entry of T where first appear as second or second appear as a first
            std::string_view first_ts(&(D_prime.at(first)[0]), window_length);
            std::string_view second_ts(&(D_prime.at(second)[D_prime.at(second).size() - window_length]) , window_length);
            
            std::vector<std::pair<std::pair<size_type, size_type>, size_type>> to_remove, to_add;
            if (T_table.contains(first_ts))
            {
                for (auto& pair : T_table.at(first_ts))
                {
                    if (pair.first.second == first)
                    {
                        to_remove.push_back(pair);
                        to_add.push_back(std::pair(std::pair(pair.first.first, merged_phrase_id), pair.second));
                    }
                }
            }

            if (T_table.contains(second_ts))
            {
                for (auto& pair : T_table.at(second_ts))
                {
                    if (pair.first.first == second)
                    {
                        to_remove.push_back(pair);
                        to_add.push_back(std::pair(std::pair(merged_phrase_id, pair.first.second), pair.second));
                    }
                }
            }

            // Perform updates
            for (auto& pair : to_remove)
            {
                std::string_view ts(&(D_prime.at(pair.first.second)[0]), window_length);
                T_table[ts].erase(pair.first);
            }

            for (auto& pair : to_add)
            {
                std::string_view ts(&(D_prime.at(pair.first.second)[0]), window_length);
                T_table[ts].insert(pair);
                
                to_update_cost.insert(ts);
            }
        }

        // Update Priority queue
        for (auto& ts : to_update_cost)
        {
            int ts_index = this->trigger_string_pq_ids.at(ts);
            this->priority_queue.push(ts_index, cost_of_removing_trigger_string(ts));
        }

        this->priority_queue.push(max_cost_trigger_string.second, 0);

        // Keep iterating
        max_cost_trigger_string = priority_queue.get_max();
    }

    return bytes_removed;
}

//------------------------------------------------------------------------------