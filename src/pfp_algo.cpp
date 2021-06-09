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
            hash_type hash = this->dictionary->check_and_add(phrase);
        
            out_file.write((char*) (&hash), sizeof(hash_type)); this->parse_size += 1;
    
            //if (phrase[0] != DOLLAR_PRIME)
            //{
            spdlog::debug("------------------------------------------------------------");
            spdlog::debug("Parsed phrase [{}] {}", phrase.size(), phrase);
            spdlog::debug("------------------------------------------------------------");
            //}
            
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

void
vcfbwt::pfp::AuPair::init(std::vector<std::string>& dictionary, std::vector<size_type>& parse)
{
    spdlog::info("Initializing AuPair structures");
    for (std::size_t i = 0; i < parse.size() - 1; i++)
    {
        std::string& phrase_1 = dictionary[parse[i] - 1];
        hash_type phrase_1_hash = string_hash(phrase_1.c_str(), phrase_1.size());

        std::string& phrase_2 = dictionary[parse[i + 1] - 1];
        hash_type phrase_2_hash = string_hash(phrase_2.c_str(), phrase_2.size());

        // update D'
        this->d_prime.insert(std::pair(phrase_1_hash, phrase_1));
        if (i == parse.size() - 2) { this->d_prime.insert(std::pair(phrase_2_hash, phrase_2)); }

        // update P'
        this->p_prime.emplace_back(phrase_1_hash);
        auto phrase_1_iterator = std::prev(p_prime.end());
        if (i == parse.size() - 2) { this->p_prime.emplace_back(phrase_2_hash); }

        // update T
        std::string trigger_string = phrase_1.substr(phrase_1.size() - window_length, window_length);
        this->T_table[trigger_string][phrase_1_hash].push_back(phrase_1_iterator);
    }
}

vcfbwt::size_type
vcfbwt::pfp::AuPair::compress(int threshold)
{
    spdlog::info("Start compressing with threshold {}", threshold);
    int bytes_removed = 0;
    
    std::set<std::string> ts_to_be_erased;
    std::set<hash_type> phrases_to_be_erased;

    spdlog::info("Initial TS count: {}", T_table.size());
    if (T_table.size() <= 1) { return 0; }


    // itearate over T
    for (auto& table_entry : this->T_table)
    {
        // compute the cost of removing this trigger string -----------
        int cost_of_removing_from_D = 0, cost_of_removing_from_P = 0, cost_of_removing_tot = 0;

        bool starts_and_ends_with_same_ts = false;

        // removing from P
        for (auto& pair_first_map : table_entry.second)
        {
            cost_of_removing_from_P -= pair_first_map.second.size() * sizeof(size_type);
        }

        // removing from D
        std::set<hash_type> pair_seconds, pair_firsts;
        std::set<std::pair<hash_type, hash_type>> pairs;
        for (auto& pair_first_map : table_entry.second)
        {
            for (auto& pair_first_ptr : pair_first_map.second)
            {
                auto pair_second = std::next(pair_first_ptr);
                pair_firsts.insert(*pair_first_ptr);
                pair_seconds.insert(*pair_second);
                pairs.insert(std::pair(*pair_first_ptr, *pair_second));
            }
        }

        for (auto& pair : pairs) { cost_of_removing_from_D += (this->d_prime.at(pair.first).size() + this->d_prime.at(pair.second).size() - window_length); }
        for (auto& p : pair_firsts)
        {
            cost_of_removing_from_D -= this->d_prime.at(p).size();
            if (d_prime.at(p).substr(0, window_length) == table_entry.first) { starts_and_ends_with_same_ts = true; }
        }
        for (auto& p : pair_seconds)
        {
            cost_of_removing_from_D -= this->d_prime.at(p).size();
            if (d_prime.at(p).substr(d_prime.at(p).size() - window_length, window_length) == table_entry.first) { starts_and_ends_with_same_ts = true; }
        }


        // flip sign of the cost so that can be compared with the threshold
        cost_of_removing_tot = cost_of_removing_from_D + cost_of_removing_from_P;
        cost_of_removing_tot = cost_of_removing_tot * -1;
        spdlog::info("{}\tcost:\t{}\tbytes removed:\t{}", table_entry.first, cost_of_removing_tot, bytes_removed);

        // remove trigger string if cost over threshold
        if (cost_of_removing_tot >= threshold and not starts_and_ends_with_same_ts)
        {
            bytes_removed += cost_of_removing_tot;
            ts_to_be_erased.insert(table_entry.first);

            for (auto& pair_tuple : table_entry.second)
            {
                for (auto& pair_first : pair_tuple.second)
                {
                    auto pair_second = std::next(pair_first);
                    hash_type pi1 = *pair_first, pi2 = *pair_second;

                    phrases_to_be_erased.insert(pi1);
                    phrases_to_be_erased.insert(pi2);

                    std::string& s1 = d_prime.at(pi1); std::string& s2 = d_prime.at(pi2);

                    std::string merged_phrase = s1 + s2.substr(window_length);
                    hash_type merged_phrase_hash = string_hash(merged_phrase.c_str(), merged_phrase.size());
                    if (d_prime.find(merged_phrase_hash) == d_prime.end())
                    {
                        d_prime.insert(std::pair(merged_phrase_hash, merged_phrase));
                    }

                    std::string pair_second_ending_ts = s2.substr(s2.size() - window_length);
                    if (std::next(pair_second) != this->p_prime.end())
                    {
                        T_table.at(pair_second_ending_ts)[merged_phrase_hash].push_back(pair_first);

                        if (pair_second_ending_ts != table_entry.first)
                        {
                            T_table.at(pair_second_ending_ts).erase(*pair_second);
                        }
                        else // if starts and ends with same trigger string delete only if
                        {
                            std::cout << "Should not be here" << std::endl;
                        }
                    }

                    hash_type &pair_first_ptr = *pair_first;
                    pair_first_ptr = merged_phrase_hash; // substitute first element
                    this->p_prime.erase(pair_second); // erase second element
                }
            }
        }
    }
     
     // remove processed trigger strings from T_table
     for (auto& ts : ts_to_be_erased) { T_table.erase(ts); }
     for (auto& p : phrases_to_be_erased) { d_prime.erase(p); }

    return bytes_removed;
}

void
vcfbwt::pfp::AuPair::close()
{
    if (this->closed) { return; }
    this->closed = true;
    
    // sort dictionary
    std::vector<std::pair<std::reference_wrapper<std::string>, hash_type>> sorted_phrases;
    for (auto& d_pair : d_prime) { sorted_phrases.push_back(std::pair(std::ref(d_pair.second), d_pair.first)); }
    
    std::sort(sorted_phrases.begin(), sorted_phrases.end(), ref_smaller);
    
    std::unordered_map<hash_type, size_type> hash_to_rank;
    for (std::size_t i = 0; i < sorted_phrases.size(); i++) { hash_to_rank[sorted_phrases[i].second] = i + 1; }
    
    spdlog::info("AuPair: writing dictionary to disk NOT COMPRESSED");
    std::string dict_file_name = this->out_prefix + ".dict";
    std::ofstream dict(dict_file_name);
    
    for (size_type i = 0; i < sorted_phrases.size(); i++)
    {
        dict.write(sorted_phrases[i].first.get().c_str(), sorted_phrases[i].first.get().size());
        dict.put(ENDOFWORD);
    }
    dict.put(ENDOFDICT);
    vcfbwt::DiskWrites::update(dict.tellp()); // Disk Stats
    dict.close();
    
    spdlog::info("AuPair: writing parse to disk");
    std::string parse_file_name = this->out_prefix + ".parse";
    std::ofstream parse(parse_file_name);
    
    for (auto& parse_element : p_prime)
    {
        if (hash_to_rank.find(parse_element) == hash_to_rank.end())
        {
            std::cout << "Error " << parse_element << std::endl;
        }
        else
        {
            parse.write((char*) &(hash_to_rank.at(parse_element)), sizeof(size_type));
        }
    }
    
    vcfbwt::DiskWrites::update(parse.tellp()); // Disk Stats
    parse.close();
}

std::string
vcfbwt::pfp::AuPair::_TESTING_unparse()
{
    std::string out;
    for (auto& p : p_prime)
    {
        out += d_prime[p].substr(0, d_prime[p].size() - window_length);
    }
    out += d_prime[p_prime.back()].substr(d_prime[p_prime.back()].size() - window_length);
    return out;
}

//------------------------------------------------------------------------------

