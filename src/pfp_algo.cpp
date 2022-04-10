//
//  pfp_algo.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#include <pfp_algo.hpp>

//------------------------------------------------------------------------------

void
vcfbwt::pfp::ReferenceParse::init(const std::string& reference)
{
    if (params.ignore_ts_file.size() > 0)
    {
        spdlog::info("Not Impelmented!!!");
    }
    
    std::vector<vcfbwt::char_type> phrase;
    spdlog::info("Parsing reference");
    
    // Karp Robin Hash Function for sliding window
    Mersenne_KarpRabinHash kr_hash(this->params.w);
    
    // Reference as first sample, just one dollar to be compatible with Giovanni's pscan.cpp
    phrase.emplace_back(DOLLAR);
    
    for (std::size_t ref_it = 0; ref_it < reference.size(); ref_it++)
    {
        char c = reference[ref_it];
        
        phrase.push_back(c);
        if (phrase.size() == params.w) { kr_hash.initialize(std::string_view((char*)phrase.data(), params.w)); }
        else if (phrase.size() > params.w) { kr_hash.update(phrase[phrase.size() - params.w - 1], phrase[phrase.size() - 1]); }
        
        if ((phrase.size() > this->params.w) and ((kr_hash.get_hash() % this->params.p) == 0))
        {
            hash_type hash = this->dictionary.check_and_add(phrase);
            
            this->parse.push_back(hash);
            this->trigger_strings_position.push_back(ref_it - this->params.w + 1);
    
            phrase.erase(phrase.begin(), phrase.end() - this->params.w); // Keep the last w chars
            
            kr_hash.reset(); kr_hash.initialize(std::string_view((char*) phrase.data(), params.w));
        }
    }
    
    // Last phrase
    if (phrase.size() > this->params.w)
    {
        // Append w-1 dollar prime and 1 dollar sequence at the end, reference as the first sample
        for (size_type j = 0; j < this->params.w - 1; j++) { phrase.emplace_back(DOLLAR_PRIME); }
        phrase.emplace_back(DOLLAR_SEQUENCE);
    
        hash_type hash = this->dictionary.check_and_add(phrase);
    
        this->parse.push_back(hash);
        this->trigger_strings_position.push_back(reference.size() - 1);
    }
    else { spdlog::error("The reference doesn't have w dollar prime at the end!"); std::exit(EXIT_FAILURE); }
}

//------------------------------------------------------------------------------

void
vcfbwt::pfp::ParserVCF::init(const Params& params, const std::string& prefix, ReferenceParse& rp, std::size_t t)
{
    this->w = params.w; this->out_file_prefix = prefix; this->p = params.p; this->tags = t; this->parse_size = 0;
    if (not ((tags & MAIN) or (tags & WORKER))) { spdlog::error("A parser must be either the main parser or a worker"); std::exit(EXIT_FAILURE); }
    
    if (tags & MAIN) {this->out_file_name = out_file_prefix + EXT::PARSE; }
    if ((tags & MAIN) and params.compress_dictionary) { tags = tags | COMPRESSED; }
    this->tmp_out_file_name = TempFile::getName("parse");
    this->out_file.open(tmp_out_file_name, std::ios::binary);
    this->reference_parse = &rp;
    this->dictionary = &this->reference_parse->dictionary;
    
    this->params = params;
}

void
vcfbwt::pfp::ParserVCF::operator()(const vcfbwt::Sample& sample)
{
    Sample::iterator sample_iterator(sample, this->working_genotype);
    this->samples_processed.push_back(sample.id());
    
    std::vector<vcfbwt::char_type> phrase;
    
    // Karp Robin Hash Function for sliding window
    Mersenne_KarpRabinHash kr_hash(this->params.w);
    
    // Every sample starts with w-1 dollar prime and one dollar seq
    for (size_type j = 0; j < this->params.w - 1; j++) { phrase.emplace_back(DOLLAR_PRIME); }
    phrase.emplace_back(DOLLAR_SEQUENCE);
    kr_hash.initialize(std::string_view((char*)phrase.data(), params.w));
    
    // Shorthands
    std::vector<std::size_t>& tsp = reference_parse->trigger_strings_position;
    
    std::size_t start_window = 0, end_window = 0;
    while (not sample_iterator.end())
    {
        // Compute where we are on the reference
        std::size_t pos_on_reference = sample_iterator.get_ref_it();
        
        if ( not ((sample_iterator.get_var_it() > 0) and (sample_iterator.prev_variation_end() > (pos_on_reference - (2 * this->w)))))
        {
            // Set start postion to the position in the reference parse after the last computed phrase
            if (params.use_acceleration and ((phrase.size() == this->w) and ((pos_on_reference != 0) and (phrase[0] != DOLLAR_PRIME))))
            {
                start_window = end_window;
                while ((tsp[start_window] + this->w) <= pos_on_reference and (start_window < tsp.size() - 2))
                { start_window++; }
    
                // Iterate over the parse up to the next variation
                while (tsp[end_window + 1] < (sample_iterator.next_variation() - (this->w + 1))) { end_window++; }
                
                // If the window is not empty
                if ((start_window < end_window - 1) and (tsp[end_window] > pos_on_reference))
                {
                    spdlog::debug("------------------------------------------------------------");
                    spdlog::debug("from {}", sample.get_reference().substr(tsp[start_window - 1], this->w));
                    spdlog::debug("copied from {} to {}", tsp[start_window], tsp[end_window] + this->w);
                    spdlog::debug("next variation: {}", sample_iterator.next_variation());
                    spdlog::debug("skipped phrases: {}", end_window - start_window);
                    
                    // copy from parse[start_window : end_window]
                    out_file.write((char*) &(this->reference_parse->parse[start_window]), sizeof(hash_type) * (end_window - start_window + 1));
                    this->parse_size += end_window - start_window + 1;
            
                    // move iterators and re initialize phrase
                    sample_iterator.go_to(tsp[end_window]);
                    phrase.clear();
                    for (std::size_t i = 0; i < this->w; i++) { ++sample_iterator; phrase.push_back(*sample_iterator);}
                    
                    kr_hash.reset(); kr_hash.initialize(string_view((char*) phrase.data(), params.w));
                    
                    ++sample_iterator;
                    spdlog::debug("New phrase [{}]: {}", phrase.size(), std::string((char*) phrase.data(), phrase.size()));
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
    
            if (phrase[0] != DOLLAR_PRIME)
            {
                spdlog::debug("------------------------------------------------------------");
                spdlog::debug("Parsed phrase [{}] {}", phrase.size(), std::string((char*) phrase.data(), phrase.size()));
                spdlog::debug("------------------------------------------------------------");
            }
            
            phrase.erase(phrase.begin(), phrase.end() - this->w); // Keep the last w chars
    
            kr_hash.reset(); kr_hash.initialize(std::string_view((char*) phrase.data(), params.w));
        }
    }
    
    // Last phrase
    if (phrase.size() > this->w)
    {
        // Append w dollar prime at the end of each sample, also w DOLLAR if it's the last sample
        phrase.insert(phrase.end(), params.w - 1, DOLLAR_PRIME);
        if (sample.last(this->working_genotype)) { phrase.insert(phrase.end(), params.w, DOLLAR); }
        else { phrase.emplace_back(DOLLAR_SEQUENCE); }

        hash_type hash = this->dictionary->check_and_add(phrase);
        
        out_file.write((char*) (&hash), sizeof(hash_type));   this->parse_size += 1;
    }
    else { spdlog::error("A sample doesn't have w dollar prime at the end!"); std::exit(EXIT_FAILURE); }
}

void
vcfbwt::pfp::ParserVCF::close()
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
        
        spdlog::info("Main parser: Replacing hash values with ranks in MAIN, WORKERS and reference, wirting .last ans .sai");
        
        std::string last_file_name = out_file_prefix + EXT::LAST;
        std::ofstream last_file(last_file_name);
    
        std::string sai_file_name = out_file_prefix + EXT::SAI;
        std::ofstream sai_file(sai_file_name);
        
        std::size_t pos_for_sai = 0;
    
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
    
                const std::vector<vcfbwt::char_type>& dict_string = this->dictionary->sorted_entry_at(rank - 1);
                last_file.put(dict_string[(dict_string.size() - this->params.w) - 1]);
    
                if (pos_for_sai == 0) { pos_for_sai = dict_string.size() - 1; } // -1 is for the initial $ of the first word
                else { pos_for_sai += dict_string.size() - this->params.w; }
                sai_file.write((char*) &pos_for_sai, IBYTES);
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
    
                const std::vector<vcfbwt::char_type>& dict_string = this->dictionary->sorted_entry_at(rank - 1);
                last_file.put(dict_string[(dict_string.size() - this->params.w) - 1]);
    
                if (pos_for_sai == 0) { pos_for_sai = dict_string.size() - 1; } // -1 is for the initial $ of the first word
                else { pos_for_sai += dict_string.size() - this->params.w; }
                sai_file.write((char*) &pos_for_sai, IBYTES);
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
    
                    const std::vector<vcfbwt::char_type>& dict_string = this->dictionary->sorted_entry_at(rank - 1);
                    last_file.put(dict_string[(dict_string.size() - this->params.w) - 1]);
    
                    if (pos_for_sai == 0) { pos_for_sai = dict_string.size() - 1; } // -1 is for the initial $ of the first word
                    else { pos_for_sai += dict_string.size() - this->params.w; }
                    sai_file.write((char*) &pos_for_sai, IBYTES);
                }
                rw_mmap.unmap();
                truncate_file(worker.get().tmp_out_file_name, worker.get().parse_size * sizeof(size_type));
            }
        }
    
        vcfbwt::DiskWrites::update(last_file.tellp());
        last_file.close();
    
        vcfbwt::DiskWrites::update(sai_file.tellp());
        sai_file.close();
        
        // Merging files together
        spdlog::info("Main parser: concatenating parsings from workers and reference, reference as first");
        std::ofstream merged(out_file_name, std::ios_base::binary);
        
        // Reference
        std::size_t out_parse_size = 0;
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
            std::string dict_file_name = out_file_prefix + EXT::DICT;
            std::ofstream dict(dict_file_name);
        
            for (size_type i = 0; i < this->dictionary->size(); i++)
            {
                dict.write((char*) this->dictionary->sorted_entry_at(i).data(), this->dictionary->sorted_entry_at(i).size());
                dict.put(ENDOFWORD);
            }
        
            dict.put(ENDOFDICT);
            
            vcfbwt::DiskWrites::update(dict.tellp()); // Disk Stats
            dict.close();
        }
        
        if (tags & COMPRESSED)
        {
            spdlog::info("Main parser: writing dictionary to disk COMPRESSED");
            std::ofstream dicz(out_file_prefix + EXT::DICT_COMPRESSED);
            std::ofstream lengths(out_file_prefix + EXT::DICT_COMPRESSED_LENGTHS);

            for (size_type i = 0; i < this->dictionary->size(); i++)
            {
                std::size_t shift = 1; // skip dollar on first phrase
                if (i != 0) { shift = this->w; }
                dicz.write((char*) this->dictionary->sorted_entry_at(i).data() + shift,
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
            std::string occ_file_name = out_file_prefix + EXT::OCC;
            std::ofstream occ(occ_file_name, std::ios::out | std::ios::binary);
            occ.write((char*)&occurrences[0], occurrences.size() * sizeof(size_type));
            
            vcfbwt::DiskWrites::update(occ.tellp()); // Disk Stats
            occ.close();
        }

        spdlog::info("Main parser: closed");
    }
}

//------------------------------------------------------------------------------

void
vcfbwt::pfp::ParserFasta::init(const Params& params, const std::string& prefix)
{
    this->params = params;
    this->w = this->params.w;
    this->p = this->params.p;
    this->parse_size = 0;
    this->out_file_prefix = prefix;
    this->out_file_name = prefix + EXT::PARSE;
    this->out_file.open(out_file_name);
}

void
vcfbwt::pfp::ParserFasta::operator()()
{
    // Open input file with kseq
    gzFile fp; kseq_t *record;
    fp = gzopen(this->in_file_path.c_str(), "r");
    if (fp == 0)
    {
        spdlog::error("Failed to open input file {}", in_file_path);
        exit(EXIT_FAILURE);
    }
    
    std::vector<vcfbwt::char_type> phrase;
    spdlog::info("Parsing sequence");
    
    // Karp Robin Hash Function for sliding window
    Mersenne_KarpRabinHash kr_hash(this->params.w);
    
    // First sequence start with one dollar
    phrase.emplace_back(DOLLAR);
    
    record = kseq_init(fp);
    while(kseq_read(record) >= 0)
    {
        std::string sequence_name("<error reading sequence name>"), sequence_comment;
        if (record->name.s != NULL) { sequence_name = record->name.s; }
        if (record->comment.s != NULL) { sequence_comment = record->comment.s; }
        this->sequences_processed.push_back(sequence_name + " " + sequence_comment);
        spdlog::debug("Parsed:\t{}", sequence_name + " " + sequence_comment);
        
        // Previous last phrase
        if (phrase[0] != DOLLAR and phrase.size() >= this->params.w)
        {
            // Append w-1 dollar prime, and one dollar seq at the end of each sequence
            for (size_type j = 0; j < this->params.w - 1; j++) { phrase.emplace_back(DOLLAR_PRIME); }
            phrase.emplace_back(DOLLAR_SEQUENCE);
    
            hash_type hash = this->dictionary.check_and_add(phrase);
     
            out_file.write((char*) (&hash), sizeof(hash_type)); this->parse_size += 1;
    
            // Reset phrase
            phrase.clear();
            for (size_type j = 0; j < this->params.w - 1; j++) { phrase.emplace_back(DOLLAR_PRIME); }
            phrase.emplace_back(DOLLAR_SEQUENCE);
            kr_hash.reset(); kr_hash.initialize(std::string_view((char*) phrase.data(), params.w));
        }
        
        for (std::size_t seq_it = 0; seq_it < record->seq.l; seq_it++)
        {
            char c = record->seq.s[seq_it];
        
            phrase.push_back(c);
            if (phrase.size() == params.w) { kr_hash.initialize(std::string_view((char*) phrase.data(), params.w)); }
            else if (phrase.size() > params.w) { kr_hash.update(phrase[phrase.size() - params.w - 1], phrase[phrase.size() - 1]); }
        
            if ((phrase.size() > this->params.w) and ((kr_hash.get_hash() % this->params.p) == 0))
            {
                hash_type hash = this->dictionary.check_and_add(phrase);
    
                out_file.write((char*) (&hash), sizeof(hash_type)); this->parse_size += 1;
                
                phrase.erase(phrase.begin(), phrase.end() - this->params.w); // Keep the last w chars
                
                kr_hash.reset(); kr_hash.initialize(std::string_view((char*) phrase.data(), params.w));
            }
        }
    }
    
    // Last phrase
    if (phrase.size() >= this->params.w)
    {
        // Append w-1 dollar prime and a dollar seq
        phrase.insert(phrase.end(), this->params.w - 1, DOLLAR_PRIME);
        phrase.emplace_back(DOLLAR_SEQUENCE);

        // Append w dollar at the end
        phrase.insert(phrase.end(), this->params.w, DOLLAR);
        
        hash_type hash = this->dictionary.check_and_add(phrase);
        
        out_file.write((char*) (&hash), sizeof(hash_type)); this->parse_size += 1;
    }
    else { spdlog::error("Missing w DOLLAR at the end!"); std::exit(EXIT_FAILURE); }
    
    kseq_destroy(record);
    gzclose(fp);
}


void
vcfbwt::pfp::ParserFasta::close()
{
    if (closed) return; closed = true;
    
    vcfbwt::DiskWrites::update(out_file.tellp()); // Disk Stats
    this->out_file.close();
    
    // Occurrences
    std::vector<size_type> occurrences(this->dictionary.size(), 0);
    
    spdlog::info("Main parser: Replacing hash values with ranks, writing .last and .sai");
    
    std::string last_file_name = out_file_prefix + EXT::LAST;
    std::ofstream last_file(last_file_name);
    
    std::string sai_file_name = out_file_prefix + EXT::SAI;
    std::ofstream sai_file(sai_file_name);
    
    std::size_t pos_for_sai = 0;
    
    // mmap file and substitute
    if (this->parse_size != 0)
    {
        std::error_code error;
        mio::mmap_sink rw_mmap = mio::make_mmap_sink(out_file_name, 0, mio::map_entire_file, error);
        if (error) { spdlog::error(error.message()); std::exit(EXIT_FAILURE); }
        
        for (size_type i = 0; i < (rw_mmap.size() / sizeof(hash_type)); i++)
        {
            hash_type hash;
            std::memcpy(&hash, rw_mmap.data() + (i * sizeof(hash_type)), sizeof(hash_type));
            size_type rank = this->dictionary.hash_to_rank(hash);
            std::memcpy(rw_mmap.data() + (i * sizeof(size_type)), &rank, sizeof(size_type));
            occurrences[rank - 1] += 1;
    
            const std::vector<vcfbwt::char_type>& dict_string = this->dictionary.sorted_entry_at(rank - 1);
            last_file.put(dict_string[(dict_string.size() - this->params.w) - 1]);
    
            if (pos_for_sai == 0) { pos_for_sai = dict_string.size() - 1; } // -1 is for the initial $ of the first word
            else { pos_for_sai += dict_string.size() - this->params.w; }
            sai_file.write((char*) &pos_for_sai, IBYTES);
        }
        rw_mmap.unmap();
        truncate_file(out_file_name, this->parse_size * sizeof(size_type));
    }
    
    vcfbwt::DiskWrites::update(last_file.tellp());
    last_file.close();
    
    vcfbwt::DiskWrites::update(sai_file.tellp());
    sai_file.close();
    
    // Print dicitionary on disk
    spdlog::info("Main parser: writing dictionary on disk NOT COMPRESSED");
    std::string dict_file_name = out_file_prefix + EXT::DICT;
    std::ofstream dict(dict_file_name);
    
    for (size_type i = 0; i < this->dictionary.size(); i++)
    {
        dict.write((char*) this->dictionary.sorted_entry_at(i).data(), this->dictionary.sorted_entry_at(i).size());
        dict.put(ENDOFWORD);
    }
    
    dict.put(ENDOFDICT);
    
    vcfbwt::DiskWrites::update(dict.tellp()); // Disk Stats
    dict.close();
    
    if (this->params.compress_dictionary)
    {
        spdlog::info("Main parser: writing dictionary on disk COMPRESSED");
        std::ofstream dicz(out_file_prefix + EXT::DICT_COMPRESSED);
        std::ofstream lengths(out_file_prefix + EXT::DICT_COMPRESSED_LENGTHS);

        for (size_type i = 0; i < this->dictionary.size(); i++)
        {
            std::size_t shift = 1; // skip dollar on first phrase
            if (i != 0) { shift = this->w; }
            dicz.write((char*) this->dictionary.sorted_entry_at(i).data() + shift,
                       this->dictionary.sorted_entry_at(i).size() - shift);
            int32_t len = this->dictionary.sorted_entry_at(i).size() - shift;
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
        std::string occ_file_name = out_file_prefix + EXT::OCC;
        std::ofstream occ(occ_file_name, std::ios::out | std::ios::binary);
        occ.write((char*)&occurrences[0], occurrences.size() * sizeof(size_type));
        
        vcfbwt::DiskWrites::update(occ.tellp()); // Disk Stats
        occ.close();
    }

    spdlog::info("Main parser: closed");
}

//------------------------------------------------------------------------------

void
vcfbwt::pfp::ParserText::init(const Params& params, const std::string& prefix)
{
    this->params = params;
    this->w = this->params.w;
    this->p = this->params.p;
    this->parse_size = 0;
    this->out_file_prefix = prefix;
    this->out_file_name = prefix + EXT::PARSE;
    this->out_file.open(out_file_name);
}

void
vcfbwt::pfp::ParserText::operator()()
{
    // Open input file with kseq
    gzFile fp;
    fp = gzopen(this->in_file_path.c_str(), "r");
    if (fp == nullptr)
    {
        spdlog::error("Failed to open input file {}", in_file_path);
        exit(EXIT_FAILURE);
    }
    
    std::vector<vcfbwt::char_type> phrase;
    spdlog::info("Parsing {}", in_file_path);
    
    // Karp Robin Hash Function for sliding window
    Mersenne_KarpRabinHash kr_hash(this->params.w);
    
    // First sequence start with one dollar
    phrase.emplace_back(DOLLAR);
    
    char c;
    while(gzread(fp, &c, 1) > 0)
    {
        phrase.push_back(c);
        if (phrase.size() == params.w) { kr_hash.initialize(std::string_view((char*) phrase.data(), params.w)); }
        else if (phrase.size() > params.w) { kr_hash.update(phrase[phrase.size() - params.w - 1], phrase[phrase.size() - 1]); }
        
        if ((phrase.size() > this->params.w) and ((kr_hash.get_hash() % this->params.p) == 0))
        {
            hash_type hash = this->dictionary.check_and_add(phrase);
            
            out_file.write((char*) (&hash), sizeof(hash_type)); this->parse_size += 1;
            
            phrase.erase(phrase.begin(), phrase.end() - this->params.w); // Keep the last w chars
            
            kr_hash.reset(); kr_hash.initialize(std::string_view((char*) phrase.data(), params.w));
        }
    }
    
    // Last phrase
    if (phrase.size() >= this->params.w)
    {
        // Append w dollar at the end
       phrase.insert(phrase.end(), this->params.w, DOLLAR);
        
        hash_type hash = this->dictionary.check_and_add(phrase);
        
        out_file.write((char*) (&hash), sizeof(hash_type)); this->parse_size += 1;
    }
    else { spdlog::error("A sequence doesn't have w DOLLAR at the end!"); std::exit(EXIT_FAILURE); }
    
    gzclose(fp);
}


void
vcfbwt::pfp::ParserText::close()
{
    if (closed) return; closed = true;
    
    vcfbwt::DiskWrites::update(out_file.tellp()); // Disk Stats
    this->out_file.close();
    
    // Occurrences
    std::vector<size_type> occurrences(this->dictionary.size(), 0);
    
    spdlog::info("Main parser: Replacing hash values with ranks, writing .last and .sai");
  
    std::string last_file_name = out_file_prefix + EXT::LAST;
    std::ofstream last_file(last_file_name);
    
    std::string sai_file_name = out_file_prefix + EXT::SAI;
    std::ofstream sai_file(sai_file_name);
    
    std::size_t pos_for_sai = 0;
    
    // mmap file and substitute
    if (this->parse_size != 0)
    {
        std::error_code error;
        mio::mmap_sink rw_mmap = mio::make_mmap_sink(out_file_name, 0, mio::map_entire_file, error);
        if (error) { spdlog::error(error.message()); std::exit(EXIT_FAILURE); }
        
        for (size_type i = 0; i < (rw_mmap.size() / sizeof(hash_type)); i++)
        {
            hash_type hash;
            std::memcpy(&hash, rw_mmap.data() + (i * sizeof(hash_type)), sizeof(hash_type));
            size_type rank = this->dictionary.hash_to_rank(hash);
            std::memcpy(rw_mmap.data() + (i * sizeof(size_type)), &rank, sizeof(size_type));
            occurrences[rank - 1] += 1;
            
            const std::vector<vcfbwt::char_type>& dict_string = this->dictionary.sorted_entry_at(rank - 1);
            last_file.put(dict_string[(dict_string.size() - this->params.w) - 1]);
    
            if (pos_for_sai == 0) { pos_for_sai = dict_string.size() - 1; } // -1 is for the initial $ of the first word
            else { pos_for_sai += dict_string.size() - this->params.w; }
            sai_file.write((char*) &pos_for_sai, IBYTES);
        }
        rw_mmap.unmap();
        truncate_file(out_file_name, this->parse_size * sizeof(size_type));
    }
    
    vcfbwt::DiskWrites::update(last_file.tellp());
    last_file.close();
    
    vcfbwt::DiskWrites::update(sai_file.tellp());
    sai_file.close();
    
    // Print dicitionary on disk
    spdlog::info("Main parser: writing dictionary on disk NOT COMPRESSED");
    std::string dict_file_name = out_file_prefix + EXT::DICT;
    std::ofstream dict(dict_file_name);
    
    for (size_type i = 0; i < this->dictionary.size(); i++)
    {
        dict.write((char*) this->dictionary.sorted_entry_at(i).data(), this->dictionary.sorted_entry_at(i).size());
        dict.put(ENDOFWORD);
    }
    
    dict.put(ENDOFDICT);
    
    vcfbwt::DiskWrites::update(dict.tellp()); // Disk Stats
    dict.close();
    
    if (this->params.compress_dictionary)
    {
        spdlog::info("Main parser: writing dictionary on disk COMPRESSED");
        std::ofstream dicz(out_file_prefix + EXT::DICT_COMPRESSED);
        std::ofstream lengths(out_file_prefix + EXT::DICT_COMPRESSED_LENGTHS);
        
        for (size_type i = 0; i < this->dictionary.size(); i++)
        {
            std::size_t shift = 1; // skip dollar on first phrase
            if (i != 0) { shift = this->w; }
            dicz.write((char*) this->dictionary.sorted_entry_at(i).data() + shift,
                       this->dictionary.sorted_entry_at(i).size() - shift);
            int32_t len = this->dictionary.sorted_entry_at(i).size() - shift;
            lengths.write((char*) &len, sizeof(int32_t));
        }
        
        vcfbwt::DiskWrites::update(dicz.tellp()); // Disk Stats
        dicz.close();
        
        vcfbwt::DiskWrites::update(lengths.tellp()); // Disk Stats
        lengths.close();
    }
    
    // Output Occurrencies
    if(this->params.compute_occurrences)
    {
        spdlog::info("Main parser: writing occurrences to file");
        std::string occ_file_name = out_file_prefix + EXT::OCC;
        std::ofstream occ(occ_file_name, std::ios::out | std::ios::binary);
        occ.write((char*)&occurrences[0], occurrences.size() * sizeof(size_type));
        
        vcfbwt::DiskWrites::update(occ.tellp()); // Disk Stats
        occ.close();
    }

    spdlog::info("Main parser: closed");
}

//------------------------------------------------------------------------------

void
vcfbwt::pfp::ParserIntegers::init(const Params& params, const std::string& prefix)
{
    this->params = params;
    this->w = this->params.w;
    this->p = this->params.p;
    this->parse_size = 0;
    this->out_file_prefix = prefix;
    this->out_file_name = prefix + EXT::PARSE;
    this->out_file.open(out_file_name);
}

void
vcfbwt::pfp::ParserIntegers::operator()()
{
    // Open input file with kseq
    gzFile fp;
    fp = gzopen(this->in_file_path.c_str(), "r");
    if (fp == nullptr)
    {
        spdlog::error("Failed to open input file {}", in_file_path);
        exit(EXIT_FAILURE);
    }

    std::vector<int32_t> phrase;
    spdlog::info("Parsing {}", in_file_path);

    // Karp Robin Hash Function for sliding window
    Mersenne_KarpRabinHash4 kr_hash(this->params.w * 4);

    // First sequence start with one dollar
    phrase.emplace_back(DOLLAR);

    int32_t c;
    while(gzread(fp, &c, 4) > 0)
    {
        phrase.push_back(c + this->params.integers_shift);
        if (phrase.size() == params.w) { kr_hash.initialize(std::string_view((char*) phrase.data(), phrase.size() * sizeof(int32_t))); }
        else if (phrase.size() > params.w)
        {
            kr_hash.update((const vcfbwt::char_type*)&(phrase[phrase.size() - params.w - 1]), (const vcfbwt::char_type*) &phrase[phrase.size() - 1]);
        }

        if ((phrase.size() > this->params.w) and ((kr_hash.get_hash() % this->params.p) == 0))
        {
            hash_type hash = this->dictionary.check_and_add(phrase);

            out_file.write((char*) (&hash), sizeof(hash_type)); this->parse_size += 1;

            phrase.erase(phrase.begin(), phrase.end() - this->params.w); // Keep the last w chars

            kr_hash.reset(); kr_hash.initialize(std::string_view((char*) phrase.data(), phrase.size() * sizeof(int32_t)));
        }
    }

    // Last phrase
    if (phrase.size() >= this->params.w)
    {
        // Append w dollar at the end
        phrase.insert(phrase.end(), this->params.w, DOLLAR);

        hash_type hash = this->dictionary.check_and_add(phrase);

        out_file.write((char*) (&hash), sizeof(hash_type)); this->parse_size += 1;
    }
    else { spdlog::error("A sequence doesn't have w DOLLAR at the end!"); std::exit(EXIT_FAILURE); }

    gzclose(fp);
}


void
vcfbwt::pfp::ParserIntegers::close()
{
    if (closed) return; closed = true;

    vcfbwt::DiskWrites::update(out_file.tellp()); // Disk Stats
    this->out_file.close();

    // Occurrences
    std::vector<size_type> occurrences(this->dictionary.size(), 0);

    spdlog::info("Main parser: Replacing hash values with ranks, writing .last and .sai");

    std::string last_file_name = out_file_prefix + EXT::LAST;
    std::ofstream last_file(last_file_name);

    std::string sai_file_name = out_file_prefix + EXT::SAI;
    std::ofstream sai_file(sai_file_name);

    std::size_t pos_for_sai = 0;

    // mmap file and substitute
    if (this->parse_size != 0)
    {
        std::error_code error;
        mio::mmap_sink rw_mmap = mio::make_mmap_sink(out_file_name, 0, mio::map_entire_file, error);
        if (error) { spdlog::error(error.message()); std::exit(EXIT_FAILURE); }

        for (size_type i = 0; i < (rw_mmap.size() / sizeof(hash_type)); i++)
        {
            hash_type hash;
            std::memcpy(&hash, rw_mmap.data() + (i * sizeof(hash_type)), sizeof(hash_type));
            size_type rank = this->dictionary.hash_to_rank(hash);
            std::memcpy(rw_mmap.data() + (i * sizeof(size_type)), &rank, sizeof(size_type));
            occurrences[rank - 1] += 1;

            const std::vector<int32_t>& dict_string = this->dictionary.sorted_entry_at(rank - 1);
            last_file.write((char*) &(dict_string[(dict_string.size() - this->params.w) - 1]), sizeof(size_type));

            if (pos_for_sai == 0) { pos_for_sai = dict_string.size() - 1; } // -1 is for the initial $ of the first word
            else { pos_for_sai += dict_string.size() - this->params.w; }
            sai_file.write((char*) &pos_for_sai, IBYTES);
        }
        rw_mmap.unmap();
        truncate_file(out_file_name, this->parse_size * sizeof(size_type));
    }

    vcfbwt::DiskWrites::update(last_file.tellp());
    last_file.close();

    vcfbwt::DiskWrites::update(sai_file.tellp());
    sai_file.close();

    // Print dicitionary on disk
    spdlog::info("Main parser: writing dictionary on disk NOT COMPRESSED");
    std::string dict_file_name = out_file_prefix + EXT::DICT;
    std::ofstream dict(dict_file_name);

    for (size_type i = 0; i < this->dictionary.size(); i++)
    {
        dict.write((char*) this->dictionary.sorted_entry_at(i).data(), this->dictionary.sorted_entry_at(i).size() * sizeof(int32_t));
        int32_t eow = ENDOFWORD;
        dict.write((char*) &eow, sizeof(int32_t));
    }

    int32_t eod = ENDOFDICT;
    dict.write((char*) &eod, sizeof(int32_t));

    vcfbwt::DiskWrites::update(dict.tellp()); // Disk Stats
    dict.close();

    if (this->params.compress_dictionary)
    {
        spdlog::info("Main parser: writing dictionary on disk COMPRESSED");
        std::ofstream dicz(out_file_prefix + EXT::DICT_COMPRESSED);
        std::ofstream lengths(out_file_prefix + EXT::DICT_COMPRESSED_LENGTHS);

        for (size_type i = 0; i < this->dictionary.size(); i++)
        {
            std::size_t shift = 1; // skip dollar on first phrase
            if (i != 0) { shift = this->w; }
            dicz.write((char*) (this->dictionary.sorted_entry_at(i).data() + shift),
                       (this->dictionary.sorted_entry_at(i).size() - shift) * sizeof(int32_t));
            int32_t len = this->dictionary.sorted_entry_at(i).size() - shift;
            lengths.write((char*) &len, sizeof(int32_t));
        }

        vcfbwt::DiskWrites::update(dicz.tellp()); // Disk Stats
        dicz.close();

        vcfbwt::DiskWrites::update(lengths.tellp()); // Disk Stats
        lengths.close();
    }

    // Output Occurrencies
    if(this->params.compute_occurrences)
    {
        spdlog::info("Main parser: writing occurrences to file");
        std::string occ_file_name = out_file_prefix + EXT::OCC;
        std::ofstream occ(occ_file_name, std::ios::out | std::ios::binary);
        occ.write((char*)&occurrences[0], occurrences.size() * sizeof(size_type));

        vcfbwt::DiskWrites::update(occ.tellp()); // Disk Stats
        occ.close();
    }

    spdlog::info("Main parser: closed");
}

//------------------------------------------------------------------------------