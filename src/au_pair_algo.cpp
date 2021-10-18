//
//  au_pair_algo.cpp.c
//
//  Copyright 2021 Marco Oliva. All rights reserved.
//

#include <au_pair_algo.hpp>

//------------------------------------------------------------------------------

const std::string&
vcfbwt::pfp::AuPair::d_prime::at(std::size_t i) const
{
    if (i < d_prime_vector.size()) { return d_prime_vector.at(i); }
    else { return d_prime_map.at(i); }
}

void vcfbwt::pfp::AuPair::d_prime::remove(vcfbwt::size_type i)
{
    if (i < d_prime_vector.size()) { d_prime_vector[i] = ""; }
    else { d_prime_map.erase(i); }
}


int
vcfbwt::pfp::AuPair::cost_of_removing_trigger_string(const string_view& ts)
{
    if (ts[0] == DOLLAR or ts[0] == DOLLAR_PRIME) { return 0; } // don't want to move structural dollar signs
    if (not T_table.contains(ts)) { return 0; }
    
    auto& table_entry = T_table.at(ts);
    
    // compute the cost of removing this trigger string -----------
    int cost_of_removing_from_D = 0, cost_of_removing_from_P = 0, cost_of_removing_tot = 0;
    
    // removing from D
    std::set<hash_type> pair_seconds, pair_firsts;
    std::set<std::pair<hash_type, hash_type>> pairs;
    std::size_t erased_elements = 0;
    for (std::size_t i = 0; i < table_entry.size(); i++)
    {
        size_type* pair_first_ptr = table_entry[i];
        
        // Clean the vector to reduce memory usage, we are already iterating anyway
        if (parse.removed(pair_first_ptr)) { erased_elements += 1; continue; }
        
        size_type pair_first_v = (*pair_first_ptr) - 1;
        size_type pair_second_v = 0;
        
        if (this->parse.next(pair_first_ptr) != this->parse.end()) { pair_second_v = (*(parse.next(pair_first_ptr))) - 1; }
        else { continue; }
        
        // check ts
        std::string_view check_f(&(D_prime.at(pair_first_v)[D_prime.at(pair_first_v).size() - window_length]) , window_length);
        std::string_view check_s(&(D_prime.at(pair_second_v)[0]), window_length);
        if (check_f != check_s) { return 0; }
        
        pair_firsts.insert(pair_first_v);
        pair_seconds.insert(pair_second_v);
        pairs.insert(std::make_pair(pair_first_v, pair_second_v));
        
        std::string_view f_ts(&(D_prime.at(pair_first_v)[0]), window_length);
        std::string_view l_ts(&(D_prime.at(pair_second_v)[D_prime.at(pair_second_v).size() - window_length]), window_length);
        
        if (f_ts[0] == DOLLAR or f_ts[0] == DOLLAR_PRIME) { return 0; }
        if (l_ts[0] == DOLLAR or l_ts[0] == DOLLAR_PRIME) { return 0; }
        if (f_ts == ts or l_ts == ts) { return 0; }
    }
    
    // removing from P, no empty elements here!
    cost_of_removing_from_P -= (table_entry.size() - erased_elements) * sizeof(size_type);
    
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
    
    return cost_of_removing_tot;
}

void
vcfbwt::pfp::AuPair::init_structures()
{
    spdlog::info("Reading Dictionary");
    vcfbwt::pfp::ParserUtils::read_dictionary(in_prefix + EXT::DICT, this->D_prime.d_prime_vector);
    this->curr_id = this->D_prime.d_prime_vector.size() + 10;
    
    spdlog::info("Reading Parse");
    mio::basic_mmap_source<char> in_parse(in_prefix + EXT::PARSE);
    this->parse.init((const size_type*) in_parse.data(), in_parse.size() / sizeof(size_type));
    
    // initial size
    this->initial_size = this->parse.size() * sizeof(size_type);
    for (auto& s : this->D_prime.d_prime_vector) { this->initial_size += s.size(); }
    
    spdlog::info("Initializing data structures to compute costs");
    for (std::size_t i = 0; i < this->parse.size() - 1; i++)
    {
        const std::string& phrase_1 = D_prime.at(this->parse.at(i) - 1);
        
        // update T
        std::string_view trigger_string(&(phrase_1[phrase_1.size() - window_length]), window_length);
        this->T_table[trigger_string].push_back(&(this->parse.at(i)));
        
        if (i == 0)
        {
            std::string_view trigger_string_first(&(D_prime.at(0))[0], window_length);
            auto& default_adding = this->T_table[trigger_string_first];
        }
        if (i == (this->parse.size() - 2))
        {
            const std::string& phrase_2 = D_prime.at(this->parse.at(i + 1) - 1);
            std::string_view trigger_string_last(&(phrase_2[phrase_2.size() - window_length]), window_length);
            auto& default_adding = this->T_table[trigger_string_last];
        }
    }
}

void
vcfbwt::pfp::AuPair::init_costs()
{
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

std::size_t
vcfbwt::pfp::AuPair::remove_simple(std::set<std::string_view>& removed_trigger_strings)
{
    spdlog::info("Removing simple trigger strings");
    
    std::size_t bytes_removed = 0;
    std::set<size_type> removed_phrases;
    for (auto& ts_pair : T_table)
    {
        // Remove trigger strings with only one pair in the list
        std::string_view current_trigger_string = ts_pair.first;
    
        // checks for integrity
        if ((removed_trigger_strings.contains(current_trigger_string)) or
        (current_trigger_string[0] == DOLLAR or current_trigger_string[0] == DOLLAR_PRIME))
        { continue; }
    
        std::set<size_type> pair_elements;
    
        bool checks_failed = false;
        for (auto pair_first_ptr : T_table.at(current_trigger_string))
        {
            if (parse.removed(pair_first_ptr)) { continue; }
        
            size_type pair_first_v = (*pair_first_ptr) - 1;
            size_type pair_second_v = 0;
        
            if (this->parse.next(pair_first_ptr) != this->parse.end()) { pair_second_v = (*(parse.next(pair_first_ptr))) - 1; }
            else { continue; }
        
            std::string_view f_ts(&(D_prime.at(pair_first_v)[0]), window_length);
            std::string_view l_ts(&(D_prime.at(pair_second_v)[D_prime.at(pair_second_v).size() - window_length]), window_length);
    
            if (f_ts[0] == DOLLAR
            or (f_ts[0] == DOLLAR_PRIME
            or (l_ts[0] == DOLLAR
            or (l_ts[0] == DOLLAR_PRIME
            or (f_ts == current_trigger_string or l_ts == current_trigger_string)))))
            { checks_failed = true; break; }
        
            pair_elements.insert(pair_first_v);
            pair_elements.insert(pair_second_v);
            if (pair_elements.size() > 2) { break; }
        }
    
        if (pair_elements.size() == 2 and not checks_failed)
        {
            removed_trigger_strings.insert(current_trigger_string);
            bytes_removed += sizeof(size_type) * T_table.at(current_trigger_string).size(); // form P
            bytes_removed += window_length; // form D
    
            // remove this trigger string
            spdlog::debug("{}\tbytes removed:\t{}\tremoved ts: {}\tT size: {}",
                         current_trigger_string, bytes_removed, removed_trigger_strings.size(), T_table.size());
            
            // update parse and dict
            std::map<std::pair<size_type, size_type>, hash_type> merged_pairs;
            for (auto pair_first_ptr : T_table.at(current_trigger_string))
            {
                // All checks already performed, just check if removed
                if (parse.removed(pair_first_ptr)) { continue; }
    
                size_type pair_first_v = (*pair_first_ptr) - 1;
                size_type pair_second_v = 0;
    
                size_type* pair_second_ptr;
    
                if (this->parse.next(pair_first_ptr) != this->parse.end())
                { pair_second_ptr = parse.next(pair_first_ptr); pair_second_v = *pair_second_ptr - 1; }
                else { continue; }
    
                // update entry of T where first appear as second or second appear as a first
                std::string_view first_ts(&(D_prime.at(pair_first_v)[0]), window_length);
                std::string_view second_ts(&(D_prime.at(pair_second_v)[D_prime.at(pair_second_v).size() - window_length]) , window_length);
    
                // new phrase
                size_type merged_phrase_id = 0;
                if (not merged_pairs.contains(std::pair(pair_first_v, pair_second_v)))
                {
                    std::string merged_phrase = D_prime.at(pair_first_v) + D_prime.at(pair_second_v).substr(window_length);
                    if (curr_id > (std::numeric_limits<size_type>::max() - 10)) { spdlog::error("Current id size (bytes) not big enough"); std::exit(EXIT_FAILURE); }
        
                    merged_phrase_id = curr_id++;
                    merged_phrase_id = merged_phrase_id - 1; // compatibility with values from the parse
                    D_prime.d_prime_map.insert(std::pair(merged_phrase_id, merged_phrase));
                    merged_pairs.insert(std::make_pair(std::make_pair(pair_first_v, pair_second_v), merged_phrase_id));
        
                    removed_phrases.insert(pair_first_v);
                    removed_phrases.insert(pair_second_v);
                }
                else
                {
                    merged_phrase_id = merged_pairs.at(std::make_pair(pair_first_v, pair_second_v));
                }
    
                // update parse, second pointer
                *pair_second_ptr = merged_phrase_id + 1; // compatibility with rest of values
    
                // update parse, delete first
                parse.remove(pair_first_ptr);
            }
        }
    }
    
    spdlog::info("Removed {} SIMPLE ts out of {} total", removed_trigger_strings.size(), T_table.size());
    
    // Put cost of removed trigger strings to 0
    for (auto& ts : removed_trigger_strings) { this->T_table.erase(ts); }
    
    // remove phrases from dictionary
    for (auto phrase_id : removed_phrases) { this->D_prime.remove(phrase_id); }
    
    return bytes_removed;
}

std::size_t
vcfbwt::pfp::AuPair::remove_by_cost(std::set<std::string_view>& removed_trigger_strings, int threshold)
{
    spdlog::info("Start compressing with threshold {}", threshold);
    std::size_t bytes_removed = 0;
    
    spdlog::info("Initial TS count: {}", T_table.size());
    if (T_table.size() <= 1) { return 0; }
    
    if (threshold == 0)
    {
        std::pair<int, int> max_cost_trigger_string = priority_queue.get_max();
        threshold = this->window_length + 1;
        spdlog::info("Setting threshold to {}", threshold);
    }
    
    // To update the dictionary at the end
    std::set<size_type> removed_phrases;
    
    // iterate over priority queue
    std::pair<int, int> max_cost_trigger_string = priority_queue.get_max();
    while (max_cost_trigger_string.first > threshold)
    {
        std::string_view& current_trigger_string = this->trigger_string_pq_ids_inv.at(max_cost_trigger_string.second);
        
        // checks for integrity
        if ((removed_trigger_strings.contains(current_trigger_string)) or
        (current_trigger_string[0] == DOLLAR or current_trigger_string[0] == DOLLAR_PRIME))
        {
            // keep iterating
            this->priority_queue.push(max_cost_trigger_string.second, 0);
            max_cost_trigger_string = priority_queue.get_max();
            continue;
        }
        
        // check if any starts and ends with same ts, check for dollars
        bool checks_failed = false;
        for (auto pair_first_ptr : T_table.at(current_trigger_string))
        {
            if (parse.removed(pair_first_ptr)) { continue; }
            
            size_type pair_first_v = (*pair_first_ptr) - 1;
            size_type pair_second_v = 0;
            
            if (this->parse.next(pair_first_ptr) != this->parse.end()) { pair_second_v = (*(parse.next(pair_first_ptr))) - 1; }
            else { continue; }
            
            std::string_view f_ts(&(D_prime.at(pair_first_v)[0]), window_length);
            std::string_view l_ts(&(D_prime.at(pair_second_v)[D_prime.at(pair_second_v).size() - window_length]), window_length);
            
            if (f_ts[0] == DOLLAR
            or (f_ts[0] == DOLLAR_PRIME
            or (l_ts[0] == DOLLAR
            or (l_ts[0] == DOLLAR_PRIME
            or (f_ts == current_trigger_string or l_ts == current_trigger_string)))))
            {
                this->priority_queue.push(max_cost_trigger_string.second, 0);
                max_cost_trigger_string = priority_queue.get_max();
                checks_failed = true;
                break;
            }
        }
        if (checks_failed) { continue; }
        
        removed_trigger_strings.insert(current_trigger_string);
        
        bytes_removed += max_cost_trigger_string.first;
        
        // remove trigger string if cost over threshold
        spdlog::debug("{}\tcost:\t{}\tbytes removed:\t{}\tremoved ts: {}\tT size: {}",
                     current_trigger_string, max_cost_trigger_string.first, bytes_removed, removed_trigger_strings.size(), T_table.size());

        std::map<std::pair<size_type, size_type>, hash_type> merged_pairs;
        std::map<std::string_view, std::set<std::pair<size_type, size_type>>> updates_l_1_add, updates_l_1_sub, updates_r_1_add, updates_r_1_sub;
        std::map<std::string_view, std::set<size_type>> updates_l_3_add, updates_l_3_sub, updates_r_2_add, updates_r_2_sub;
        std::map<std::string_view, int> update_value;
        for (auto pair_first_ptr : T_table.at(current_trigger_string))
        {
            // All checks already performed, just check if removed
            if (parse.removed(pair_first_ptr)) { continue; }

            size_type pair_first_v = (*pair_first_ptr) - 1;
            size_type pair_second_v = 0;

            size_type* pair_second_ptr;

            if (this->parse.next(pair_first_ptr) != this->parse.end())
            { pair_second_ptr = parse.next(pair_first_ptr); pair_second_v = *pair_second_ptr - 1; }
            else { continue; }

            // update entry of T where first appear as second or second appear as a first
            std::string_view first_ts(&(D_prime.at(pair_first_v)[0]), window_length);
            std::string_view second_ts(&(D_prime.at(pair_second_v)[D_prime.at(pair_second_v).size() - window_length]) , window_length);

            // new phrase
            size_type merged_phrase_id = 0;
            if (not merged_pairs.contains(std::pair(pair_first_v, pair_second_v)))
            {
                std::string merged_phrase = D_prime.at(pair_first_v) + D_prime.at(pair_second_v).substr(window_length);
                if (curr_id > (std::numeric_limits<size_type>::max() - 10)) { spdlog::error("Current id size (bytes) not big enough"); std::exit(EXIT_FAILURE); }

                merged_phrase_id = curr_id++;
                merged_phrase_id = merged_phrase_id - 1; // compatibility with values from the parse
                D_prime.d_prime_map.insert(std::pair(merged_phrase_id, merged_phrase));
                merged_pairs.insert(std::make_pair(std::make_pair(pair_first_v, pair_second_v), merged_phrase_id));

                removed_phrases.insert(pair_first_v);
                removed_phrases.insert(pair_second_v);
            }
            else
            {
                merged_phrase_id = merged_pairs.at(std::make_pair(pair_first_v, pair_second_v));
            }

            // update parse, second pointer
            *pair_second_ptr = merged_phrase_id + 1; // compatibility with rest of values

            // update parse, delete first
            parse.remove(pair_first_ptr);

            // prevs
            size_type prev_id = *(parse.prev(pair_second_ptr)) - 1;
            if ( ((first_ts[0] != DOLLAR and first_ts[0] != DOLLAR_PRIME)
            and parse.prev(pair_second_ptr) != parse.end())
            and (not updates_l_1_add[first_ts].contains(std::pair(prev_id, merged_phrase_id))))
            {
                updates_l_1_add[first_ts].insert(std::pair(prev_id, merged_phrase_id));
                update_value[first_ts] += D_prime.at(prev_id).size() + D_prime.at(merged_phrase_id).size() - window_length;

                if (not updates_l_1_sub[first_ts].contains(std::pair(prev_id, pair_first_v)))
                {
                    updates_l_1_sub[first_ts].insert(std::pair(prev_id, pair_first_v));
                    update_value[first_ts] -= D_prime.at(prev_id).size() + D_prime.at(pair_first_v).size() - window_length;
                }

                if (not updates_l_3_add[first_ts].contains(merged_phrase_id))
                {
                    updates_l_3_add[first_ts].insert(merged_phrase_id);
                    update_value[first_ts] -= D_prime.at(merged_phrase_id).size(); // subtract
                }

                if (not updates_l_3_sub[first_ts].contains(pair_first_v))
                {
                    updates_l_3_sub[first_ts].insert(pair_first_v);
                    update_value[first_ts] += D_prime.at(pair_first_v).size(); // add
                }
            }

            // nexts
            size_type next_id = *(parse.next(pair_second_ptr)) - 1;
            if ( ((second_ts[0] != DOLLAR and second_ts[0] != DOLLAR_PRIME)
            and parse.next(pair_second_ptr) != parse.end())
            and (not updates_r_1_add[second_ts].contains(std::pair(merged_phrase_id, next_id))))
            {
                updates_r_1_add[second_ts].insert(std::pair(merged_phrase_id, next_id));
                update_value[second_ts] += D_prime.at(merged_phrase_id).size() + D_prime.at(next_id).size() - window_length;

                if (not updates_r_1_sub[second_ts].contains(std::pair(pair_second_v, next_id)))
                {
                    updates_r_1_sub[second_ts].insert(std::pair(pair_second_v, next_id));
                    update_value[second_ts] -= D_prime.at(pair_second_v).size() + D_prime.at(next_id).size() - window_length;
                }

                if (not updates_r_2_add[second_ts].contains(merged_phrase_id))
                {
                    updates_r_2_add[second_ts].insert(merged_phrase_id);
                    update_value[second_ts] -= D_prime.at(merged_phrase_id).size(); // subtract
                }

                if (not updates_r_2_sub[second_ts].contains(pair_second_v))
                {
                    updates_r_2_sub[second_ts].insert(pair_second_v);
                    update_value[second_ts] += D_prime.at(pair_second_v).size(); //add
                }
            }
        }

        // apply ts cost updates
        for (auto& update : update_value)
        {
            string_view ts = update.first;
            if (removed_trigger_strings.contains(ts)) { continue; }

            int value = update.second;

            int ts_index = this->trigger_string_pq_ids.at(ts);
            int old_cost = priority_queue.get_key(ts_index);
//            int sb_cost = cost_of_removing_trigger_string(ts);
//            if (old_cost != 0 and sb_cost != old_cost - value)
//            {
//                spdlog::error("should be: {} instead is {}", sb_cost, old_cost - value);
//            }
            if (old_cost != 0) { this->priority_queue.push(ts_index, old_cost - value); }
        }
        
        this->priority_queue.push(max_cost_trigger_string.second, 0);
        
        // keep iterating
        max_cost_trigger_string = priority_queue.get_max();
        
        this->T_table.erase(current_trigger_string);
    }
    
    // remove phrases from dictionary
    for (auto phrase_id : removed_phrases) { this->D_prime.remove(phrase_id); }
    
    return bytes_removed;
}

void
vcfbwt::pfp::AuPair::close()
{
    if (this->closed) { return; }
    this->closed = true;
    
    // statistics
    std::size_t final_size = 0;
    
    // sort dictionary
    std::vector<std::pair<std::reference_wrapper<std::string>, hash_type>> sorted_phrases;
    
    for (std::size_t i = 0; i < D_prime.d_prime_vector.size(); i++)
    {
        if (!D_prime.d_prime_vector[i].empty())
        {
            sorted_phrases.emplace_back(std::ref(D_prime.d_prime_vector[i]), i);
        }
    }
    for (auto& d_pair : D_prime.d_prime_map)
    {
        sorted_phrases.emplace_back(std::pair(std::ref(d_pair.second), d_pair.first));
    }
    
    std::sort(sorted_phrases.begin(), sorted_phrases.end(), ref_smaller);
    
    std::unordered_map<hash_type, size_type> hash_to_rank;
    for (std::size_t i = 0; i < sorted_phrases.size(); i++) { hash_to_rank[sorted_phrases[i].second] = i + 1; }
    
    spdlog::info("AuPair: writing dictionary to disk NOT COMPRESSED");
    std::string dict_file_name = this->in_prefix + EXT::N_DICT;
    std::ofstream dict(dict_file_name);
    
    for (auto& sorted_phrase : sorted_phrases)
    {
        dict.write(sorted_phrase.first.get().c_str(), sorted_phrase.first.get().size());
        dict.put(ENDOFWORD);
        
        final_size += sorted_phrase.first.get().size();
    }
    dict.put(ENDOFDICT);
    vcfbwt::DiskWrites::update(dict.tellp()); // Disk Stats
    dict.close();

    if (this->compress_dictionary)
    {
        spdlog::info("Main parser: writing dictionary on disk COMPRESSED");
        std::ofstream dicz(this->in_prefix + EXT::N_DICT_COMPRESSED);
        std::ofstream lengths(this->in_prefix + EXT::N_DICT_COMPRESSED_LENGTHS);

        for (size_type i = 0; i < sorted_phrases.size(); i++)
        {
            std::size_t shift = 1; // skip dollar on first phrase
            if (i != 0) { shift = this->window_length; }
            dicz.write(sorted_phrases.at(i).first.get().c_str() + shift,
                       sorted_phrases.at(i).first.get().size() - shift);
            int32_t len = sorted_phrases.at(i).first.get().size() - shift;
            lengths.write((char*) &len, sizeof(int32_t));
        }

        vcfbwt::DiskWrites::update(dicz.tellp()); // Disk Stats
        dicz.close();

        vcfbwt::DiskWrites::update(lengths.tellp()); // Disk Stats
        lengths.close();
    }
    
    spdlog::info("AuPair: writing parse and .last to disk");
    std::string parse_file_name = this->in_prefix + EXT::N_PARSE;
    std::ofstream parse_file(parse_file_name);
    std::string last_file_name = this->in_prefix + EXT::N_LAST;
    std::ofstream last_file(last_file_name);
    std::string sai_file_name = this->in_prefix + EXT::N_SAI;
    std::ofstream sai_file(sai_file_name);
    
    std::size_t pos_for_sai = 0;
    
    std::vector<size_type> occurrences(sorted_phrases.size(), 0);
    std::string occ_file_name = this->in_prefix + EXT::N_OCC;
    std::ofstream occ_file(occ_file_name);
    
    auto* parse_it = this->parse.begin();
    while (parse_it != this->parse.end())
    {
        if (hash_to_rank.find((*parse_it) - 1) == hash_to_rank.end())
        {
            spdlog::debug("Phrase {} not in the dictionary after compressing", (*parse_it) - 1);
        }
        else
        {
            parse_file.write((char*) &(hash_to_rank.at((*parse_it) - 1)), sizeof(size_type));
            occurrences[hash_to_rank.at((*parse_it) - 1) - 1] += 1;
    
            const std::string& dict_string = sorted_phrases[hash_to_rank.at((*parse_it) - 1) - 1].first.get();
            last_file.put(dict_string[(dict_string.size() - window_length) - 1]);
    
            if (pos_for_sai == 0) { pos_for_sai = dict_string.size() - 1; } // -1 is for the initial $ of the first word
            else { pos_for_sai += dict_string.size() - window_length; }
            sai_file.write((char*) &pos_for_sai, IBYTES);
            
            final_size += sizeof(size_type);
        }
        parse_it = this->parse.next(parse_it);
    }
    
    vcfbwt::DiskWrites::update(parse_file.tellp()); // Disk Stats
    parse_file.close();
    
    vcfbwt::DiskWrites::update(last_file.tellp());
    last_file.close();
    
    vcfbwt::DiskWrites::update(sai_file.tellp());
    sai_file.close();
    
    occ_file.write(reinterpret_cast<const char*>(occurrences.data()), sizeof(size_type) * occurrences.size());
    vcfbwt::DiskWrites::update(occ_file.tellp()); // Disk Stats
    occ_file.close();
    
    spdlog::info("Compression ratio: {}%", (100 - std::size_t(((double(final_size) / double(this->initial_size)) * 100))));
}

//------------------------------------------------------------------------------