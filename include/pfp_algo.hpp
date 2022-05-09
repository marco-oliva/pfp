//
//  pfp_algo.hpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#ifndef pfp_algo_hpp
#define pfp_algo_hpp

#include <vector>
#include <unordered_map>
#include <set>
#include <iostream>
#include <fstream>
#include <vcf.hpp>
#include <utils.hpp>
#include <internals.hpp>

namespace vcfbwt
{
namespace pfp
{

//------------------------------------------------------------------------------

enum SPECIAL_TYPES
{
    ENDOFDICT = 0,
    ENDOFWORD = 1,
    DOLLAR = 2,
    DOLLAR_SEQUENCE = 4,
    DOLLAR_PRIME = 5
};

//------------------------------------------------------------------------------

template <typename data_type>
class Dictionary
{
public:

    vcfbwt::size_type insertions_safe_guard = 1000;

    std::mutex dictionary_mutex;

    bool sorted = false;
    
    struct DictionaryEntry
    {
        std::vector<data_type> phrase;
    
        DictionaryEntry() = default;
        explicit DictionaryEntry(const std::vector<data_type>& s) : phrase(s) {}
    };
    
    std::unordered_map<hash_type, DictionaryEntry> hash_string_map;
    std::vector<std::pair<std::reference_wrapper<std::vector<data_type>>, hash_type>> sorted_phrases;
    std::unordered_map<hash_type, size_type> hash_to_ranks;
    
    void sort()
    {
        // lock the dictionary
        std::lock_guard<std::mutex> guard(dictionary_mutex);

        // sort the dictionary
        for (auto& entry : this->hash_string_map)
        {
            this->sorted_phrases.emplace_back(std::ref(entry.second.phrase), entry.first);
        }
        std::sort(sorted_phrases.begin(), sorted_phrases.end(), ref_smaller<data_type>);

        // insert in hashmap
        this->hash_to_ranks.clear();
        for (size_type i = 0; i < sorted_phrases.size(); i++)
        {
            hash_to_ranks.insert(std::make_pair(sorted_phrases[i].second, i + 1)); // 1 based
        }

        this->sorted = true;
    }
    
    Dictionary() = default;
    
    hash_type add(const std::vector<data_type>& phrase)
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

        if (this->size() >= (std::numeric_limits<size_type>::max() - insertions_safe_guard))
        { spdlog::error("Dictionary too big for type {}", typeid(size_type).name()); std::exit(EXIT_FAILURE); }

        return phrase_hash;
    }

    hash_type check_and_add(const std::vector<data_type>& phrase)
    {
        // lock the dictionary
        std::lock_guard<std::mutex> guard(dictionary_mutex);

        // Check if present
        hash_type phrase_hash = string_hash((const char*) &(phrase[0]), phrase.size() * sizeof(data_type));
        const auto& ptr = hash_string_map.find(phrase_hash);

        if ((ptr != hash_string_map.end()) and (ptr->second.phrase != phrase))
        {
            spdlog::error("Hash collision! Hash already in the dictionary for a different phrase");
            std::exit(EXIT_FAILURE);
        }
        else if (ptr != hash_string_map.end()) { return phrase_hash; }

        this->sorted = false;

        DictionaryEntry entry(phrase);
        hash_string_map.insert(std::make_pair(phrase_hash, entry));

        if (this->size() >= (std::numeric_limits<size_type>::max() - insertions_safe_guard))
        { spdlog::error("Dictionary too big for type {}", typeid(size_type).name()); std::exit(EXIT_FAILURE); }

        return phrase_hash;
    }

    hash_type get(const std::vector<data_type>& phrase) const
    {
        return string_hash(&(phrase[0]), phrase.size());
    }

    bool contains(const std::vector<data_type>& phrase)
    {
        // lock the dictionary
        std::lock_guard<std::mutex> guard(dictionary_mutex);

        hash_type phrase_hash = string_hash((const char*) &(phrase[0]), phrase.size() * sizeof(data_type));
        const auto& ptr = hash_string_map.find(phrase_hash);

        return ((ptr != hash_string_map.end()) and (ptr->second.phrase == phrase));
    }

    size_type hash_to_rank(hash_type hash)
    {
        if (not this->sorted) { sort(); }
        auto it = hash_to_ranks.find(hash);
        if (it != hash_to_ranks.end()) { return it->second; }
        else { spdlog::error("Something went wrong"); std::exit(EXIT_FAILURE); }
    }
    
    size_type size() const { return hash_string_map.size(); }
    
    const std::vector<data_type>& sorted_entry_at(std::size_t i)
    {
        std::lock_guard<std::mutex> guard(dictionary_mutex);
        if (not this->sorted) { sort(); } return sorted_phrases[i].first.get();
    }

    friend class ParserVCF;
};

//------------------------------------------------------------------------------

struct Params
{
    hash_type p = 100;
    hash_type w =  10;
    bool compress_dictionary = false;
    bool use_acceleration = false;
    bool print_out_statistics_csv = false;
    bool compute_occurrences = true;
    bool auPair = false;
    std::string ignore_ts_file;
    int32_t integers_shift = 10;
};

struct Statistics
{
    std::size_t parse_length = 0;
    std::size_t total_dictionary_length = 0;
    std::size_t num_of_phrases_dictionary = 0;
};

class ReferenceParse
{

public :
    Dictionary<vcfbwt::char_type> dictionary;
    std::vector<hash_type> parse;
    std::vector<std::size_t> trigger_strings_position; // position of first char of each trigger string
    std::set<hash_type> to_ignore_ts_hash;
    
    const Params& params;
    
    void init(const std::string& reference);
    
    ReferenceParse(const std::string& reference, const Params& pms) : params(pms) { this->init(reference); }
    const hash_type& operator[](size_type i) const { return this->parse[i]; }
};

// Create a parse on disk
class ParserVCF
{

private:
    
    std::ofstream out_file;
    std::string out_file_prefix;
    std::string out_file_name;
    std::string tmp_out_file_name;
    std::vector<std::string> samples_processed;
    
    Params params;
    Statistics statistics;
    
    ReferenceParse* reference_parse = nullptr;
    Dictionary<vcfbwt::char_type>* dictionary = nullptr;

    // Shorthands
    hash_type w, p;
    std::size_t parse_size = 0;
    size_type tags;
    
    bool closed = false;
    
    std::vector<std::reference_wrapper<ParserVCF>> registered_workers;
    
    std::size_t working_genotype = 0;
    
public:
    
    enum tags
    {
        MAIN = 1,
        WORKER = 2,
        COMPRESSED = 4,
        UNCOMPRESSED = 8
    };
    
    
    void init(const Params& params, const std::string& out_prefix, ReferenceParse& rp, std::size_t t = MAIN | UNCOMPRESSED);
    
    ParserVCF(const Params& params, const std::string& out_prefix, ReferenceParse& rp, std::size_t t = MAIN | UNCOMPRESSED)
    {
        if (out_prefix.empty()) { this->init(params, "out", rp, t); }
        else { this->init(params, out_prefix, rp, t); }
    }
    
    ParserVCF() = default;
    
    ~ParserVCF()
    {
        close();
        size_type total_length = 0;
        for (auto& entry : dictionary->hash_string_map) { total_length += entry.second.phrase.size(); }
    
        // Fill out statistics
        this->statistics.parse_length = this->parse_size;
        this->statistics.num_of_phrases_dictionary = this->dictionary->sorted_phrases.size();
        this->statistics.total_dictionary_length = total_length;
        std::string name;
        if (tags & WORKER) { name = "Worker"; } else { name = "Main"; }
        spdlog::info("{} -\tParse size: {}\tDic Size: {} Dic Total Length: {}", name ,parse_size, dictionary->size(), total_length);
        
        if (params.print_out_statistics_csv and (tags & MAIN))
        {
            std::ofstream csv(out_file_prefix + ".csv");
            csv << "w,p,parse_lenght,dict_phrases,dict_tot_length\n";
            csv << params.w << ",";
            csv << params.p << ",";
            csv << statistics.parse_length << ",";
            csv << statistics.num_of_phrases_dictionary << ",";
            csv << statistics.total_dictionary_length;
            csv << std::endl;
            csv.close();
        }
    }
    
    const std::string& get_file_name() const { return this->out_file_name; }
    const Statistics& get_statistics() const { return this->statistics; }
    void register_worker(ParserVCF& parser) { this->registered_workers.push_back(std::ref(parser)); }
    void set_working_genotype(std::size_t genotype) { this->working_genotype = genotype; }
    
    void operator()(const Sample& sample);
    void close();
};

//------------------------------------------------------------------------------

class ParserFasta
{

private:
    
    std::ofstream out_file;
    std::string out_file_prefix;
    std::string out_file_name;
    std::string in_file_path;
    std::vector<std::string> sequences_processed;
    
    Params params;
    Statistics statistics;
    
    Dictionary<vcfbwt::char_type> dictionary;
    
    hash_type w, p;
    std::size_t parse_size = 0;
    
    bool closed = false;
    
public:
    
    void init(const Params& params, const std::string& prefix);
    
    ParserFasta(const Params& params, const std::string& file_path, const std::string& out_prefix)
    {
        this->w = params.w; this->p = params.p;
        this->in_file_path = file_path;
        if (out_prefix.empty()) { this->init(params, "out"); }
        else { this->init(params, out_prefix); }
    }
    
    ParserFasta() = default;
    
    ~ParserFasta()
    {
        close();
        size_type total_length = 0;
        for (auto& entry : dictionary.hash_string_map) { total_length += entry.second.phrase.size(); }
        
        // Fill out statistics
        this->statistics.parse_length = this->parse_size;
        this->statistics.num_of_phrases_dictionary = this->dictionary.sorted_phrases.size();
        this->statistics.total_dictionary_length = total_length;
        std::string name = "Parser Fasta";
        spdlog::info("{} -\tParse size: {}\tDic Size: {} Dic Total Length: {}", name ,parse_size, dictionary.size(), total_length);
        
        if (params.print_out_statistics_csv)
        {
            std::ofstream csv(out_file_prefix + ".csv");
            csv << "w,p,parse_lenght,dict_phrases,dict_tot_length\n";
            csv << params.w << ",";
            csv << params.p << ",";
            csv << statistics.parse_length << ",";
            csv << statistics.num_of_phrases_dictionary << ",";
            csv << statistics.total_dictionary_length;
            csv << std::endl;
            csv.close();
        }
    }
    
    const std::string& get_file_name() const { return this->out_file_name; }
    const Statistics& get_statistics() const { return this->statistics; }
    
    void operator()();
    void close();
};

//------------------------------------------------------------------------------

class ParserText
{

private:
    
    std::ofstream out_file;
    std::string out_file_prefix;
    std::string out_file_name;
    std::string in_file_path;
    std::vector<std::string> sequences_processed;
    
    Params params;
    Statistics statistics;
    
    Dictionary<vcfbwt::char_type> dictionary;
    
    hash_type w, p;
    std::size_t parse_size = 0;
    
    bool closed = false;

public:
    
    void init(const Params& params, const std::string& prefix);
    
    ParserText(const Params& params, const std::string& file_path, const std::string& out_prefix)
    {
        this->w = params.w; this->p = params.p;
        this->in_file_path = file_path;
        if (out_prefix.empty()) { this->init(params, "out"); }
        else { this->init(params, out_prefix); }
    }
    
    ParserText() = default;
    
    ~ParserText()
    {
        close();
        size_type total_length = 0;
        for (auto& entry : dictionary.hash_string_map) { total_length += entry.second.phrase.size(); }
        
        // Fill out statistics
        this->statistics.parse_length = this->parse_size;
        this->statistics.num_of_phrases_dictionary = this->dictionary.sorted_phrases.size();
        this->statistics.total_dictionary_length = total_length;
        std::string name = "Parser Text";
        spdlog::info("{} -\tParse size: {}\tDic Size: {} Dic Total Length: {}", name ,parse_size, dictionary.size(), total_length);
        
        if (params.print_out_statistics_csv)
        {
            std::ofstream csv(out_file_prefix + ".csv");
            csv << "w,p,parse_lenght,dict_phrases,dict_tot_length\n";
            csv << params.w << ",";
            csv << params.p << ",";
            csv << statistics.parse_length << ",";
            csv << statistics.num_of_phrases_dictionary << ",";
            csv << statistics.total_dictionary_length;
            csv << std::endl;
            csv.close();
        }
    }
    
    const std::string& get_file_name() const { return this->out_file_name; }
    const Statistics& get_statistics() const { return this->statistics; }
    
    void operator()();
    void close();
};

//------------------------------------------------------------------------------


class ParserIntegers
{

private:

    std::ofstream out_file;
    std::string out_file_prefix;
    std::string out_file_name;
    std::string in_file_path;
    std::vector<std::string> sequences_processed;

    Params params;
    Statistics statistics;

    Dictionary<int32_t> dictionary;

    hash_type w, p;
    std::size_t parse_size = 0;

    bool closed = false;

public:

    void init(const Params& params, const std::string& prefix);

    ParserIntegers(const Params& params, const std::string& file_path, const std::string& out_prefix)
    {
        this->w = params.w; this->p = params.p;
        this->in_file_path = file_path;
        if (out_prefix.empty()) { this->init(params, "out"); }
        else { this->init(params, out_prefix); }
    }

    ParserIntegers() = default;

    ~ParserIntegers()
    {
        close();
        size_type total_length = 0;
        for (auto& entry : dictionary.hash_string_map) { total_length += entry.second.phrase.size() * sizeof(int32_t); }

        // Fill out statistics
        this->statistics.parse_length = this->parse_size;
        this->statistics.num_of_phrases_dictionary = this->dictionary.sorted_phrases.size();
        this->statistics.total_dictionary_length = total_length;
        std::string name = "Parser Integers";
        spdlog::info("{} -\tParse size: {}\tDic Size: {} Dic Total Length: {}", name ,parse_size, dictionary.size(), total_length);

        if (params.print_out_statistics_csv)
        {
            std::ofstream csv(out_file_prefix + ".csv");
            csv << "w,p,parse_lenght,dict_phrases,dict_tot_length\n";
            csv << params.w << ",";
            csv << params.p << ",";
            csv << statistics.parse_length << ",";
            csv << statistics.num_of_phrases_dictionary << ",";
            csv << statistics.total_dictionary_length;
            csv << std::endl;
            csv.close();
        }
    }

    const std::string& get_file_name() const { return this->out_file_name; }
    const Statistics& get_statistics() const { return this->statistics; }

    void operator()();
    void close();
};

//------------------------------------------------------------------------------

template <typename data_type>
class ParserUtils
{
public:
    
    static void read_parse(std::string parse_file_name, std::vector<size_type>& parse)
    {
        std::ifstream parse_file(parse_file_name, std::ios::binary);
        if (not parse_file.is_open()) { spdlog::error("Error opening file: {}", parse_file_name); std::exit(EXIT_FAILURE); }

        while (not parse_file.eof()) { size_type i; parse_file.read((char*) &i, sizeof(size_type)); parse.push_back(i); }

        // remove last entry, apparently it reads the last twice, strano
        parse.pop_back();
    }

    static void read_dictionary(std::string dic_file_name, std::vector<std::vector<data_type>>& dictionary_vector)
    {
        std::ifstream dic_file(dic_file_name);
        if (not dic_file.is_open()) { spdlog::error("Error opening file: {}", dic_file_name); std::exit(EXIT_FAILURE); }

        data_type c = 5;
        while (c != ENDOFDICT)
        {
            std::vector<data_type> phrase;
            while ((dic_file.read((char*) &c, sizeof(data_type)).good()) and (c != ENDOFWORD and c != ENDOFDICT))
            {
                phrase.emplace_back(c);
            }
            if (phrase.size() != 0) { dictionary_vector.push_back(phrase); }
        }
    }

    static void merge(const std::string& left_prefix, const std::string& right_prefix, const std::string& out_prefix, const Params& params)
    {
        // Merged Dictionary and Parse
        Dictionary<data_type> dictionary;

        // Read Left and Right Dictionary
        spdlog::info("Loading dictionaries from disk");
        std::vector<std::vector<data_type>> left_dictionary;   read_dictionary(left_prefix + EXT::DICT, left_dictionary);
        std::vector<std::vector<data_type>> right_dictionary; read_dictionary(right_prefix + EXT::DICT, left_dictionary);

        // Read Left parse and substitute ranks with hash again storing in the merged dictionary
        spdlog::info("Adjusting ranks");
        std::error_code error;
        mio::mmap_sink rw_mmap = mio::make_mmap_sink(left_prefix + EXT::PARSE, 0, mio::map_entire_file, error);
        if (error) { spdlog::error(error.message()); std::exit(EXIT_FAILURE); }

        for (size_type i = 0; i < (rw_mmap.size() / sizeof(hash_type)); i++)
        {
            hash_type rank;
            std::memcpy(&rank, rw_mmap.data() + (i * sizeof(hash_type)), sizeof(hash_type));

            std::vector<data_type>& phrase = left_dictionary[rank];
            hash_type hash = 0;
            if (dictionary.contains(phrase))    { hash = dictionary.get(phrase); }
            else                                { hash = dictionary.add(phrase); }

            std::memcpy(rw_mmap.data() + (i * sizeof(hash_type)), &hash, sizeof(hash_type));
        }
        rw_mmap.unmap(); left_dictionary.resize(0);

        // Read right parse, changing first phrase;
        rw_mmap = mio::make_mmap_sink(left_prefix + EXT::PARSE, 0, mio::map_entire_file, error);
        if (error) { spdlog::error(error.message()); std::exit(EXIT_FAILURE); }

        for (size_type i = 0; i < (rw_mmap.size() / sizeof(hash_type)); i++)
        {
            hash_type rank;
            std::memcpy(&rank, rw_mmap.data() + (i * sizeof(hash_type)), sizeof(hash_type));

            std::vector<data_type>& phrase = right_dictionary[rank];

            if (rank == 1) { phrase[0] = DOLLAR_PRIME; phrase.insert(phrase.begin(), params.w - 1, DOLLAR_PRIME); }

            hash_type hash = 0;
            if (dictionary.contains(phrase))    { hash = dictionary.get(phrase); }
            else                                { hash = dictionary.add(phrase); }

            std::memcpy(rw_mmap.data() + (i * sizeof(hash_type)), &hash, sizeof(hash_type));
        }
        rw_mmap.unmap(); right_dictionary.resize(0);

        // Sort the Dictionary
        spdlog::info("Sorting merged dictionary");
        std::vector<std::pair<std::vector<data_type>, hash_type>> entries;
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
        std::ofstream out_parse(out_prefix + EXT::PARSE, std::ios::binary);

        // Left
        mio::mmap_source rl_mmap; rl_mmap.map(left_prefix + EXT::PARSE, error);
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
        mio::mmap_source rr_mmap; rr_mmap.map(left_prefix + EXT::PARSE, error);
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
        std::string dict_file_name = out_prefix + EXT::DICT;
        std::ofstream dict(dict_file_name);

        for (auto& entry : entries)
        {
            dict.write(entry.first.data(), entry.first.size());
            dict.write((char*) ENDOFWORD, sizeof(data_type));
        }

        dict.write((char*) ENDOFDICT, sizeof(data_type));
        dict.close();
    }
};

} // end namespace pfp
} // end namespace vcfbwt

#endif //pfp_algo_hpp
