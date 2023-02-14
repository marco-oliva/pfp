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

    std::size_t insertions_safe_guard = 1000;

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
            spdlog::error("Dictionary::addHash collision! Hash already in the dictionary");
            std::exit(EXIT_FAILURE);
        }

        DictionaryEntry entry(phrase);
        hash_string_map.insert(std::make_pair(phrase_hash, entry));

        if (this->size() >= (std::numeric_limits<size_type>::max() - insertions_safe_guard))
        { spdlog::error("Dictionary::add Dictionary too big for type {}", typeid(size_type).name()); std::exit(EXIT_FAILURE); }

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
            spdlog::error("Dictionary::check_and_add Hash collision! Hash already in the dictionary for a different phrase");
            std::exit(EXIT_FAILURE);
        }
        else if (ptr != hash_string_map.end()) { return phrase_hash; }

        this->sorted = false;

        DictionaryEntry entry(phrase);
        hash_string_map.insert(std::make_pair(phrase_hash, entry));

        if (this->size() >= (std::numeric_limits<size_type>::max() - insertions_safe_guard))
        { spdlog::error("Dictionary::check_and_add Dictionary too big for type {}", typeid(size_type).name()); std::exit(EXIT_FAILURE); }

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
        else { spdlog::error("Dictionary::hash_to_rank hash requested not in the dictionary. hash: {}", hash); std::exit(EXIT_FAILURE); }
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
    bool vcf_acgt_only = false;
    bool print_out_statistics_csv = false;
    bool output_occurrences = false;
    bool output_sai = false;
    bool output_last = false;
    uint32_t integers_shift = 10;
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
    
    const Params& params;
    
    void init(const std::string& reference);
    
    ReferenceParse(const std::string& reference, const Params& pms) : params(pms) { this->init(reference); }
    const hash_type& operator[](std::size_t i) const { return this->parse[i]; }
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
    hash_type w = params.w, p = params.p;
    std::size_t parse_size = 0;
    std::size_t tags;
    
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
        std::size_t total_length = 0;
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
    std::string tmp_out_file_name;
    std::string in_file_path;
    std::vector<std::string> sequences_processed;
    
    Params params;
    Statistics statistics;
    
    Dictionary<vcfbwt::char_type> dictionary;
    
    hash_type w = params.w, p = params.p;
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
        std::size_t total_length = 0;
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
    std::string tmp_out_file_name;
    std::string in_file_path;
    std::vector<std::string> sequences_processed;
    
    Params params;
    Statistics statistics;
    
    Dictionary<vcfbwt::char_type> dictionary;
    
    hash_type w = params.w, p = params.p;
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
        std::size_t total_length = 0;
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
    std::string tmp_out_file_name;
    std::string in_file_path;
    std::vector<std::string> sequences_processed;

    Params params;
    Statistics statistics;

    Dictionary<uint32_t> dictionary;
    
    hash_type w = params.w, p = params.p;
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
        std::size_t total_length = 0;
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

#include <sys/stat.h>

template <typename data_type>
class ParserUtils
{
public:
    
    static void read_parse(std::string parse_file_name, std::vector<size_type>& parse)
    {
        // reserve memory for the parse
        struct stat64 stat_buf;
        int rc = stat64(parse_file_name.c_str(), &stat_buf);
        if (rc == 0) { parse.reserve((stat_buf.st_size / sizeof(size_type)) + 1); }
        
        std::ifstream parse_file(parse_file_name, std::ios::binary);
        if (not parse_file.is_open()) { spdlog::error("Error opening file: {}", parse_file_name); std::exit(EXIT_FAILURE); }
        
        while (not parse_file.eof()) { size_type i; parse_file.read((char*) &i, sizeof(size_type)); parse.push_back(i); }

        // remove last entry, apparently it reads the last twice
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

        // Out parse, tmp
        std::string tmp_out_file_name = TempFile::getName("parse");
        std::ofstream tmp_out_parse(tmp_out_file_name);
        std::size_t parse_size = 0;

        // Read Left and Right Dictionary
        spdlog::info("Loading dictionaries from disk");
        std::vector<std::vector<data_type>> left_dictionary;
        read_dictionary(left_prefix + EXT::DICT, left_dictionary);
        std::vector<std::vector<data_type>> right_dictionary;
        read_dictionary(right_prefix + EXT::DICT, right_dictionary);

        // Read Left parse and substitute ranks with hash again storing in the merged dictionary
        spdlog::info("Iterating over left parse");
        std::ifstream left_parse(left_prefix + EXT::PARSE);
        if (not left_parse.is_open()) { spdlog::error("Failed to open {}", left_prefix + EXT::PARSE); std::exit(EXIT_FAILURE); }
        
        size_type pel = 0;
        hash_type last_phrase_hash;
        while(left_parse.read((char*) &pel, sizeof(size_type)))
        {
            std::vector<data_type>& phrase = left_dictionary[pel - 1];
            
            if (phrase[phrase.size() - 1] == DOLLAR)
            {
                // get rid of the dollars, it already has DOLLAR_PRIME and DOLLAR_SEQUENCE
                phrase.resize(phrase.size() - params.w);
            }
            
            hash_type hash = dictionary.check_and_add(phrase);
            last_phrase_hash = hash;

            tmp_out_parse.write((char*) (&hash), sizeof(hash_type));
            parse_size += 1;
        }
        left_parse.close();
        left_dictionary.resize(0);

        // Read right parse, changing first phrase
        spdlog::info("Iterating over right parse");
        std::ifstream right_parse(right_prefix + EXT::PARSE);
        if (not right_parse.is_open()) { spdlog::error("Failed to open {}", right_prefix + EXT::PARSE); std::exit(EXIT_FAILURE); }

        size_type per = 0;
        while(right_parse.read((char*) &per, sizeof(size_type)))
        {
            std::vector<data_type>& phrase = right_dictionary[per - 1];
            if (per == 1)
            {
                phrase[0] = DOLLAR_SEQUENCE;
                phrase.insert(phrase.begin(), params.w - 1, DOLLAR_PRIME);
            }
            hash_type hash = dictionary.check_and_add(phrase);

            tmp_out_parse.write((char*) (&hash), sizeof(hash_type));
            parse_size += 1;
        }
        right_parse.close();
        right_dictionary.resize(0);

        vcfbwt::DiskWrites::update(tmp_out_parse.tellp());
        tmp_out_parse.close();

        spdlog::info("Merge: Replacing hash values with ranks.");
        
        if (parse_size != 0)
        {
            std::ofstream out_ranks(out_prefix + EXT::PARSE);
            if (not out_ranks.is_open()) { spdlog::error("Can't open {}", out_prefix + EXT::PARSE); std::exit(EXIT_FAILURE); }

            std::ifstream in_hash(tmp_out_file_name);
            if (not in_hash.is_open()) { spdlog::error("Can't open {}", tmp_out_file_name); std::exit(EXIT_FAILURE); }

            for (std::size_t i = 0; i < parse_size; i++)
            {
                hash_type hash;
                in_hash.read((char*) &hash, sizeof(hash_type));
                size_type rank = dictionary.hash_to_rank(hash);
                out_ranks.write((char*) &rank, sizeof(size_type));
            }
            in_hash.close();
            vcfbwt::DiskWrites::update(out_ranks.tellp());
            out_ranks.close();
        
        }
        // Print dicitionary on disk
        spdlog::info("Merge: writing dictionary on disk NOT COMPRESSED");
        std::string dict_file_name = out_prefix + EXT::DICT;
        std::ofstream dict(dict_file_name);

        for (size_type i = 0; i < dictionary.size(); i++)
        {
            dict.write((char*) dictionary.sorted_entry_at(i).data(), dictionary.sorted_entry_at(i).size());
            dict.put(ENDOFWORD);
        }
        dict.put(ENDOFDICT);

        vcfbwt::DiskWrites::update(dict.tellp()); // Disk Stats
        dict.close();
    }
};


template <typename data_type>
class PropertiesWriter
{
private:
    const Params& params;
    std::string pfp_prefix;

public:
    PropertiesWriter(const std::string& prefix, const Params& pms) : pfp_prefix(prefix), params(pms) { }
    
    void write()
    {
        std::string parse_path = this->pfp_prefix + EXT::PARSE;
        std::string dict_path = this->pfp_prefix + EXT::DICT;
        
        if (not (params.output_occurrences or params.output_last or params.output_sai or params.compress_dictionary))
        { spdlog::info("No properties requested."); return; }
        
        // read in dictionary
        spdlog::info("Loading dictionary from disk.");
        std::vector<std::vector<data_type>> dictionary;
        ParserUtils<data_type>::read_dictionary(dict_path, dictionary);
        
        // output compressed dictionary if needed
        if (params.compress_dictionary)
        {
            spdlog::info("Writing dictionary on disk COMPRESSED");
            std::ofstream dicz(this->pfp_prefix + EXT::DICT_COMPRESSED);
            std::ofstream lengths(this->pfp_prefix + EXT::DICT_COMPRESSED_LENGTHS);
            
            for (std::size_t i = 0; i < dictionary.size(); i++)
            {
                std::size_t shift = 1; // skip dollar on first phrase
                if (i != 0) { shift = this->params.w; }
                dicz.write((char*) (dictionary[i].data() + shift),
                           (dictionary[i].size() - shift) * sizeof(data_type));
                uint32_t len = dictionary[i].size() - shift;
                lengths.write((char*) &len, sizeof(uint32_t));
            }
            
            vcfbwt::DiskWrites::update(dicz.tellp()); // Disk Stats
            dicz.close();
            
            vcfbwt::DiskWrites::update(lengths.tellp()); // Disk Stats
            lengths.close();
        }
        
        std::string last_file_name = this->pfp_prefix + EXT::LAST;
        std::ofstream last_file;
        if (params.output_last) { last_file.open(last_file_name); }
        
        std::string sai_file_name = this->pfp_prefix + EXT::SAI;
        std::ofstream sai_file;
        if (params.output_sai) { sai_file.open(sai_file_name); }
        
        std::vector<long_type> occurrences(dictionary.size(), 0);
        std::string occ_file_name = this->pfp_prefix + EXT::OCC;
        std::ofstream occ_file;
        if (params.output_occurrences) { occ_file.open(occ_file_name); }
        
        // read in parse and output .occ, .last and .sai if needed
        spdlog::info("Read in parse and output properties");
        std::vector<size_type> parse;
        ParserUtils<data_type>::read_parse(parse_path, parse);
    
        std::size_t pos_for_sai = 0;
        if (not parse.empty())
        {
            for (std::size_t i = 0; i < parse.size(); i++)
            {
                size_type rank = parse[i];
                occurrences[rank - 1] += 1;
            
                const std::vector<data_type>& dict_string = dictionary[rank - 1];
                if (params.output_last)
                {
                    last_file.put(dict_string[(dict_string.size() - this->params.w) - 1]);
                }
            
                if (params.output_sai)
                {
                    if (pos_for_sai == 0) { pos_for_sai = dict_string.size() - 1; } // -1 is for the initial $ of the first word
                    else { pos_for_sai += dict_string.size() - this->params.w; }
                    sai_file.write((char*) &pos_for_sai, IBYTES);
                }
            }
        }
    
        if (params.output_last)
        {
            vcfbwt::DiskWrites::update(last_file.tellp());
            last_file.close();
        }
    
        if (params.output_sai)
        {
            vcfbwt::DiskWrites::update(sai_file.tellp());
            sai_file.close();
        }
    
        if(params.output_occurrences)
        {
            spdlog::info("Writing occurrences to file");
        
            for (std::size_t i = 0; i < occurrences.size(); i++)
            {
                if (parse.size() < std::numeric_limits<short_type>::max())
                {
                    short_type to_write = occurrences[i];
                    occ_file.write((char*)&to_write, sizeof(short_type));
                }
                else
                {
                    long_type to_write = occurrences[i];
                    occ_file.write((char*)&to_write, sizeof(long_type));
                }
            }
        
        
            vcfbwt::DiskWrites::update(occ_file.tellp()); // Disk Stats
            occ_file.close();
        }
    }
};

} // end namespace pfp
} // end namespace vcfbwt

#endif //pfp_algo_hpp
