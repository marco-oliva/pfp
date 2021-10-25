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

class Dictionary
{
public:

    vcfbwt::size_type insertions_safe_guard = 1000;

    std::mutex dictionary_mutex;

    bool sorted = false;
    
    struct DictionaryEntry
    {
        std::string phrase;
    
        DictionaryEntry() = default;
        DictionaryEntry(const std::string& s) : phrase(s) {}
    };
    
    std::unordered_map<hash_type, DictionaryEntry> hash_string_map;
    std::vector<std::pair<std::reference_wrapper<std::string>, hash_type>> sorted_phrases;
    std::unordered_map<hash_type, size_type> hash_to_ranks;
    
    void sort();
    
    Dictionary() = default;
    
    hash_type add(const std::string& phrase);
    hash_type check_and_add(const std::string& phrase);
    hash_type get(const std::string& phrase) const;
    bool contains(const std::string& phrase);

    size_type hash_to_rank(hash_type hash);
    
    size_type size() const { return hash_string_map.size(); }
    
    const std::string& sorted_entry_at(std::size_t i);
    
    static void merge(Dictionary& destination, const Dictionary& source);
    
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
    Dictionary dictionary;
    std::vector<hash_type> parse;
    std::vector<size_type> trigger_strings_position; // position of first char of each trigger string
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
    Dictionary* dictionary = nullptr;

    // Shorthands
    hash_type w, p;
    size_type parse_size, tags;
    
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
            csv << "w,p,f,parse_lenght,dict_phrases,dict_tot_length\n";
            csv << params.w << ",";
            csv << params.p << ",";
            csv << statistics.parse_length << ",";
            csv << statistics.num_of_phrases_dictionary << ",";
            csv << statistics.total_dictionary_length << ",";
            csv << "\n";
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
    
    Dictionary dictionary;
    
    hash_type w, p;
    size_type parse_size, tags;
    
    bool closed = false;
    
public:
    
    void init(const Params& params, const std::string& prefix);
    
    ParserFasta(const Params& params, const std::string& file_path, const std::string& out_prefix)
    {
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
            csv << "w,p,f,parse_lenght,dict_phrases,dict_tot_length\n";
            csv << params.w << ",";
            csv << params.p << ",";
            csv << statistics.parse_length << ",";
            csv << statistics.num_of_phrases_dictionary << ",";
            csv << statistics.total_dictionary_length << ",";
            csv << "\n";
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
    
    Dictionary dictionary;
    
    hash_type w, p;
    size_type parse_size, tags;
    
    bool closed = false;

public:
    
    void init(const Params& params, const std::string& prefix);
    
    ParserText(const Params& params, const std::string& file_path, const std::string& out_prefix)
    {
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
            csv << "w,p,f,parse_lenght,dict_phrases,dict_tot_length\n";
            csv << params.w << ",";
            csv << params.p << ",";
            csv << statistics.parse_length << ",";
            csv << statistics.num_of_phrases_dictionary << ",";
            csv << statistics.total_dictionary_length << ",";
            csv << "\n";
            csv.close();
        }
    }
    
    const std::string& get_file_name() const { return this->out_file_name; }
    const Statistics& get_statistics() const { return this->statistics; }
    
    void operator()();
    void close();
};

//------------------------------------------------------------------------------

class ParserUtils
{
public:
    
    static void read_parse(std::string parse_file_name, std::vector<size_type>& parse);
    static void read_dictionary(std::string dic_file_name, std::vector<std::string>& dictionary_vector);
    static void merge(std::string left_prefix, std::string right_prefix, std::string out_prefix, const Params& params);
    static void parse_fasta(std::string fasta_file_name, std::string out_prefix, const Params& params);
    
    static std::vector<std::size_t> compute_occurrences(std::vector<std::string>& dictionary_vector, std::vector<size_type>& parse);
};

} // end namespace pfp
} // end namespace vcfbwt

#endif //pfp_algo_hpp
