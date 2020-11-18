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
    DOLLAR_PRIME = 3
};

//------------------------------------------------------------------------------

class Dictionary
{
private:
    
    bool sorted = false;
    
    struct DictionaryEntry
    {
        size_type occurrences = 1;
        std::string phrase;
    
        DictionaryEntry() = default;
        DictionaryEntry(const std::string& s) : occurrences(1), phrase(s) {}
    };
    
    std::unordered_map<hash_type, DictionaryEntry> hash_string_map;
    std::vector<std::pair<std::reference_wrapper<std::string>, hash_type>> sorted_phrases;
    std::unordered_map<hash_type, size_type> hash_to_ranks;
    
    void sort();
    
public:
    
    Dictionary() = default;
    
    hash_type add(const std::string& phrase);
    hash_type get(const std::string& phrase) const;
    bool contains(const std::string& phrase) const;
    void update_frequency(const std::string& phrase);

    size_type hash_to_rank(hash_type hash);
    
    size_type size() { return hash_string_map.size(); }
    
    const std::string& sorted_entry_at(std::size_t i) {if (not this->sorted) { sort(); } return sorted_phrases[i].first.get(); }
    
    static void merge(Dictionary& destination, const Dictionary& source);
    
    friend class Parser;
};

//------------------------------------------------------------------------------

struct Params
{
    hash_type p = 100;
    hash_type w =  10;
    bool compute_seeded_trigger_strings = true;
    bool compress_dictionary = false;
    bool use_acceleration = false;
    bool print_out_statistics_csv = false;
    
    inline std::size_t min_seed_region_length() const { return w * 3; }
    double min_frequency = 0.001;
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
    
    const Params& params;
    
    void init(const std::string& reference, const std::unordered_set<hash_type>& ts);
    
    ReferenceParse(const std::string& reference, const std::unordered_set<hash_type>& ts, const Params& pms) : params(pms) { this->init(reference,ts); }
    const hash_type& operator[](size_type i) const { return this->parse[i]; }
};

// Create a parse on disk
class Parser
{

private:
    
    std::ofstream out_file;
    std::string out_file_prefix;
    std::string out_file_name;
    std::string tmp_out_file_name;
    
    Params params;
    Statistics statistics;
    
    ReferenceParse* reference_parse = nullptr;
    
    Dictionary dictionary;
    hash_type w, p;
    size_type parse_size, tags;
    
    bool closed = false;
    
    std::vector<std::reference_wrapper<Parser>> registered_workers;
    
public:
    
    enum tags
    {
        MAIN = 1,
        WORKER = 2,
        COMPRESSED = 4,
        UNCOMPRESSED = 8,
        LAST = 16
    };
    
    void init(const Params& params, const std::string& prefix, ReferenceParse& rp, size_type t = MAIN | UNCOMPRESSED);
    
    Parser(const Params& params, const std::string& file_path, ReferenceParse& rp, size_type t = MAIN | UNCOMPRESSED)
    {
        this->init(params, file_path, rp, t);
    }
    
    Parser() = default;
    
    ~Parser()
    {
        close();
        size_type total_length = 0;
        for (auto& entry : dictionary.hash_string_map) { total_length += entry.second.phrase.size(); }
    
        // Fill out statistics
        this->statistics.parse_length = this->parse_size;
        this->statistics.num_of_phrases_dictionary = this->dictionary.sorted_phrases.size();
        this->statistics.total_dictionary_length = total_length;
        std::string name;
        if (tags & WORKER) { name = "Worker"; } else { name = "Main"; }
        spdlog::info("{} -\tParse size: {}\tDic Size: {} Dic Total Length: {}", name ,parse_size, dictionary.size(), total_length);
        
        if (params.print_out_statistics_csv and (tags & WORKER))
        {
            std::ofstream csv(out_file_prefix + ".csv");
            csv << "w,p,f,parse_lenght,dict_phrases,dict_tot_length\n";
            csv << params.w << ",";
            csv << params.p << ",";
            csv << params.min_frequency << ",";
            csv << statistics.parse_length << ",";
            csv << statistics.num_of_phrases_dictionary << ",";
            csv << statistics.total_dictionary_length << ",";
            csv << "\n";
            csv.close();
        }
    }
    
    void close();
    const std::string& get_file_name() const { return this->out_file_name; }
    const Statistics& get_statistics() const { return this->statistics; }
    void register_worker(Parser& parser) { this->registered_workers.push_back(std::ref(parser)); }
    
    void operator()(const Sample& sample, const std::unordered_set<hash_type>& trigger_strings);
    
    static void compute_trigger_strings(VCF& vcf, const Params& params, std::unordered_set<hash_type>& trigger_string_set);
    static void read_parse(std::string parse_file_name, std::vector<size_type>& parse);
    static void read_dictionary(std::string dic_file_name, std::vector<std::string>& dictionary_vector);
    static void merge(std::string left_prefix, std::string right_prefix, std::string out_prefix, const Params& params);
    static void parse_fasta(std::string fasta_file_name, std::string out_prefix, const std::unordered_set<hash_type>& trigger_strings, const Params& params);
};

//------------------------------------------------------------------------------

} // end namespace pfp
} // end namespace vcfbwt

#endif //pfp_algo_hpp
