//
//  au_pair_algo.hpp.h
//
//  Copyright 2021 Marco Oliva. All rights reserved.
//

#ifndef au_pair_algo_hpp
#define au_pair_algo_hpp

#include <utils.hpp>
#include <internals.hpp>
#include <pfp_algo.hpp>

namespace vcfbwt
{
namespace pfp
{

//------------------------------------------------------------------------------

class AuPair
{
private:
    // Table to compute the cost of removing each trigger string.
    std::map<std::string_view, std::vector<size_type*>> T_table;
    
    // Priority queue
    indexMaxPQ priority_queue;
    std::map<std::string_view, size_type> trigger_string_pq_ids;
    std::map<size_type, std::string_view> trigger_string_pq_ids_inv;
    
    // Parse
    LinkedList<size_type> parse;
    
    size_type window_length;
    std::string in_prefix;
    size_type batch_size;
    
    std::size_t curr_id;
    
    bool closed = false;
    
    struct d_prime
    {
        std::vector<std::string> d_prime_vector;
        std::map<size_type, std::string> d_prime_map;
        
        void remove(size_type i);
        const std::string& at(std::size_t i) const;
    } D_prime;
    
    int cost_of_removing_trigger_string(const std::string_view& ts);

public:
    
    // Builds the structures end empties the vectors when done
    AuPair(std::string in, size_type w, size_type batch_s = 1) : window_length(w), in_prefix(std::move(in)), batch_size(batch_s) { this->init(); }
    
    ~AuPair() { this->close(); }
    
    void init();
    void close();
    
    std::size_t compress(std::set<std::string_view>& removed_trigger_strings, int threshold);
    std::size_t remove_simple(std::set<std::string_view>& removed_trigger_strings);
};

}
}

#endif //au_pair_algo_hpp