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

    bool closed = false;
    bool compress_dictionary = false;
    
    struct d_prime
    {
        std::vector<std::string> d_prime_vector;
        std::vector<std::string> d_prime_vector_additions;

        size_type insert(const std::string& element);
        void remove(size_type i);
        const std::string& at(size_type i) const;
    } D_prime;
    
    int cost_of_removing_trigger_string(const std::string_view& ts);
    
    // Initialize D, P and T_table
    void init_structures();
    
    // Initialize priority queue
    void init_costs();
    bool costs_initialized = false;
    
    // Compression steps
    std::size_t remove_by_cost(std::set<std::string_view>& removed_trigger_strings, int threshold);
    std::size_t remove_simple(std::set<std::string_view>& removed_trigger_strings);
    
    std::size_t initial_size = 0;

public:
    
    // Builds the structures end empties the vectors when done
    AuPair(std::string in, size_type w, bool c, size_type batch_s = 1)
    : window_length(w), in_prefix(std::move(in)), batch_size(batch_s), compress_dictionary(c)
    { this->init_structures(); }
    
    ~AuPair() { this->close(); }

    void close();
    
    std::size_t operator()(std::set<std::string_view>& removed_trigger_strings, int threshold = 0)
    {
        std::size_t removed_bytes = 0;
        
        // First remove simple ts
        removed_bytes += remove_simple(removed_trigger_strings);
        
        // Init costs
        if (not costs_initialized) { costs_initialized = true; init_costs(); }
        
        // Remove ts with cost over threshold
        removed_bytes += remove_by_cost(removed_trigger_strings, threshold);
        
        return removed_bytes;
    }
};

}
}

#endif //au_pair_algo_hpp
