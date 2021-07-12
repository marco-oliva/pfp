//
// Created by marco on 7/9/21.
//

#ifndef VCF_BWT_INTERNALS_HPP
#define VCF_BWT_INTERNALS_HPP

#include <utils.hpp>


namespace vcfbwt
{
namespace pfp
{

template <typename DataType> // intended to be an integer
class LinkedList
{
private:

    // Data array
    std::vector<DataType> data;

    // Out of Range
    DataType out_of_range;

    // Start pointer
    long_type start_pointer;

    // Bitvector to mark deleted elements
    std::vector<bool> deleted_elements;

    long_type convert(const DataType* ptr)
    {
        DataType* start = &(data[0]);
        return ptr - start;
    }
    
    bool initialized = false;
    long_type num_of_elements = 0;

public:

    LinkedList(long_type size) :
        data(size), deleted_elements(size, false), start_pointer(0), initialized(true),
        num_of_elements(size) {}
    
    LinkedList() : data(0), deleted_elements(0), start_pointer(0), initialized(false) {}
    
    void init(const DataType* in, long_type size)
    {
        initialized = true;
        deleted_elements.resize(size, false);
        num_of_elements = size;
        
        data.resize(size);
        std::memcpy((char*) &(data[0]), (char*) in, size * sizeof(DataType));
    }
    
    DataType& operator[](long_type i) { assert(not deleted_elements[i]); return data[i]; }
    DataType& at(long_type i) { assert(not deleted_elements[i]); return data.at(i); }

    DataType* next_at(long_type i)
    {
        assert(not deleted_elements[i]);

        if (deleted_elements[i + 1])
        {
            if (data[i + 1] == 0) { return &out_of_range; }
            else { return &(data[i + data[i + 1] + 1]); }
        }
        else { return &(data[i + 1]); }
    }


    DataType* prev_at(long_type i)
    {
        assert(not deleted_elements[i]);

        if (deleted_elements[i - 1])
        {
            if (data[i - 1] == 0) { return &out_of_range; }
            else { return &(data[i - data[i - 1] - 1]); }
        }
        else { return &(data[i - 1]); }
    }


    void remove_at(long_type i)
    {
        assert(not deleted_elements[i]);

        deleted_elements[i] = true;
        num_of_elements -= 1;

        // first element
        if (i == start_pointer)
        {
            if (deleted_elements[i + 1])
            {
                start_pointer = i + data[i + 1] + 1;
                data[i + data[i + 1]] = 0;
                data[i + 1] = 0;
            }
            else { start_pointer = start_pointer + 1; }
            data[i] = 0;
            return;
        }
        if (i == (data.size() - 1))
        {
            if (deleted_elements[i - 1]) { data[i - data[i - 1]] = 0; data[i - 1] = 0; }
            data[i] = 0;
            return;
        }

        // check right
        if (deleted_elements[i + 1] and (not deleted_elements[i - 1]))
        {
            data[i] = data[i + 1] + 1;
            data[i + data[i] - 1] = data[i];
            if (data[i] > 2) { data[i + 1] = 0; }
        }

        // check left
        if (deleted_elements[i - 1] and (not deleted_elements[i + 1]))
        {
            data[i] = data[i - 1] + 1;
            data[i - data[i] + 1] = data[i];
            if (data[i] > 2) { data[i - 1] = 0; }
        }

        // check left and right
        if (deleted_elements[i - 1] and deleted_elements[i + 1])
        {
            DataType gap = data[i - 1] + data[i + 1] + 1;

            // update left
            data[i - data[i - 1]] = gap;
            if (data[i - 1] != gap) { data[i - 1] = 0; }

            // update right
            data[i + data[i + 1]] = gap;
            if (data[i + 1] != gap) { data[i + 1] = 0; }

            // update here
            data[i] = 0;
        }

        // remove element
        if (not (deleted_elements[i + 1] or deleted_elements[i - 1]))
        {
            data[i] = 1;
        }
    }
    
    DataType* next(DataType* d) { return this->next_at(convert(d)); }
    DataType* prev(DataType* d) { return this->prev_at(convert(d)); }
    void remove(DataType* d) { this->remove(convert_at(d)); }
    long_type size() const { return this->num_of_elements; }

    DataType* end() { return &(out_of_range); }
    DataType* begin() { return &(data[start_pointer]); }

};

} // end namespace pfp
} // end namespace vcfbwt

#endif //VCF_BWT_INTERNALS_HPP
