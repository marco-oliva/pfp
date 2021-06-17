////////////////////////////////////////////////////////////////////////////////
// indexMaxPQ.h
//   Implementation of the maximum-oriented indexed Priority Queue based on
//   Sedgewick and Wayne's Java implementation retrivable at
//   https://algs4.cs.princeton.edu/24pq/IndexMaxPQ.java.html.
////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2019 Zsuzsanna Lipták, Simon J. Puglisi and Massimiliano Rossi
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
////////////////////////////////////////////////////////////////////////////////

#ifndef INDEXMAXPQ_H
#define INDEXMAXPQ_H

#include <vector>
#include <string>
#include <sstream> // to use stringstream
#include <iomanip> // to use stringstream
#include <assert.h>

using namespace std;
// Arguments:
//   X[0..n-1] = input string,
//   C[0..n-1] = input string colors,
//   y = a color,
//   F = a pointer (can to be NULL) to a container storing the output
//     y-unique minimal patterns as a sequence of pairs (pattern, distance)
//     where pattern is a substring (pos,len) of X and distance stores the
//     position of the color y with respect to pos + len -1.
// Returns:
//   the number of minimal pattern in the parsing of X.

class indexMaxPQ{
public:
    // Index Max PQ constructor
    indexMaxPQ():
            n(0),
            nge0(0),
            pq(),
            ipq(),
            keys()
    {
        //Ntd
    }

    // Index Max PQ distructor
    ~indexMaxPQ(){
        // Ntd
    }

    // Init the Index Max PQ
    // Arguments:
    //  n = number of elements in the priority queue
    void init(int _n){
        n = _n;
        nge0 = 0;
        pq.resize(n+1);
        ipq.resize(_n);
        keys.resize(_n,-1);
        for(int i = 1; i < n+1; ++i){
            pq[i] = i-1;
            ipq[pq[i]] = i;
        }
    }

    // Initialize a new key in the PQ
    // Arguments:
    //  index = the index on the key list
    //  key   = the priority of the index
    void push(int index, int key){
        //assert(key < n);
        if(key > keys[index]) promote(index,key);
        else demote(index,key);
    }

    // Increase the priority of index to key
    // Arguments:
    //  index = the index on the key list
    //  key   = the priority of the index
    void promote(int index, int key){
        if(keys[index] < 0 && key >= 0) nge0++;
        keys[index] = key;
        index = ipq[index];
        while( index > 1 && greater(index,index/2)){
            swape(index,index/2);
            index/=2;
        }
    }


    // Decrease the priority of index to key
    // Arguments:
    //  index = the index on the key list
    //  key   = the priority of the index
    void demote(int index, int key){
        if(keys[index] >= 0 && key < 0) nge0--;
        keys[index] = key;
        index = ipq[index];
        while(2*index <= n){
            int tmp = 2*index;
            if(tmp < n && greater(tmp+1,tmp)) tmp++;
            if(greater(index,tmp)) break;
            swape(index,tmp);
            index = tmp;
        }
    }

    // Check whethere there are non negative keys
    bool is_empty(){
        return (nge0 == 0);
    }

    pair<int,int> get_max(){
        return make_pair(keys[pq[1]],pq[1]);
    }

    // Give the key value at that index
    // Arguments:
    //  index = the index in the keys list
    // Returns:
    //  the value of the key at index.
    int get_key(int index){
        return keys[index];
    }

    // Print the priority queue
    // Returns:
    //  the string with the values of the priority queue.
    string print(){
        stringstream ss;
        for(int i = 1; i < pq.size(); ++i){
            ss << setw(2) << pq[i] << " ";
        }
        ss << endl;
        for(int i = 1; i < pq.size(); ++i){
            ss << setw(2) << keys[pq[i]] << " ";
        }
        ss << endl;
        return ss.str();
    }

protected:

    // Compare two elements in the PQ
    // Arguments:
    //  i = the first index in the priority queue
    //  j = the second index in the priority queue
    // Returns:
    //  true   if keys[pq[i]] > keys[pq[j]] OR (keys[pq[i]] == keys[pq[j]] AND pq[i] > pq[j])
    //  false  otherwise
    bool greater(int i,int j){
        return (keys[pq[i]] > keys[pq[j]] || (keys[pq[i]] == keys[pq[j]] && pq[i] > pq[j]));
    }

    // Compare two elements in the PQ
    // Arguments:
    //  i = the first index in the priority queue
    //  j = the second index in the priority queue
    // Returns:
    //  true   if keys[pq[i]] < keys[pq[j]] OR (keys[pq[i]] == keys[pq[j]] AND pq[i] < pq[j])
    //  false  otherwise
    bool smaller(int i,int j){
        return (keys[pq[i]] < keys[pq[j]] || (keys[pq[i]] == keys[pq[j]] && pq[i] < pq[j]));
    }

    // Exchange two elements in the PQ
    // Arguments:
    //  i = the first index in the priority queue
    //  j = the secondirst index in the priority queue
    void swape(int i,int j){
        int tmp = pq[i];
        pq[i] = pq[j];
        pq[j] = tmp;
        ipq[pq[i]] = i;
        ipq[pq[j]] = j;
    }

private:
    int n;               // The number of elements in the priority queue.
    int nge0;            // The number of elements >= 0 in the priority queue.
    vector<int> pq;      // binary heap using 1-based indexing.
    vector<int> ipq;     // inverse of pq. i.e. pq[ipq[i]] = ipq[pq[i]]= i.
    vector<int> keys;    // the keys of the priority queue. i.e. keys[i] = priority of i.

};



#endif /* INDEXMAXPQ_H */
