//
//  version.hpp.in
//

#ifndef version_hpp
#define version_hpp

#include <utils.hpp>

namespace vcfbwt
{

struct Version
{
    static std::string VCFBWT_GIT_BRANCH;
    static std::string VCFBWT_GIT_COMMIT_HASH;
    static int VCFBWT_MAJOR;
    static int VCFBWT_MINOR;
    static int VCFBWT_PATCH;
    
    static void print();
};

}

#endif //version_hpp