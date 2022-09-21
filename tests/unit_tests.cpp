//
//  unit_tests.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>

#include <vcf.hpp>
#include <utils.hpp>
#include <pfp_algo.hpp>

//------------------------------------------------------------------------------

struct listener : Catch::EventListenerBase
{
    using EventListenerBase::EventListenerBase;

    virtual void testCaseStarting(Catch::TestCaseInfo const& testInfo) override
    {
        std::cout << testInfo.tagsAsString() << " " << testInfo.name << std::endl;
    }
};
CATCH_REGISTER_LISTENER(listener)

//------------------------------------------------------------------------------

std::string testfiles_dir = "../tests/files";

std::size_t p_global = 75;
std::size_t w_global = 20;

//------------------------------------------------------------------------------
template <typename data_type>
bool
unparse_and_check(std::string& in_prefix, std::vector<data_type>& what_it_should_be, std::size_t window_length, char DOLLAR, bool n = false)
{
    // Unparse
    std::vector<vcfbwt::size_type> parse;
    std::string parse_ext = n ? vcfbwt::EXT::N_PARSE : vcfbwt::EXT::PARSE;
    std::string dictionary_ext = n ? vcfbwt::EXT::N_DICT : vcfbwt::EXT::DICT;
    vcfbwt::pfp::ParserUtils<data_type>::read_parse(in_prefix + parse_ext, parse);
    std::vector<std::vector<data_type>> dictionary;
    vcfbwt::pfp::ParserUtils<data_type>::read_dictionary(in_prefix + dictionary_ext, dictionary);
    
    std::vector<data_type> unparsed;
    std::vector<vcfbwt::long_type> occ_computed(dictionary.size(), 0);
    for (auto& p : parse)
    {
        occ_computed[p - 1] += 1;

        if (p > dictionary.size()) { spdlog::error("Something wrong in the parse"); exit(EXIT_FAILURE); }
        std::vector<data_type> dict_string(dictionary[p - 1].begin(), dictionary[p - 1].begin() + dictionary[p - 1].size() - window_length);
        unparsed.insert(unparsed.end(), dict_string.begin(), dict_string.end());
    }
    unparsed.insert(unparsed.end(), window_length, DOLLAR);


    // Check the occ file
    bool occ_good = true;
    std::string occ_ext = n ? vcfbwt::EXT::N_OCC : vcfbwt::EXT::OCC;
    std::ifstream occ_stream(in_prefix + occ_ext);
    if (parse.size() < std::numeric_limits<vcfbwt::short_type>::max())
    {
        std::vector<vcfbwt::short_type> occ(dictionary.size(), 0);
        occ_stream.read((char*) occ.data(), sizeof(vcfbwt::short_type) * occ.size());
        int should_be_eof = occ_stream.get(); occ_good = occ_good and occ_stream.eof();

        for (std::size_t i = 0; i < occ.size(); i++) { occ_good = occ_good and (occ[i] == occ_computed[i]); }
    }
    else
    {
        std::vector<vcfbwt::long_type> occ(dictionary.size(), 0);
        occ_stream.read((char*) occ.data(), sizeof(vcfbwt::long_type) * occ.size());
        int should_be_eof = occ_stream.get(); occ_good = occ_good and occ_stream.eof();

        for (std::size_t i = 0; i < occ.size(); i++) { occ_good = occ_good and (occ[i] == occ_computed[i]); }
    }


    // Compare the two strings
    std::size_t i = 0;
    while ( ((i < unparsed.size()) and (i < what_it_should_be.size()))
    and (unparsed[i] == what_it_should_be[i])) { i++; }
    spdlog::info("First missmatch: {}", i);
    return (occ_good and ((i == (unparsed.size())) and (i == (what_it_should_be.size()))));
}

//------------------------------------------------------------------------------
TEST_CASE( "Initialization", "[LinkedList]" )
{
    vcfbwt::pfp::LinkedList<vcfbwt::size_type> linked_list(10);

    for (vcfbwt::size_type i = 0; i < 10; i++) { linked_list[i] = i; }

    bool all_set = true;
    for (vcfbwt::size_type i = 0; i < 10; i++) { all_set = all_set and (linked_list[i] == i); }

    REQUIRE(all_set);
}

TEST_CASE( "Delete from the middle", "[LinkedList]" )
{
    vcfbwt::pfp::LinkedList<vcfbwt::size_type> linked_list(10);

    for (vcfbwt::size_type i = 0; i < 10; i++) { linked_list[i] = i + 100; }

    linked_list.remove_at(4);
    linked_list.remove_at(5);
    linked_list.remove_at(6);
    linked_list.remove_at(3);

    REQUIRE(*(linked_list.next_at(2)) == 107);
}

TEST_CASE( "Meet deletions", "[LinkedList]" )
{
    vcfbwt::pfp::LinkedList<vcfbwt::size_type> linked_list(10);

    for (vcfbwt::size_type i = 0; i < 10; i++) { linked_list[i] = i + 100; }

    linked_list.remove_at(7);
    linked_list.remove_at(8);
    linked_list.remove_at(5);
    linked_list.remove_at(6);

    REQUIRE(*(linked_list.next_at(4)) == 109);
}

TEST_CASE( "Remove first 3 left to right", "[LinkedList]" )
{
    vcfbwt::pfp::LinkedList<vcfbwt::size_type> linked_list(10);

    for (vcfbwt::size_type i = 0; i < 10; i++) { linked_list[i] = i + 100; }

    linked_list.remove_at(0);
    linked_list.remove_at(1);
    linked_list.remove_at(2);

    REQUIRE(*(linked_list.begin()) == 103);
}

TEST_CASE( "Remove first 3 right to left", "[LinkedList]" )
{
    vcfbwt::pfp::LinkedList<vcfbwt::size_type> linked_list(10);

    for (vcfbwt::size_type i = 0; i < 10; i++) { linked_list[i] = i + 100; }

    linked_list.remove_at(2);
    linked_list.remove_at(1);
    linked_list.remove_at(0);

    REQUIRE(*(linked_list.begin()) == 103);
}

TEST_CASE( "Remove first 3, 021", "[LinkedList]" )
{
    vcfbwt::pfp::LinkedList<vcfbwt::size_type> linked_list(10);

    for (vcfbwt::size_type i = 0; i < 10; i++) { linked_list[i] = i + 100; }

    linked_list.remove_at(0);
    linked_list.remove_at(2);
    linked_list.remove_at(1);

    REQUIRE(*(linked_list.begin()) == 103);
}

TEST_CASE( "Remove last 3 in order, left to right", "[LinkedList]" )
{
    vcfbwt::pfp::LinkedList<vcfbwt::size_type> linked_list(10);

    for (vcfbwt::size_type i = 0; i < 10; i++) { linked_list[i] = i + 100; }

    linked_list.remove_at(7);
    linked_list.remove_at(8);
    linked_list.remove_at(9);

    REQUIRE(linked_list.next_at(6) == linked_list.end());
}

TEST_CASE( "Remove last 3 in order, right to left", "[LinkedList]" )
{
    vcfbwt::pfp::LinkedList<vcfbwt::size_type> linked_list(10);

    for (vcfbwt::size_type i = 0; i < 10; i++) { linked_list[i] = i + 100; }

    linked_list.remove_at(9);
    linked_list.remove_at(8);
    linked_list.remove_at(7);

    REQUIRE(linked_list.next_at(6) == linked_list.end());
}

TEST_CASE( "Meet deletions at the end", "[LinkedList]" )
{
    vcfbwt::pfp::LinkedList<vcfbwt::size_type> linked_list(10);

    for (vcfbwt::size_type i = 0; i < 10; i++) { linked_list[i] = i + 100; }

    linked_list.remove_at(5);
    linked_list.remove_at(6);
    linked_list.remove_at(9);
    linked_list.remove_at(7);
    linked_list.remove_at(8);

    REQUIRE(linked_list.next_at(4) == linked_list.end());
}

TEST_CASE( "Meet deletions and prev,next", "[LinkedList]" )
{
    vcfbwt::pfp::LinkedList<vcfbwt::size_type> linked_list(10);

    for (vcfbwt::size_type i = 0; i < 10; i++) { linked_list[i] = i + 100; }

    linked_list.remove_at(7);
    linked_list.remove_at(8);
    linked_list.remove_at(5);
    linked_list.remove_at(6);

    REQUIRE(*(linked_list.next_at(4)) == 109);
    REQUIRE(*(linked_list.prev(linked_list.next_at(4))) == 104);
}

TEST_CASE( "Remove last and meet", "[LinkedList]" )
{
    vcfbwt::pfp::LinkedList<vcfbwt::size_type> linked_list(10);

    for (vcfbwt::size_type i = 0; i < 10; i++) { linked_list[i] = i + 100; }

    linked_list.remove_at(7);
    linked_list.remove_at(8);
    linked_list.remove_at(5);
    linked_list.remove_at(6);
    linked_list.remove_at(9);

    REQUIRE(linked_list.next_at(4) == linked_list.end());
}

TEST_CASE( "Remove first and meet", "[LinkedList]" )
{
    vcfbwt::pfp::LinkedList<vcfbwt::size_type> linked_list(10);

    for (vcfbwt::size_type i = 0; i < 10; i++) { linked_list[i] = i + 100; }

    linked_list.remove_at(3);
    linked_list.remove_at(1);
    linked_list.remove_at(2);
    linked_list.remove_at(0);

    REQUIRE(*(linked_list.begin()) == 104);
    REQUIRE(*(linked_list.next(linked_list.begin())) == 105);
}

TEST_CASE( "Remove last but don't meet", "[LinkedList]" )
{
    vcfbwt::pfp::LinkedList<vcfbwt::size_type> linked_list(10);

    for (vcfbwt::size_type i = 0; i < 10; i++) { linked_list[i] = i + 100; }

    linked_list.remove_at(7);
    linked_list.remove_at(5);
    linked_list.remove_at(6);
    linked_list.remove_at(9);

    REQUIRE(*(linked_list.next_at(4)) == 108);
    REQUIRE(linked_list.next(linked_list.next_at(4)) == linked_list.end());
}

TEST_CASE( "Remove first but don't meet", "[LinkedList]" )
{
    vcfbwt::pfp::LinkedList<vcfbwt::size_type> linked_list(10);

    for (vcfbwt::size_type i = 0; i < 10; i++) { linked_list[i] = i + 100; }

    linked_list.remove_at(3);
    linked_list.remove_at(2);
    linked_list.remove_at(0);

    REQUIRE(*(linked_list.begin()) == 101);
    REQUIRE(*(linked_list.next(linked_list.begin())) == 104);
}

TEST_CASE( "Remove everything but first and last", "[LinkedList]" )
{
    vcfbwt::pfp::LinkedList<vcfbwt::size_type> linked_list(10);

    for (vcfbwt::size_type i = 0; i < 10; i++) { linked_list[i] = i + 100; }

    linked_list.remove_at(7);
    linked_list.remove_at(2);
    linked_list.remove_at(5);
    linked_list.remove_at(6);
    linked_list.remove_at(3);
    linked_list.remove_at(1);
    linked_list.remove_at(4);
    linked_list.remove_at(8);

    REQUIRE(*(linked_list.begin()) == 100);
    REQUIRE(*(linked_list.next(linked_list.begin())) == 109);
}

TEST_CASE( "Remove everything but second and second to last", "[LinkedList]" )
{
    vcfbwt::pfp::LinkedList<vcfbwt::size_type> linked_list(10);

    for (vcfbwt::size_type i = 0; i < 10; i++) { linked_list[i] = i + 100; }

    linked_list.remove_at(7);
    linked_list.remove_at(2);
    linked_list.remove_at(6);
    linked_list.remove_at(3);
    linked_list.remove_at(0);
    linked_list.remove_at(5);
    linked_list.remove_at(4);
    linked_list.remove_at(9);

    REQUIRE(*(linked_list.begin()) == 101);
    REQUIRE(*(linked_list.next(linked_list.begin())) == 108);
}

TEST_CASE( "Remove everything but third and third to last", "[LinkedList]" )
{
    vcfbwt::pfp::LinkedList<vcfbwt::size_type> linked_list(10);

    for (vcfbwt::size_type i = 0; i < 10; i++) { linked_list[i] = i + 100; }

    linked_list.remove_at(8);
    linked_list.remove_at(2);
    linked_list.remove_at(6);
    linked_list.remove_at(1);
    linked_list.remove_at(0);
    linked_list.remove_at(5);
    linked_list.remove_at(4);
    linked_list.remove_at(9);

    REQUIRE(*(linked_list.begin()) == 103);
    REQUIRE(*(linked_list.next(linked_list.begin())) == 107);
}

//------------------------------------------------------------------------------
TEST_CASE( "Initialization KR", "[KR Window]" )
{
    std::string test_string = "12345";

    vcfbwt::KarpRabinHash kr_window(5);
    kr_window.set_constant(256); kr_window.set_prime(2147483647);
    kr_window.initialize(test_string);

    REQUIRE(kr_window.get_hash() == 842216599);
}

TEST_CASE( "Update 1 charachter", "[KR Window]" )
{
    std::string test_string = "12345";

    vcfbwt::KarpRabinHash kr_window(5);
    kr_window.set_constant(256); kr_window.set_prime(2147483647);
    kr_window.initialize(test_string);

    kr_window.update('1', '6');

    REQUIRE(kr_window.get_hash() == 859059610);
}

TEST_CASE( "Update 2 charachter", "[KR Window]" )
{
    std::string test_string = "12345";

    vcfbwt::KarpRabinHash kr_window(5);
    kr_window.set_constant(256); kr_window.set_prime(2147483647);
    kr_window.initialize(test_string);

    kr_window.update('1', '6');
    kr_window.update('2', '7');

    REQUIRE(kr_window.get_hash() == 875902621);
}

TEST_CASE( "Periodic string", "[KR Window]" )
{
    std::string test_string = "11111";

    vcfbwt::KarpRabinHash kr_window(5);
    kr_window.set_constant(256); kr_window.set_prime(2147483647);
    kr_window.initialize(test_string);

    vcfbwt::hash_type before = kr_window.get_hash();
    kr_window.update('1', '1');
    vcfbwt::hash_type after = kr_window.get_hash();

    REQUIRE(before == 825307539);
    REQUIRE(after  == 825307539);
}

TEST_CASE( "Periodic string of Ns", "[KR Window]" )
{
    std::string test_string = "NNNNNNNNNNNNNNNNNNNN";

    vcfbwt::KarpRabinHash kr_window(test_string.size());
    kr_window.set_constant(256); kr_window.set_prime(2147483647);
    kr_window.initialize(test_string);

    vcfbwt::hash_type before = kr_window.get_hash();
    kr_window.update('N', 'N');
    vcfbwt::hash_type after = kr_window.get_hash();

    REQUIRE(before == 2071690116);
    REQUIRE(after  == 2071690116);
}

TEST_CASE( "String Hash", "[KR Window]" )
{
    std::string test_string = "12345";

    vcfbwt::KarpRabinHash kr_window(5);
    kr_window.initialize(test_string);

    REQUIRE(kr_window.get_hash() == vcfbwt::KarpRabinHash::string_hash(test_string));
}

TEST_CASE( "String Hash 2", "[KR Window]" )
{
    std::string test_string = "12345";

    vcfbwt::KarpRabinHash kr_window(5);
    kr_window.initialize(test_string);
    kr_window.update('1', '6');

    REQUIRE(kr_window.get_hash() == vcfbwt::KarpRabinHash::string_hash("23456"));
}

TEST_CASE( "Reproducing bug", "[KR Window]" )
{
    std::string test_string_1 = ".a.little.late;You.fo";
    std::string test_string_2 = "\5\5\5\5\4You.fo";

    vcfbwt::KarpRabinHash kr_window_1(5, true);
    kr_window_1.initialize(test_string_1.substr(0,5));
    for (std::size_t i = 0; i <= 15; i++) { kr_window_1.update(test_string_1[i], test_string_1[i + 5]); }

    vcfbwt::KarpRabinHash kr_window_2(5, true);
    kr_window_2.initialize(test_string_2.substr(0,5));
    for (std::size_t i = 0; i <= 5; i++) { kr_window_2.update(test_string_2[i], test_string_2[i + 5]); }

    REQUIRE(kr_window_1.get_hash() == vcfbwt::KarpRabinHash::string_hash("ou.fo"));
    REQUIRE(kr_window_2.get_hash() == vcfbwt::KarpRabinHash::string_hash("ou.fo"));
    REQUIRE(kr_window_1.get_hash() % 30 == 0);
    REQUIRE(kr_window_2.get_hash() % 30 == 0);
}

//------------------------------------------------------------------------------

TEST_CASE( "Initialization KR Mersenne", "[KR Mersenne Window]" )
{
    std::string test_string = "12345";

    vcfbwt::Mersenne_KarpRabinHash kr_window(5);
    // base = 660162925935593667, prime = 2305843009213693951
    kr_window.initialize(test_string);

    REQUIRE(kr_window.get_hash() == 1337084880462018802);
}

TEST_CASE( "Update 1 charachter Mersenne", "[KR Mersenne Window]" )
{
    std::string test_string = "12345";

    vcfbwt::Mersenne_KarpRabinHash kr_window(5);
    // base = 660162925935593667, prime = 2305843009213693951
    kr_window.initialize(test_string);

    kr_window.update('1', '6');

    REQUIRE(kr_window.get_hash() == 255739482431644834);
}

TEST_CASE( "Update 2 charachters Mersenne", "[KR Mersenne Window]" )
{
    std::string test_string = "12345";

    vcfbwt::Mersenne_KarpRabinHash kr_window(5);
    // base = 660162925935593667, prime = 2305843009213693951
    kr_window.initialize(test_string);

    kr_window.update('1', '6');
    kr_window.update('2', '7');

    REQUIRE(kr_window.get_hash() == 1480237093614964817);
}

TEST_CASE( "Periodic string Mersenne", "[KR Mersenne Window]" )
{
    std::string test_string = "11111";

    vcfbwt::Mersenne_KarpRabinHash kr_window(5);
    // base = 660162925935593667, prime = 2305843009213693951
    kr_window.initialize(test_string);

    vcfbwt::hash_type before = kr_window.get_hash();
    kr_window.update('1', '1');
    vcfbwt::hash_type after = kr_window.get_hash();

    REQUIRE(before == 48464708426636441);
    REQUIRE(after  == 48464708426636441);
}

//------------------------------------------------------------------------------

TEST_CASE( "Initialization KR Mersenne4", "[KR Mersenne4 Window]" )
{
    std::vector<int32_t> test_data = {1,2,3,4,5};

    vcfbwt::Mersenne_KarpRabinHash4 kr_window(5 * 4, true);
    // base = 660162925935593667, prime = 2305843009213693951
    kr_window.initialize(std::string_view((char*) test_data.data(), test_data.size() * sizeof(int32_t)));

    REQUIRE(kr_window.get_hash() == 207274774005008393);
}

TEST_CASE( "Update 1 charachter Mersenne4", "[KR Mersenne4 Window]" )
{
    std::vector<int32_t> test_data = {1,2,3,4,5};

    vcfbwt::Mersenne_KarpRabinHash4 kr_window(5 * 4, true);
    // base = 660162925935593667, prime = 2305843009213693951
    kr_window.initialize(std::string_view((char*) test_data.data(), test_data.size() * sizeof(int32_t)));

    int32_t out = 1, in = 6;
    kr_window.update((const vcfbwt::char_type*) &out, (const vcfbwt::char_type*) &in);

    REQUIRE(kr_window.get_hash() == 1431772385188328376);
}

TEST_CASE( "Update 2 charachters Mersenne4", "[KR Mersenne4 Window]" )
{
    std::vector<int32_t> test_data = {1,2,3,4,5};

    vcfbwt::Mersenne_KarpRabinHash4 kr_window(5 * 4, true);
    // base = 660162925935593667, prime = 2305843009213693951
    kr_window.initialize(std::string_view((char*) test_data.data(), test_data.size() * sizeof(int32_t)));

    int32_t out = 1, in = 6;
    kr_window.update((const vcfbwt::char_type*) &out, (const vcfbwt::char_type*) &in);

    out = 2, in = 7;
    kr_window.update((const vcfbwt::char_type*) &out, (const vcfbwt::char_type*) &in);

    REQUIRE(kr_window.get_hash() == 350426987157954408);
}

TEST_CASE( "Periodic string Mersenne4", "[KR Mersenne4 Window]" )
{
    std::vector<int32_t> test_data = {1,1,1,1,1};

    vcfbwt::Mersenne_KarpRabinHash4 kr_window(5 * 4, true);
    // base = 660162925935593667, prime = 2305843009213693951
    kr_window.initialize(std::string_view((char*) test_data.data(), test_data.size() * sizeof(int32_t)));

    int32_t out = 1, in = 1;
    vcfbwt::hash_type before = kr_window.get_hash();
    kr_window.update((const vcfbwt::char_type*) &out, (const vcfbwt::char_type*) &in);
    vcfbwt::hash_type after = kr_window.get_hash();

    REQUIRE(before == 1224497611183319983);
    REQUIRE(after  == 1224497611183319983);
}

//------------------------------------------------------------------------------

TEST_CASE( "Dictionary size", "[Dictionary]")
{
    vcfbwt::pfp::Dictionary<vcfbwt::char_type> dictionary;

    vcfbwt::size_type elem, tot_elem = 100000;
    for (elem = 0; elem < tot_elem; elem++)
    {
        std::string elem_string = std::to_string(elem);
        dictionary.check_and_add(std::vector<vcfbwt::char_type>(elem_string.begin(), elem_string.end()));
    }

    bool all_elements_in_dict = true;
    for (elem = 0; elem < tot_elem; elem++)
    {
        std::string elem_string = std::to_string(elem);
        all_elements_in_dict + all_elements_in_dict and dictionary.contains(std::vector<vcfbwt::char_type>(elem_string.begin(), elem_string.end()));
    }
    REQUIRE(all_elements_in_dict);
    REQUIRE(dictionary.size() == elem);
}

//------------------------------------------------------------------------------

TEST_CASE( "Constructor with samples specified", "[VCF parser]" )
{
    std::string vcf_file_name = testfiles_dir + "/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz";
    std::string ref_file_name = testfiles_dir + "/Y.fa.gz";
    vcfbwt::VCF vcf(ref_file_name, vcf_file_name, "");

    // read_samples list from file
    std::ifstream samples_file(testfiles_dir + "/samples_list.txt");
    std::vector<std::string> samples_list;
    std::copy(std::istream_iterator<std::string>(samples_file),
              std::istream_iterator<std::string>(),
              std::back_inserter(samples_list));

    bool all_match = true;
    for (std::size_t i = 0; i < vcf.size(); i++)
    {
        all_match = all_match  & (vcf[i].id() == samples_list[i]);
    }

    REQUIRE(vcf.size() == samples_list.size());
    REQUIRE(all_match);
}

TEST_CASE( "Constructor", "[VCF parser]" )
{
    std::string vcf_file_name = testfiles_dir + "/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz";
    std::string ref_file_name = testfiles_dir + "/Y.fa.gz";
    vcfbwt::VCF vcf(ref_file_name, vcf_file_name, "");

    // read_samples list from file
    std::ifstream samples_file(testfiles_dir + "/samples_list.txt");
    std::vector<std::string> samples_list;
    std::copy(std::istream_iterator<std::string>(samples_file),
              std::istream_iterator<std::string>(),
              std::back_inserter(samples_list));

    bool all_match = true;
    for (std::size_t i = 0; i < vcf.size(); i++)
    {
        all_match = all_match  & (vcf[i].id() == samples_list[i]);
    }

    REQUIRE(vcf.size() == samples_list.size());
    REQUIRE(all_match);
}

TEST_CASE("Sample: HG00101", "[VCF parser]")
{
    std::string vcf_file_name = testfiles_dir + "/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz";
    std::string ref_file_name = testfiles_dir + "/Y.fa.gz";
    vcfbwt::VCF vcf(ref_file_name, vcf_file_name, "", 2);

    REQUIRE(vcf[1].id() == "HG00101");

    std::string test_sample_path = testfiles_dir + "/HG00101_chrY_H1.fa.gz";
    std::ifstream in_stream(test_sample_path);

    REQUIRE(vcfbwt::is_gzipped(in_stream));

    zstr::istream is(in_stream);
    std::string line, from_fasta;
    while (getline(is, line))
    {
        if (not(line.empty() or line[0] == '>'))
        {
            from_fasta.append(line);
        }
    }

    vcfbwt::Sample::iterator it(vcf[1]);
    std::string from_vcf;
    while (not it.end())
    {
        from_vcf.push_back(*it);
        ++it;
    }

    std::size_t i = 0;
    while (((i < from_vcf.size()) and (i < from_fasta.size())) and (from_vcf[i] == from_fasta[i]))
    {
        i++;
    }
    REQUIRE(((i == (from_vcf.size())) and (i == (from_fasta.size()))));
}

TEST_CASE( "Sample: HG00103", "[VCF parser]" )
{
    std::string vcf_file_name = testfiles_dir + "/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz";
    std::string ref_file_name = testfiles_dir + "/Y.fa.gz";
    vcfbwt::VCF vcf(ref_file_name, vcf_file_name, "", 3);

    REQUIRE(vcf[2].id() == "HG00103");

    std::string test_sample_path = testfiles_dir + "/HG00103_chrY_H1.fa.gz";
    std::ifstream in_stream(test_sample_path);

    REQUIRE(vcfbwt::is_gzipped(in_stream));

    zstr::istream is(in_stream);
    std::string line, from_fasta;
    while (getline(is, line)) { if ( not (line.empty() or line[0] == '>') ) { from_fasta.append(line); } }

    vcfbwt::Sample::iterator it(vcf[2]);
    std::string from_vcf;
    while (not it.end()) { from_vcf.push_back(*it); ++it; }

    std::size_t i = 0;
    while ( ((i < from_vcf.size()) and (i < from_fasta.size()))
    and (from_vcf[i] == from_fasta[i])) { i++; }
    REQUIRE(((i == (from_vcf.size())) and (i == (from_fasta.size()))));
}

TEST_CASE( "Selecting only Sample: HG00103", "[VCF parser]" )
{
    std::string vcf_file_name = testfiles_dir + "/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz";
    std::string ref_file_name = testfiles_dir + "/Y.fa.gz";
    std::string samples_file_name = testfiles_dir + "/allowed_samples_list.txt";
    vcfbwt::VCF vcf(ref_file_name, vcf_file_name, samples_file_name);

    REQUIRE(vcf[0].id() == "HG00103");

    std::string test_sample_path = testfiles_dir + "/HG00103_chrY_H1.fa.gz";
    std::ifstream in_stream(test_sample_path);

    REQUIRE(vcfbwt::is_gzipped(in_stream));

    zstr::istream is(in_stream);
    std::string line, from_fasta;
    while (getline(is, line)) { if ( not (line.empty() or line[0] == '>') ) { from_fasta.append(line); } }

    vcfbwt::Sample::iterator it(vcf[0]);
    std::string from_vcf;
    while (not it.end()) { from_vcf.push_back(*it); ++it; }

    std::size_t i = 0;
    while ( ((i < from_vcf.size()) and (i < from_fasta.size()))
            and (from_vcf[i] == from_fasta[i])) { i++; }
    REQUIRE(((i == (from_vcf.size())) and (i == (from_fasta.size()))));
}

TEST_CASE( "Reference + Sample HG00096, No acceleration", "[PFP algorithm]" )
{
    std::string vcf_file_name = testfiles_dir + "/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz";
    std::string ref_file_name = testfiles_dir + "/Y.fa.gz";
    vcfbwt::VCF vcf(ref_file_name, vcf_file_name, "", 1);

    // Produce dictionary and parsing
    vcfbwt::pfp::Params params;
    params.w = w_global; params.p = p_global;
    params.use_acceleration = false;
    params.compute_occurrences = true;
    vcfbwt::pfp::ReferenceParse reference_parse(vcf.get_reference(), params);

    std::string out_prefix = testfiles_dir + "/parser_out";
    vcfbwt::pfp::ParserVCF main_parser(params, out_prefix, reference_parse);

    vcfbwt::pfp::ParserVCF worker;
    std::size_t tag = 0;
    tag = tag | vcfbwt::pfp::ParserVCF::WORKER;
    tag = tag | vcfbwt::pfp::ParserVCF::UNCOMPRESSED;

    worker.init(params, out_prefix, reference_parse, tag);
    main_parser.register_worker(worker);

    // Run
    worker(vcf[0]);

    // Close the main parser
    main_parser.close();

    // Generate the desired outcome from the test files, reference first
    std::vector<vcfbwt::char_type> what_it_should_be;
    what_it_should_be.insert(what_it_should_be.end(), 1, vcfbwt::pfp::DOLLAR);
    what_it_should_be.insert(what_it_should_be.end(), vcf.get_reference().begin(), vcf.get_reference().end());
    what_it_should_be.insert(what_it_should_be.end(), params.w - 1, vcfbwt::pfp::DOLLAR_PRIME);
    what_it_should_be.insert(what_it_should_be.end(), 1, vcfbwt::pfp::DOLLAR_SEQUENCE);

    std::string test_sample_path = testfiles_dir + "/HG00096_chrY_H1.fa.gz";
    std::ifstream in_stream(test_sample_path);
    zstr::istream is(in_stream);
    std::string line, from_fasta;
    while (getline(is, line)) { if ( not (line.empty() or line[0] == '>') ) { from_fasta.append(line); } }

    what_it_should_be.insert(what_it_should_be.end(), from_fasta.begin(), from_fasta.end());
    what_it_should_be.insert(what_it_should_be.end(), params.w - 1, vcfbwt::pfp::DOLLAR_PRIME);
    //what_it_should_be.append(1, vcfbwt::pfp::DOLLAR_SEQUENCE);
    what_it_should_be.insert(what_it_should_be.end(), params.w, vcfbwt::pfp::DOLLAR);

    // Check
    bool check = unparse_and_check<vcfbwt::char_type>(out_prefix, what_it_should_be, params.w, vcfbwt::pfp::DOLLAR);
    REQUIRE(check);
}

TEST_CASE( "Reference + Sample HG00096, WITH acceleration", "[PFP algorithm]" )
{
    std::string vcf_file_name = testfiles_dir + "/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz";
    std::string ref_file_name = testfiles_dir + "/Y.fa.gz";
    vcfbwt::VCF vcf(ref_file_name, vcf_file_name, "", 1);

    // Produce dictionary and parsing
    vcfbwt::pfp::Params params;
    params.w = w_global; params.p = p_global;
    params.use_acceleration = true;
    params.compute_occurrences = true;
    vcfbwt::pfp::ReferenceParse reference_parse(vcf.get_reference(), params);

    std::string out_prefix = testfiles_dir + "/parser_out";
    vcfbwt::pfp::ParserVCF main_parser(params, out_prefix, reference_parse);

    vcfbwt::pfp::ParserVCF worker;
    std::size_t tag = 0;
    tag = tag | vcfbwt::pfp::ParserVCF::WORKER;
    tag = tag | vcfbwt::pfp::ParserVCF::UNCOMPRESSED;

    worker.init(params, out_prefix, reference_parse, tag);
    main_parser.register_worker(worker);

    // Run
    worker(vcf[0]);

    // Close the main parser
    main_parser.close();

    // Generate the desired outcome from the test files, reference first
    std::vector<vcfbwt::char_type> what_it_should_be;
    what_it_should_be.insert(what_it_should_be.end(),1, vcfbwt::pfp::DOLLAR);
    what_it_should_be.insert(what_it_should_be.end(), vcf.get_reference().begin(), vcf.get_reference().end());
    what_it_should_be.insert(what_it_should_be.end(), params.w - 1, vcfbwt::pfp::DOLLAR_PRIME);
    what_it_should_be.insert(what_it_should_be.end(), 1, vcfbwt::pfp::DOLLAR_SEQUENCE);

    std::string test_sample_path = testfiles_dir + "/HG00096_chrY_H1.fa.gz";
    std::ifstream in_stream(test_sample_path);
    zstr::istream is(in_stream);
    std::string line, from_fasta;
    while (getline(is, line)) { if ( not (line.empty() or line[0] == '>') ) { from_fasta.append(line); } }

    what_it_should_be.insert(what_it_should_be.end(), from_fasta.begin(), from_fasta.end());
    what_it_should_be.insert(what_it_should_be.end(), params.w - 1, vcfbwt::pfp::DOLLAR_PRIME);
    //what_it_should_be.append(1, vcfbwt::pfp::DOLLAR_SEQUENCE);
    what_it_should_be.insert(what_it_should_be.end(), params.w, vcfbwt::pfp::DOLLAR);

    // Check
    bool check = unparse_and_check<vcfbwt::char_type>(out_prefix, what_it_should_be, params.w, vcfbwt::pfp::DOLLAR);
    REQUIRE(check);
}

TEST_CASE( "Sample: HG00096, twice chromosome Y", "[VCF parser]" )
{
    std::vector<std::string> vcf_file_names =
            {
                    testfiles_dir + "/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz",
                    testfiles_dir + "/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz"
            };

    std::vector<std::string> ref_file_names =
            {
                    testfiles_dir + "/Y.fa.gz",
                    testfiles_dir + "/Y.fa.gz"
            };

    vcfbwt::VCF vcf(ref_file_names, vcf_file_names, "", 1);

    REQUIRE(vcf[0].id() == "HG00096");

    std::string test_sample_path = testfiles_dir + "/HG00096_chrY_H1.fa.gz";
    std::ifstream in_stream(test_sample_path);

    REQUIRE(vcfbwt::is_gzipped(in_stream));

    zstr::istream is(in_stream);
    std::string line, from_fasta;
    while (getline(is, line)) { if ( not (line.empty() or line[0] == '>') ) { from_fasta.append(line); } }
    from_fasta.push_back(vcfbwt::pfp::DOLLAR_PRIME);
    from_fasta.append(from_fasta);
    from_fasta.pop_back();


    vcfbwt::Sample::iterator it(vcf[0]);
    std::string from_vcf;
    while (not it.end()) { from_vcf.push_back(*it); ++it; }

    std::size_t i = 0;
    while ( ((i < from_vcf.size()) and (i < from_fasta.size()))
            and (from_vcf[i] == from_fasta[i])) { i++; }
    REQUIRE(((i == (from_vcf.size())) and (i == (from_fasta.size()))));
}

TEST_CASE( "Sample: HG00096, fasta", "[PFP Algo]" )
{
    // Produce dictionary and parsing
    vcfbwt::pfp::Params params;
    params.w = w_global; params.p = p_global;
    params.compute_occurrences = true;

    std::string test_sample_path = testfiles_dir + "/HG00096_chrY_H1.fa.gz";
    std::string out_prefix = testfiles_dir + "/HG00096_chrY_H1_tpfa";
    vcfbwt::pfp::ParserFasta main_parser(params, test_sample_path, out_prefix);

    // Run
    main_parser();

    // Close the main parser
    main_parser.close();

    // Generate the desired outcome from the test files, reference first
    std::vector<vcfbwt::char_type> what_it_should_be;
    what_it_should_be.insert(what_it_should_be.end(), 1, vcfbwt::pfp::DOLLAR);

    std::ifstream in_stream(test_sample_path);
    zstr::istream is(in_stream);
    std::string line, from_fasta;
    while (getline(is, line)) { if ( not (line.empty() or line[0] == '>') ) { from_fasta.append(line); } }

    what_it_should_be.insert(what_it_should_be.end(), from_fasta.begin(), from_fasta.end());
    what_it_should_be.insert(what_it_should_be.end(), params.w - 1, vcfbwt::pfp::DOLLAR_PRIME);
    what_it_should_be.insert(what_it_should_be.end(), 1, vcfbwt::pfp::DOLLAR_SEQUENCE);
    what_it_should_be.insert(what_it_should_be.end(), params.w, vcfbwt::pfp::DOLLAR);

    // Check
    bool check = unparse_and_check<vcfbwt::char_type>(out_prefix, what_it_should_be, params.w, vcfbwt::pfp::DOLLAR);
    REQUIRE(check);
}

TEST_CASE( "Sample: HG00096, text", "[PFP Algo]" )
{
    // Produce dictionary and parsing
    vcfbwt::pfp::Params params;
    params.w = w_global; params.p = p_global;
    params.compute_occurrences = true;

    std::string test_sample_path = testfiles_dir + "/HG00096_chrY_H1.fa.gz";
    std::string out_prefix = testfiles_dir + "/HG00096_chrY_H1_tptxt";
    vcfbwt::pfp::ParserText main_parser(params, test_sample_path, out_prefix);

    // Run
    main_parser();

    // Close the main parser
    main_parser.close();

    // Generate the desired outcome from the test files, reference first
    std::vector<vcfbwt::char_type> what_it_should_be;
    what_it_should_be.insert(what_it_should_be.end(), 1, vcfbwt::pfp::DOLLAR);

    std::ifstream in_stream(test_sample_path);
    zstr::istream is(in_stream);
    std::string from_text_file((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());

    what_it_should_be.insert(what_it_should_be.end(), from_text_file.begin(), from_text_file.end());
    what_it_should_be.insert(what_it_should_be.end(), params.w, vcfbwt::pfp::DOLLAR);

    // Check
    bool check = unparse_and_check<vcfbwt::char_type>(out_prefix, what_it_should_be, params.w, vcfbwt::pfp::DOLLAR);
    REQUIRE(check);
}

TEST_CASE( "Sample: HG00096, integers", "[PFP Algo]" )
{
    // Produce dictionary and parsing
    vcfbwt::pfp::Params params;
    params.w = w_global; params.p = p_global;
    params.compute_occurrences = true;
    params.integers_shift = 0;

    // Test file
    std::ofstream test_file(testfiles_dir + "/repetitive_int32_t.bin");
    std::vector<int32_t> what_it_should_be;
    what_it_should_be.insert(what_it_should_be.end(), 1, vcfbwt::pfp::DOLLAR);
    for (std::size_t i = 0; i < 10; i++)
    {
        for (int32_t j = 10; j < 1000; j++)
        {
            test_file.write((char*) &j, sizeof(int32_t));
            what_it_should_be.emplace_back(j);
        }
    }
    test_file.close();
    what_it_should_be.insert(what_it_should_be.end(), params.w, vcfbwt::pfp::DOLLAR);

    std::string test_sample_path = testfiles_dir + "/repetitive_int32_t.bin";
    std::string out_prefix = testfiles_dir + "/repetitive_int32_t_tpintegers";
    vcfbwt::pfp::ParserIntegers main_parser(params, test_sample_path, out_prefix);

    // Run
    main_parser();

    // Close the main parser
    main_parser.close();

    // Check
    bool check = unparse_and_check<int32_t>(out_prefix, what_it_should_be, params.w, vcfbwt::pfp::DOLLAR);
    REQUIRE(check);
}

//------------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
    Catch::Session session;

    using namespace Catch::Clara;

    auto cli = session.cli() |
    Opt( testfiles_dir, "dir" ) ["--test-dir"] ("specify the directory containing the test dna sequences files") |
    Opt( w_global, "int" ) ["-W"] ("specify w") |
    Opt( p_global, "int" ) ["-P"] ("specify p");

    session.cli(cli);

    int returnCode = session.applyCommandLine(argc, argv);

    if( returnCode != 0 ) return returnCode;

    spdlog::info("Tests running with w: {}\tp: {}", w_global, p_global);

    session.run();
}

