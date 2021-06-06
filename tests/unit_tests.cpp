//
//  unit_tests.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>

#include <vcf.hpp>
#include <utils.hpp>
#include <pfp_algo.hpp>

//------------------------------------------------------------------------------

struct listener : Catch::TestEventListenerBase
{
    using TestEventListenerBase::TestEventListenerBase;
    
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

////------------------------------------------------------------------------------
//TEST_CASE( "Initializaiton", "[KR Window]" )
//{
//    std::string test_string = "12345";
//
//    vcfbwt::KarpRabinHash kr_window(5);
//    kr_window.set_constant(256); kr_window.set_prime(2147483647);
//    kr_window.initialize(test_string);
//
//    REQUIRE(kr_window.get_hash() == 842216599);
//}
//
//TEST_CASE( "Update 1 charachter", "[KR Window]" )
//{
//    std::string test_string = "12345";
//
//    vcfbwt::KarpRabinHash kr_window(5);
//    kr_window.set_constant(256); kr_window.set_prime(2147483647);
//    kr_window.initialize(test_string);
//
//    kr_window.update('1', '6');
//
//    REQUIRE(kr_window.get_hash() == 859059610);
//}
//
//TEST_CASE( "Update 2 charachter", "[KR Window]" )
//{
//    std::string test_string = "12345";
//
//    vcfbwt::KarpRabinHash kr_window(5);
//    kr_window.set_constant(256); kr_window.set_prime(2147483647);
//    kr_window.initialize(test_string);
//
//    kr_window.update('1', '6');
//    kr_window.update('2', '7');
//
//    REQUIRE(kr_window.get_hash() == 875902621);
//}
//
//TEST_CASE( "Periodic string", "[KR Window]" )
//{
//    std::string test_string = "11111";
//
//    vcfbwt::KarpRabinHash kr_window(5);
//    kr_window.set_constant(256); kr_window.set_prime(2147483647);
//    kr_window.initialize(test_string);
//
//    vcfbwt::hash_type before = kr_window.get_hash();
//    kr_window.update('1', '1');
//    vcfbwt::hash_type after = kr_window.get_hash();
//
//    REQUIRE(before == 825307539);
//    REQUIRE(after  == 825307539);
//}
//
//TEST_CASE( "Periodic string of Ns", "[KR Window]" )
//{
//    std::string test_string = "NNNNNNNNNNNNNNNNNNNN";
//
//    vcfbwt::KarpRabinHash kr_window(test_string.size());
//    kr_window.set_constant(256); kr_window.set_prime(2147483647);
//    kr_window.initialize(test_string);
//
//    vcfbwt::hash_type before = kr_window.get_hash();
//    kr_window.update('N', 'N');
//    vcfbwt::hash_type after = kr_window.get_hash();
//
//    REQUIRE(before == 2071690116);
//    REQUIRE(after  == 2071690116);
//}
//
////------------------------------------------------------------------------------
//
//TEST_CASE( "Constructor", "[VCF parser]" )
//{
//    std::string vcf_file_name = testfiles_dir + "/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz";
//    std::string ref_file_name = testfiles_dir + "/Y.fa.gz";
//    vcfbwt::VCF vcf(ref_file_name, vcf_file_name);
//
//    // read_samples list from file
//    std::ifstream samples_file(testfiles_dir + "/samples_list.txt");
//    std::vector<std::string> samples_list;
//    std::copy(std::istream_iterator<std::string>(samples_file),
//              std::istream_iterator<std::string>(),
//              std::back_inserter(samples_list));
//
//    bool all_match = true;
//    for (std::size_t i = 0; i < vcf.size(); i++)
//    {
//        all_match = all_match  & (vcf[i].id() == samples_list[i]);
//    }
//
//    REQUIRE(vcf.size() == samples_list.size());
//    REQUIRE(all_match);
//}
//
//TEST_CASE( "Sample: HG00096", "[VCF parser]" )
//{
//    std::string vcf_file_name = testfiles_dir + "/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz";
//    std::string ref_file_name = testfiles_dir + "/Y.fa.gz";
//    vcfbwt::VCF vcf(ref_file_name, vcf_file_name, 1);
//
//    REQUIRE(vcf[0].id() == "HG00096");
//
//    std::string test_sample_path = testfiles_dir + "/HG00096_chrY_H1.fa.gz";
//    std::ifstream in_stream(test_sample_path);
//
//    REQUIRE(vcfbwt::is_gzipped(in_stream));
//
//    zstr::istream is(in_stream);
//    std::string line, from_fasta;
//    while (getline(is, line)) { if ( not (line.empty() or line[0] == '>') ) { from_fasta.append(line); } }
//
//    vcfbwt::Sample::iterator it(vcf[0]);
//    std::string from_vcf;
//    while (not it.end()) { from_vcf.push_back(*it); ++it; }
//
//    std::size_t i = 0;
//    while ( ((i < from_vcf.size()) and (i < from_fasta.size()))
//            and (from_vcf[i] == from_fasta[i])) { i++; }
//    REQUIRE(((i == (from_vcf.size())) and (i == (from_fasta.size()))));
//}
//
//TEST_CASE( "Sample: HG00103", "[VCF parser]" )
//{
//    std::string vcf_file_name = testfiles_dir + "/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz";
//    std::string ref_file_name = testfiles_dir + "/Y.fa.gz";
//    vcfbwt::VCF vcf(ref_file_name, vcf_file_name, 3);
//
//    REQUIRE(vcf[2].id() == "HG00103");
//
//    std::string test_sample_path = testfiles_dir + "/HG00103_chrY_H1.fa.gz";
//    std::ifstream in_stream(test_sample_path);
//
//    REQUIRE(vcfbwt::is_gzipped(in_stream));
//
//    zstr::istream is(in_stream);
//    std::string line, from_fasta;
//    while (getline(is, line)) { if ( not (line.empty() or line[0] == '>') ) { from_fasta.append(line); } }
//
//    vcfbwt::Sample::iterator it(vcf[2]);
//    std::string from_vcf;
//    while (not it.end()) { from_vcf.push_back(*it); ++it; }
//
//    std::size_t i = 0;
//    while ( ((i < from_vcf.size()) and (i < from_fasta.size()))
//    and (from_vcf[i] == from_fasta[i])) { i++; }
//    REQUIRE(((i == (from_vcf.size())) and (i == (from_fasta.size()))));
//}
//
//TEST_CASE( "Reference + Sample HG00096, No acceleration", "[PFP algorithm]" )
//{
//    std::string vcf_file_name = testfiles_dir + "/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz";
//    std::string ref_file_name = testfiles_dir + "/Y.fa.gz";
//    vcfbwt::VCF vcf(ref_file_name, vcf_file_name, 1);
//
//    // Only work on sample HG00096
//    vcf.set_max_samples(1);
//
//    // Produce dictionary and parsing
//    vcfbwt::pfp::Params params;
//    params.w = w_global; params.p = p_global;
//    params.use_acceleration = false;
//    vcfbwt::pfp::ReferenceParse reference_parse(vcf.get_reference(), params);
//
//    std::string out_prefix = testfiles_dir + "/parser_out";
//    vcfbwt::pfp::Parser main_parser(params, out_prefix, reference_parse);
//
//    vcfbwt::pfp::Parser worker;
//    std::size_t tag = 0;
//    tag = tag | vcfbwt::pfp::Parser::WORKER;
//    tag = tag | vcfbwt::pfp::Parser::UNCOMPRESSED;
//    tag = tag | vcfbwt::pfp::Parser::LAST;
//
//    worker.init(params, out_prefix, reference_parse, tag);
//    main_parser.register_worker(worker);
//
//    // Run
//    worker(vcf[0]);
//
//    // Close the main parser
//    main_parser.close();
//
//    // Generate the desired outcome from the test files, reference first
//    std::string what_it_should_be;
//    what_it_should_be.append(1, vcfbwt::pfp::DOLLAR);
//    what_it_should_be.insert(what_it_should_be.end(), vcf.get_reference().begin(), vcf.get_reference().end());
//    what_it_should_be.append(params.w, vcfbwt::pfp::DOLLAR_PRIME);
//
//    std::string test_sample_path = testfiles_dir + "/HG00096_chrY_H1.fa.gz";
//    std::ifstream in_stream(test_sample_path);
//    zstr::istream is(in_stream);
//    std::string line, from_fasta;
//    while (getline(is, line)) { if ( not (line.empty() or line[0] == '>') ) { from_fasta.append(line); } }
//
//    what_it_should_be.insert(what_it_should_be.end(), from_fasta.begin(), from_fasta.end());
//    what_it_should_be.append(params.w, vcfbwt::pfp::DOLLAR);
//
//    // Unparse
//    std::vector<vcfbwt::size_type>  parse;
//    vcfbwt::pfp::Parser::read_parse(out_prefix + ".parse", parse);
//    std::vector<std::string> dict;
//    vcfbwt::pfp::Parser::read_dictionary(out_prefix + ".dict", dict);
//
//    std::string unparsed;
//    for (auto& p : parse)
//    {
//        if (p > dict.size()) { spdlog::error("Something wrong in the parse"); exit(EXIT_FAILURE); }
//        std::string dict_string = dict[p - 1].substr(0, dict[p - 1].size() - params.w);
//        unparsed.insert(unparsed.end(), dict_string.begin(), dict_string.end());
//    }
//    unparsed.append(params.w, vcfbwt::pfp::DOLLAR);
//
//
//    // Compare the two strings
//    std::size_t i = 0;
//    while ( ((i < unparsed.size()) and (i < what_it_should_be.size()))
//    and (unparsed[i] == what_it_should_be[i])) { i++; }
//    REQUIRE(((i == (unparsed.size())) and (i == (what_it_should_be.size()))));
//}
//
//TEST_CASE( "Reference + Sample HG00096, WITH acceleration", "[PFP algorithm]" )
//{
//    std::string vcf_file_name = testfiles_dir + "/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz";
//    std::string ref_file_name = testfiles_dir + "/Y.fa.gz";
//    vcfbwt::VCF vcf(ref_file_name, vcf_file_name, 1);
//
//    // Only work on sample HG00096
//    vcf.set_max_samples(1);
//
//    // Produce dictionary and parsing
//    vcfbwt::pfp::Params params;
//    params.w = w_global; params.p = p_global;
//    params.use_acceleration = true;
//    vcfbwt::pfp::ReferenceParse reference_parse(vcf.get_reference(), params);
//
//    std::string out_prefix = testfiles_dir + "/parser_out";
//    vcfbwt::pfp::Parser main_parser(params, out_prefix, reference_parse);
//
//    vcfbwt::pfp::Parser worker;
//    std::size_t tag = 0;
//    tag = tag | vcfbwt::pfp::Parser::WORKER;
//    tag = tag | vcfbwt::pfp::Parser::UNCOMPRESSED;
//    tag = tag | vcfbwt::pfp::Parser::LAST;
//
//    worker.init(params, out_prefix, reference_parse, tag);
//    main_parser.register_worker(worker);
//
//    // Run
//    worker(vcf[0]);
//
//    // Close the main parser
//    main_parser.close();
//
//    // Generate the desired outcome from the test files, reference first
//    std::string what_it_should_be;
//    what_it_should_be.append(1, vcfbwt::pfp::DOLLAR);
//    what_it_should_be.insert(what_it_should_be.end(), vcf.get_reference().begin(), vcf.get_reference().end());
//    what_it_should_be.append(params.w, vcfbwt::pfp::DOLLAR_PRIME);
//
//    std::string test_sample_path = testfiles_dir + "/HG00096_chrY_H1.fa.gz";
//    std::ifstream in_stream(test_sample_path);
//    zstr::istream is(in_stream);
//    std::string line, from_fasta;
//    while (getline(is, line)) { if ( not (line.empty() or line[0] == '>') ) { from_fasta.append(line); } }
//
//    what_it_should_be.insert(what_it_should_be.end(), from_fasta.begin(), from_fasta.end());
//    what_it_should_be.append(params.w, vcfbwt::pfp::DOLLAR);
//
//    // Unparse
//    std::vector<vcfbwt::size_type>  parse;
//    vcfbwt::pfp::Parser::read_parse(out_prefix + ".parse", parse);
//    std::vector<std::string> dict;
//    vcfbwt::pfp::Parser::read_dictionary(out_prefix + ".dict", dict);
//
//    std::string unparsed;
//    for (auto& p : parse)
//    {
//        if (p > dict.size()) { spdlog::error("Something wrong in the parse"); exit(EXIT_FAILURE); }
//        std::string dict_string = dict[p - 1].substr(0, dict[p - 1].size() - params.w);
//        unparsed.insert(unparsed.end(), dict_string.begin(), dict_string.end());
//    }
//    unparsed.append(params.w, vcfbwt::pfp::DOLLAR);
//
//    // Compare the two strings
//    std::size_t i = 0;
//    while ( ((i < unparsed.size()) and (i < what_it_should_be.size()))
//            and (unparsed[i] == what_it_should_be[i])) { i++; }
//    spdlog::info("First missmatch: {}", i);
//    REQUIRE(((i == (unparsed.size())) and (i == (what_it_should_be.size()))));
//}
//
//TEST_CASE( "Sample: HG00096, twice chromosome Y", "[VCF parser]" )
//{
//    std::vector<std::string> vcf_file_names =
//            {
//                    testfiles_dir + "/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz",
//                    testfiles_dir + "/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz"
//            };
//
//    std::vector<std::string> ref_file_names =
//            {
//                    testfiles_dir + "/Y.fa.gz",
//                    testfiles_dir + "/Y.fa.gz"
//            };
//
//    vcfbwt::VCF vcf(ref_file_names, vcf_file_names, 1);
//
//    REQUIRE(vcf[0].id() == "HG00096");
//
//    std::string test_sample_path = testfiles_dir + "/HG00096_chrY_H1.fa.gz";
//    std::ifstream in_stream(test_sample_path);
//
//    REQUIRE(vcfbwt::is_gzipped(in_stream));
//
//    zstr::istream is(in_stream);
//    std::string line, from_fasta;
//    while (getline(is, line)) { if ( not (line.empty() or line[0] == '>') ) { from_fasta.append(line); } }
//    from_fasta.push_back(vcfbwt::pfp::SPECIAL_TYPES::DOLLAR_PRIME);
//    from_fasta.append(from_fasta);
//    from_fasta.pop_back();
//
//
//    vcfbwt::Sample::iterator it(vcf[0]);
//    std::string from_vcf;
//    while (not it.end()) { from_vcf.push_back(*it); ++it; }
//
//    std::size_t i = 0;
//    while ( ((i < from_vcf.size()) and (i < from_fasta.size()))
//            and (from_vcf[i] == from_fasta[i])) { i++; }
//    REQUIRE(((i == (from_vcf.size())) and (i == (from_fasta.size()))));
//}

TEST_CASE( "AuPair Reference + Sample HG00096, No acceleration", "[AuPair]" )
{
    std::string vcf_file_name = testfiles_dir + "/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz";
    std::string ref_file_name = testfiles_dir + "/Y.fa.gz";
    vcfbwt::VCF vcf(ref_file_name, vcf_file_name, 1);

    // Only work on sample HG00096
    vcf.set_max_samples(1);

    // Produce dictionary and parsing
    vcfbwt::pfp::Params params;
    params.w = w_global; params.p = p_global;
    params.use_acceleration = false;
    vcfbwt::pfp::ReferenceParse reference_parse(vcf.get_reference(), params);

    std::string out_prefix = testfiles_dir + "/parser_out";
    vcfbwt::pfp::Parser main_parser(params, out_prefix, reference_parse);

    vcfbwt::pfp::Parser worker;
    std::size_t tag = 0;
    tag = tag | vcfbwt::pfp::Parser::WORKER;
    tag = tag | vcfbwt::pfp::Parser::UNCOMPRESSED;
    tag = tag | vcfbwt::pfp::Parser::LAST;

    worker.init(params, out_prefix, reference_parse, tag);
    main_parser.register_worker(worker);

    // Run
    worker(vcf[0]);

    // Close the main parser
    main_parser.close();
    
    // Generate the desired outcome from the test files, reference first
    std::string what_it_should_be;
    what_it_should_be.append(1, vcfbwt::pfp::DOLLAR);
    what_it_should_be.insert(what_it_should_be.end(), vcf.get_reference().begin(), vcf.get_reference().end());
    what_it_should_be.append(params.w, vcfbwt::pfp::DOLLAR_PRIME);

    std::string test_sample_path = testfiles_dir + "/HG00096_chrY_H1.fa.gz";
    std::ifstream in_stream(test_sample_path);
    zstr::istream is(in_stream);
    std::string line, from_fasta;
    while (getline(is, line)) { if ( not (line.empty() or line[0] == '>') ) { from_fasta.append(line); } }

    what_it_should_be.insert(what_it_should_be.end(), from_fasta.begin(), from_fasta.end());
    what_it_should_be.append(params.w, vcfbwt::pfp::DOLLAR);
    
    // Unparse for AuPair
    std::vector<vcfbwt::size_type>  parse;
    vcfbwt::pfp::Parser::read_parse(out_prefix + ".parse", parse);
    std::vector<std::string> dict;
    vcfbwt::pfp::Parser::read_dictionary(out_prefix + ".dict", dict);
    
    std::string compressed_out_prefix = testfiles_dir + "/compressed";
    vcfbwt::pfp::AuPair au_pair_algo(dict, parse, w_global, compressed_out_prefix);
    
    std::size_t removed = au_pair_algo.compress(100);
    spdlog::info("Removed: {} bytes", removed);
    
    au_pair_algo.close();
    
    // Unparse to check
    std::vector<vcfbwt::size_type>  parse_compressed;
    vcfbwt::pfp::Parser::read_parse(compressed_out_prefix + ".parse", parse_compressed);
    std::vector<std::string> dict_compressed;
    vcfbwt::pfp::Parser::read_dictionary(compressed_out_prefix + ".dict", dict_compressed);
    
    std::string unparsed;
    for (auto& p : parse_compressed)
    {
        if (p > dict_compressed.size()) { spdlog::error("Something wrong in the parse"); exit(EXIT_FAILURE); }
        std::string dict_string = dict_compressed[p - 1].substr(0, dict_compressed[p - 1].size() - params.w);
        unparsed.insert(unparsed.end(), dict_string.begin(), dict_string.end());
    }
    unparsed.append(params.w, vcfbwt::pfp::DOLLAR);


    // Compare the two strings
    std::size_t i = 0;
    while ( ((i < unparsed.size()) and (i < what_it_should_be.size()))
    and (unparsed[i] == what_it_should_be[i])) { i++; }
    REQUIRE(((i == (unparsed.size())) and (i == (what_it_should_be.size()))));
}

//------------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
    Catch::Session session;

    using namespace Catch::clara;

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