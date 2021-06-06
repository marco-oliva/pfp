//
//  aupair.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#include <CLI/CLI.hpp>
#include <version.hpp>
#include <utils.hpp>
#include <pfp_algo.hpp>

vcfbwt::size_type deleted_element = 0; // parse elements start at 1

int main(int argc, char **argv)
{
    CLI::App app("Aupair");

    std::string out_file;
    std::string dict_file;
    std::string parse_file;
    std::string occ_file;
    std::size_t window_size;

    app.add_option("-o,--out-file", out_file, "Output file")->check(CLI::NonexistentPath)->required();
    app.add_option("-d,--dictionary", dict_file, "Dictionary file")->check(CLI::ExistingFile)->required();
    app.add_option("-p,--parse", parse_file, "Parse file")->check(CLI::ExistingFile)->required();
    app.add_option("-f,--occurrences", occ_file, "Occurrences file")->check(CLI::ExistingFile)->required();
    app.add_option("-w, --window", window_size, "Window size")->required();
    app.add_flag_callback("--version",vcfbwt::Version::print,"Version");
    app.allow_windows_style_options();

    CLI11_PARSE(app, argc, argv);

    // Print out configurations
    spdlog::info("Current Configuration:\n{}", app.config_to_str(true,true));

    // Read Input
    spdlog::info("Reading Dictionary");
    std::vector<std::string> dict;
    vcfbwt::pfp::Parser::read_dictionary(dict_file, dict);

//    spdlog::info("Reading Occurrences");
//    std::vector<vcfbwt::size_type> occurrences(dict.size());
//    std::ifstream occurrences_stream(occ_file, std::ios::binary);
//    occurrences_stream.read((char*) &occurrences[0], occurrences.size() * sizeof(vcfbwt::size_type));

    spdlog::info("Reading Parse");
    std::vector<vcfbwt::size_type> parse;
    vcfbwt::pfp::Parser::read_parse(parse_file, parse);

    spdlog::info("Initializing Trigger strings map");
    std::unordered_map<std::string, std::vector<std::size_t>> trigger_strings_info;

}