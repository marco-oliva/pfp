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
    std::size_t threshold;

    app.add_option("-o,--out-file", out_file, "Output file")->check(CLI::NonexistentPath)->required();
    app.add_option("-d,--dictionary", dict_file, "Dictionary file")->check(CLI::ExistingFile)->required();
    app.add_option("-p,--parse", parse_file, "Parse file")->check(CLI::ExistingFile)->required();
    app.add_option("-t,--threshold", threshold, "Threshold")->required();
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


    spdlog::info("Reading Parse");
    std::vector<vcfbwt::size_type> parse;
    vcfbwt::pfp::Parser::read_parse(parse_file, parse);

    std::size_t dict_size = 0;
    for (auto& d : dict) { dict_size += d.size(); }
    spdlog::info("Parse + Dicionary: {} bytes", (parse.size() * sizeof(vcfbwt::size_type)) + (dict_size));

    vcfbwt::pfp::AuPair au_pair_algo(dict, parse, window_size, out_file);

    std::size_t iteration = 0, max_iterations = 10, removed = 1;
    while (removed != 0 and iteration < max_iterations)
    {
        removed = au_pair_algo.compress(threshold);
        spdlog::info("[{}]Removed: {} bytes", iteration, removed);
        iteration++;
    }
}