//
//  check_integrity.cpp
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
    CLI::App app("AuPair");

    std::string input_dict_path;
    std::string input_parse_path;
    std::size_t window_size;

    app.add_option("-d,--dictionary", input_dict_path, "Dictionary")->required();
    app.add_option("-p,--parse", input_parse_path, "Parse")->required();
    app.add_option("-w, --window", window_size, "Window size")->required();
    app.add_flag_callback("--version",vcfbwt::Version::print,"Version");
    app.allow_windows_style_options();

    CLI11_PARSE(app, argc, argv);

    // Print out configurations
    spdlog::info("Current Configuration:\n{}", app.config_to_str(true,true));

    // Unparse
    spdlog::info("Reading dictionary");
    std::vector<std::string> dict;
    vcfbwt::pfp::ParserUtils::read_dictionary(input_dict_path, dict);

    std::ifstream parse_stream(input_parse_path, std::ios::binary);
    vcfbwt::size_type curr_parse_element = 0, prev_parse_element = 0;

    // Read first element
    parse_stream.read((char*)&prev_parse_element, sizeof(vcfbwt::size_type));
    if (prev_parse_element > dict.size()) { spdlog::error("Something wrong in the parse"); std::exit(EXIT_FAILURE); }
    if (prev_parse_element != 1) { spdlog::error("First element of the parse sould be 1"); std::exit(EXIT_FAILURE); }

    spdlog::info("Checking parse");
    while (not parse_stream.eof())
    {
        parse_stream.read((char*)&curr_parse_element, sizeof(vcfbwt::size_type));

        if (curr_parse_element > dict.size()) { spdlog::error("Something wrong in the parse"); exit(EXIT_FAILURE); }
        if (curr_parse_element == 0) { spdlog::error("Element in the parse is 0, skipping"); continue; }

        // Check Trigger Strings
        std::string_view ts_prev(&(dict[prev_parse_element - 1][dict[prev_parse_element - 1].size() - window_size]), window_size);
        std::string_view ts_curr(&(dict[curr_parse_element - 1][0]), window_size);

        if (ts_prev != ts_curr)
        {
            spdlog::error("\n[{}, {}] {}\n[{}, {}] {}",
                          prev_parse_element, vcfbwt::string_hash(dict[prev_parse_element - 1].c_str(), dict[prev_parse_element - 1].size()),
                          dict[prev_parse_element - 1].substr(dict[prev_parse_element - 1].size() - (2 * window_size), (2 * window_size)),
                          curr_parse_element, vcfbwt::string_hash(dict[curr_parse_element - 1].c_str(), dict[curr_parse_element - 1].size()),
                          dict[curr_parse_element - 1].substr(0, (2 * window_size)));
        }

        prev_parse_element = curr_parse_element;
    }

    parse_stream.close();
}