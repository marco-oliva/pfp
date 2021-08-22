//
//  unparse.cpp.c
//
//  Copyright 2021 Marco Oliva. All rights reserved.
//

#include <CLI/CLI.hpp>
#include <version.hpp>
#include <utils.hpp>
#include <pfp_algo.hpp>

int main(int argc, char **argv)
{
    CLI::App app("unparse");
    
    std::string out_file;
    std::string dict_file;
    std::string parse_file;
    std::size_t window_size;
    bool size_only = false;
    
    app.add_option("-o,--out-file", out_file, "Output file")->check(CLI::NonexistentPath);
    app.add_option("-d,--dictionary", dict_file, "Dictionary file")->check(CLI::ExistingFile)->required();
    app.add_option("-p,--parse", parse_file, "Parse file")->check(CLI::ExistingFile)->required();
    app.add_option("-w, --window", window_size, "Window size")->required();
    app.add_flag("--size-only", size_only, "Only compute size of unparsed text");
    app.add_flag_callback("--version",vcfbwt::Version::print,"Version");
    app.allow_windows_style_options();
    
    CLI11_PARSE(app, argc, argv);
    
    // Print out configurations
    spdlog::info("Current Configuration:\n{}", app.config_to_str(true,true));
    
    // Unparse
    spdlog::info("Reading dictionary");
    std::vector<std::string> dict;
    vcfbwt::pfp::ParserUtils::read_dictionary(dict_file, dict);
    
    std::ifstream parse_stream(parse_file, std::ios::binary);
    
    spdlog::info("Start unparsing");
    
    std::size_t unparsed_text_size = 0;
    std::ofstream unparsed;
    if (not size_only) { unparsed.open(out_file); }
    while (not parse_stream.eof())
    {
        vcfbwt::size_type p = 0;
        parse_stream.read((char*)&p, sizeof(vcfbwt::size_type));
        
        if (p > dict.size()) { spdlog::error("Something wrong in the parse"); exit(EXIT_FAILURE); }
        if (p == 0) { spdlog::error("Element in the parse is 0, skipping"); continue; }
        if (dict[p - 1].size() <= window_size) { spdlog::error("Phrase shorter than {} in the parse",window_size); exit(EXIT_FAILURE);  }
        std::string_view dict_string(&(dict[p - 1][0]), dict[p - 1].size() - window_size);
        
        unparsed_text_size += dict_string.size();
        if (not size_only) { unparsed << dict_string; }
       
    }
    if (not size_only) { unparsed << vcfbwt::pfp::DOLLAR; }
    
    spdlog::info("Unparsed text size: {}", unparsed_text_size);
}

