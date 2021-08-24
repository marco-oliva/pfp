//
//  aupair.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#include <CLI/CLI.hpp>
#include <version.hpp>
#include <utils.hpp>
#include <pfp_algo.hpp>
#include <au_pair_algo.hpp>

vcfbwt::size_type deleted_element = 0; // parse elements start at 1

int main(int argc, char **argv)
{
    CLI::App app("AuPair");

    std::string out_file;
    std::string input_prefix;
    std::size_t window_size;
    std::size_t threshold = 0;

    app.add_option("-o,--out-file", out_file, "Output file")->required();
    app.add_option("-i,--input", input_prefix, "Input Prefix")->required();
    app.add_option("-t,--threshold", threshold, "Threshold");
    app.add_option("-w, --window", window_size, "Window size")->required();
    app.add_flag_callback("--version",vcfbwt::Version::print,"Version");
    app.allow_windows_style_options();

    CLI11_PARSE(app, argc, argv);

    // Print out configurations
    spdlog::info("Current Configuration:\n{}", app.config_to_str(true,true));

    vcfbwt::pfp::AuPair au_pair_algo(input_prefix, window_size);

    std::set<std::string_view> removed_trigger_strings;
    //std::size_t removed_bytes = au_pair_algo.remove_simple(removed_trigger_strings);
    std::size_t removed_bytes = au_pair_algo.compress(removed_trigger_strings, threshold);

    spdlog::info("Removed {} bytes, can be inaccurate if --batch-size > 1", removed_bytes);
    spdlog::info("Removed {} trigger strings", removed_trigger_strings.size());
    
    std::ofstream out_ts(out_file);
    for (auto& rts : removed_trigger_strings)
    {
        for (auto& c : rts) { out_ts.put(c); }
        out_ts.put(vcfbwt::pfp::ENDOFWORD);
    }
    out_ts.put(vcfbwt::pfp::ENDOFDICT);
    
    vcfbwt::DiskWrites::update(out_ts.tellp());
    out_ts.close();
}