//
//  merge.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#include <CLI/CLI.hpp>
#include <version.hpp>
#include <utils.hpp>
#include <pfp_algo.hpp>

int main(int argc, char **argv)
{
    CLI::App app("PFP++ from fasta");
    
    std::string left_file;
    std::string right_file;
    std::string out_prefix;
    
    vcfbwt::pfp::Params params;
    
    app.add_option("-l,--left-prefix", left_file, "Left Prefix")->configurable()->required();
    app.add_option("-r,--right-prefix", right_file, "Right Prefix")->configurable()->required();
    app.add_option("-o,--out-prefix", out_prefix, "Output prefix")->configurable()->required();
    app.add_option("-w, --window-size", params.w, "Sliding window size")->check(CLI::Range(0, 100))->configurable();
    app.add_option("-p, --module", params.p, "Module used during parisng")->check(CLI::Range(0, 1000))->configurable();
    app.add_flag_callback("--version",vcfbwt::Version::print,"Version");
    app.set_config("--configure");
    app.allow_windows_style_options();
    
    CLI11_PARSE(app, argc, argv);
    
    // Print out configurations
    spdlog::info("Current Configuration:\n{}", app.config_to_str(true,true));
    
    // Merge parsings
    vcfbwt::pfp::ParserUtils<char>::merge(left_file, right_file, out_prefix, params);
    
    return 0;
}
