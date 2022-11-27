//
//  exprop.cpp
//
//  Copyright 2022 Marco Oliva. All rights reserved.
//

#include <CLI/CLI.hpp>
#include <version.hpp>
#include <utils.hpp>
#include <pfp_algo.hpp>

int main(int argc, char **argv)
{
    CLI::App app("Extract properties from PFPs");
    
    std::string pfp_prefix;
    bool integers_pfp = false;
    
    vcfbwt::pfp::Params params;
    
    app.add_option("--pfp-prefix", pfp_prefix, "PFP Prefix")->configurable()->required();
    app.add_option("-w, --window-size", params.w, "Sliding window size")->check(CLI::Range(0, 100))->configurable();
    app.add_option("-p, --module", params.p, "Module used during parisng")->check(CLI::Range(0, 1000))->configurable();
    app.add_flag("--output-occurrences", params.output_occurrences, "Output count for each dictionary phrase.")->configurable();
    app.add_flag("--output-sai", params.output_sai, "Output sai array.")->configurable();
    app.add_flag("--output-last", params.output_last, "Output last array.")->configurable();
    app.add_flag("--output-compressed-dict", params.compress_dictionary, "Output compressed dictionary.")->configurable();
    app.add_flag("--integers", integers_pfp, "Integer (uint32_t) PFP");
    app.add_flag_callback("--version",vcfbwt::Version::print,"Version");
    app.set_config("--configure");
    app.allow_windows_style_options();
    
    CLI11_PARSE(app, argc, argv);
    
    // Print out configurations
    spdlog::info("Current Configuration:\n{}", app.config_to_str(true,true));
    
    if (not integers_pfp)
    {
        vcfbwt::pfp::PropertiesWriter<vcfbwt::char_type> properties_out(pfp_prefix, params);
        properties_out.write();
    }
    else
    {
        vcfbwt::pfp::PropertiesWriter<uint32_t> properties_out(pfp_prefix, params);
        properties_out.write();
    }
    
    return 0;
    
    return 0;
}