//
//  vcf_to_fa.cpp
//
//  Copyright 2020 Marco Oliva. All rights reserved.
//

#include <CLI/CLI.hpp>
#include <version.hpp>
#include <utils.hpp>
#include <vcf.hpp>
#include <pfp_algo.hpp>

int main(int argc, char **argv)
{
    CLI::App app("PFP++");
    
    std::vector<std::string> vcfs_file_names;
    std::vector<std::string> refs_file_names;
    std::string out_file;
    std::size_t max_samples = 0;
    std::size_t threads = 1;
    
    vcfbwt::pfp::Params params;
    
    app.add_option("-v,--vcf", vcfs_file_names, "List of vcf files. Assuming in genome order!")->allow_extra_args(true)->expected(-1)->configurable();
    app.add_option("-r,--ref", refs_file_names, "List of reference files. Assuming in genome order!")->allow_extra_args(true)->expected(-1)->configurable();
    app.add_option("-o,--out-file", out_file, "Output prefix")->configurable();
    app.add_option("-m, --max", max_samples, "Max number of samples to analyze")->configurable();
    app.add_option("-t, --threads", threads, "Number of threads")->configurable();
    app.add_flag_callback("--version",vcfbwt::Version::print,"Version");
    app.set_config("--configure");
    app.allow_windows_style_options();
    
    CLI11_PARSE(app, argc, argv);
    
    // Clean file name vectors
    vcfs_file_names.erase(std::remove_if(vcfs_file_names.begin(), vcfs_file_names.end(),
                                         [] (std::string& s) {return s.size() == 0; } ), vcfs_file_names.end());
    refs_file_names.erase(std::remove_if(refs_file_names.begin(), refs_file_names.end(),
                                         [] (std::string& s) {return s.size() == 0; } ), refs_file_names.end());
    
    // Print out configurations
    spdlog::info("Current Configuration:\n{}", app.config_to_str(true,true));
    
    // Parse the VCF
    vcfbwt::VCF vcf(refs_file_names, vcfs_file_names);
    if (max_samples != 0) { vcf.set_max_samples(max_samples); }
    
    // Generate fasta file, reference first
    std::ofstream samples(out_file);
    std::string reference;
    reference.insert(reference.end(), vcf.get_reference().begin(), vcf.get_reference().end());
    samples << "> Reference" << std::endl;
    samples.write(reference.c_str(), reference.size());
    samples.put('\n');

    for (std::size_t i = 0; i < vcf.size(); i++)
    {
        vcfbwt::Sample::iterator it(vcf[i]);
        std::string sample;
        
        while (not it.end()) { sample.push_back(*it); ++it; }
        samples << "> " + vcf[i].id() + "\n";
        samples.write(sample.c_str(), sample.size());
        samples.put('\n');
    }
    samples.close();
    
    return 0;
    
}
