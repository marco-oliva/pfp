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
    CLI::App app("VCF to Fasta");
    
    std::vector<std::string> vcfs_file_names;
    std::vector<std::string> refs_file_names;
    std::string out_file;
    std::size_t max_samples = 0;
    std::size_t threads = 1;
    std::string samples_file_name;
    std::string haplotype_string = "1";

    vcfbwt::pfp::Params params;
    
    app.add_option("-v,--vcf", vcfs_file_names, "List of comma ',' separated vcf files. Assuming in genome order!")->allow_extra_args(true)->configurable()->delimiter(',');
    app.add_option("-r,--ref", refs_file_names, "List of comma ',' separated reference files. Assuming in genome order!")->allow_extra_args(true)->configurable()->delimiter(',');
    app.add_option("-H,--haplotype", haplotype_string, "Haplotype: [1,2,12].")->configurable();
    app.add_option("-o,--out-file", out_file, "Output prefix")->configurable();
    app.add_option("-m, --max", max_samples, "Max number of samples to analyze")->configurable();
    app.add_option("-S, --samples", samples_file_name, "File containing the list of samples to parse")->configurable();
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
    vcfbwt::VCF vcf(refs_file_names, vcfs_file_names, samples_file_name, max_samples);

    // Generate fasta file, reference first
    std::ofstream samples(out_file);
    std::string reference;
    reference.insert(reference.end(), vcf.get_reference().begin(), vcf.get_reference().end());
    samples << "> Reference" << std::endl;
    samples.write(reference.c_str(), reference.size());
    samples.put('\n');
    
    
    for (std::size_t i = 0; i < vcf.size(); i++)
    {
        if (haplotype_string == "1" or haplotype_string == "2")
        {
            std::size_t sample_genotype;
            if (haplotype_string == "1") { sample_genotype = 0; }
            else { sample_genotype = 1; }
    
            vcfbwt::Sample::iterator it(vcf[i], sample_genotype);
            std::string sample;
    
            while (not it.end()) { sample.push_back(*it); ++it; }
            samples << "> " + vcf[i].id() + "\n";
            samples.write(sample.c_str(), sample.size());
            samples.put('\n');
        }
        else
        {
            // first haplotype
            vcfbwt::Sample::iterator it_h1(vcf[i], 0);
            std::string sample_h1;
            while (not it_h1.end()) { sample_h1.push_back(*it_h1); ++it_h1; }
            samples << "> " + vcf[i].id() + "H1 \n";
            samples.write(sample_h1.c_str(), sample_h1.size());
            samples.put('\n');
    
            // second haplotype
            vcfbwt::Sample::iterator it_h2(vcf[i], 0);
            std::string sample_h2;
            while (not it_h2.end()) { sample_h2.push_back(*it_h2); ++it_h2; }
            samples << "> " + vcf[i].id() + "H2 \n";
            samples.write(sample_h2.c_str(), sample_h2.size());
            samples.put('\n');
        }
    }
    samples.close();
    
    return 0;
    
}
