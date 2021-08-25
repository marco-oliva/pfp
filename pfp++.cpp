//
//  pfp++.cpp
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
    //spdlog::set_pattern("%v");
    
    CLI::App app("PFP++");
    
    std::vector<std::string> vcfs_file_names;
    std::vector<std::string> refs_file_names;
    std::string fasta_file_path;
    std::string text_file_path;
    std::string out_prefix;
    std::string tmp_dir;
    std::size_t max_samples = 0;
    std::size_t threads = 1;
    bool check = false;
    bool only_trigger_strings = false;
    
    vcfbwt::pfp::Params params;
    
    app.add_option("-v,--vcf", vcfs_file_names, "List of vcf files. Assuming in genome order!")->allow_extra_args(true)->configurable();
    app.add_option("-r,--ref", refs_file_names, "List of reference files. Assuming in genome order!")->allow_extra_args(true)->configurable();
    app.add_option("-f,--fasta", fasta_file_path, "Fasta file to parse.")->configurable()->check(CLI::ExistingFile);
    app.add_option("-t,--text", text_file_path, "Text file to parse.")->configurable()->check(CLI::ExistingFile);
    app.add_option("-o,--out-prefix", out_prefix, "Output prefix")->configurable();
    app.add_option("-m, --max", max_samples, "Max number of samples to analyze")->configurable();
    app.add_option("-w, --window-size", params.w, "Sliding window size")->check(CLI::Range(3, 200))->configurable();
    app.add_option("-p, --modulo", params.p, "Module used during parisng")->check(CLI::Range(5, 20000))->configurable();
    app.add_option("-j, --threads", threads, "Number of threads")->configurable();
    app.add_option("--ignore-ts", params.ignore_ts_file, "Ignore Trigger Strings in file")->configurable();
    app.add_option("--tmp-dir", tmp_dir, "Tmp file directory")->check(CLI::ExistingDirectory)->configurable();
    app.add_flag("-c, --compression", params.compress_dictionary, "Compress the dictionary")->configurable();
    app.add_flag("--use-acceleration", params.use_acceleration, "Use reference parse to avoid re-parsing")->configurable();
    app.add_flag("--print-statistics", params.print_out_statistics_csv, "Print out csv containing stats")->configurable();
    app.add_flag("--occurrences", params.compute_occurrences, "Compute the .occ file")->configurable();
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
    
    if (not fasta_file_path.empty())
    {
        if (out_prefix.empty()) { spdlog::error("If parsing fasta -o,--out-prefix required"); return; }
        vcfbwt::pfp::ParserFasta main_parser(params, fasta_file_path, out_prefix);
    
        // Run
        main_parser();
    
        // Close the main parser
        main_parser.close();
    }
    else if (not text_file_path.empty())
    {
        if (out_prefix.empty()) { spdlog::error("If parsing text -o,--out-prefix required"); return; }
        vcfbwt::pfp::ParserText main_parser(params, text_file_path, out_prefix);
    
        // Run
        main_parser();
    
        // Close the main parser
        main_parser.close();
    }
    else
    {
        // Set tmp file dir
        if (tmp_dir != "") { vcfbwt::TempFile::setDirectory(tmp_dir); }
    
        // Parse the VCF
        vcfbwt::VCF vcf(refs_file_names, vcfs_file_names, max_samples);
    
        // Set threads accordingly to configuration
        omp_set_num_threads(threads);
    
        vcfbwt::pfp::ReferenceParse reference_parse(vcf.get_reference(), params);
    
        vcfbwt::pfp::ParserVCF main_parser(params, out_prefix, reference_parse);
    
        std::vector<vcfbwt::pfp::ParserVCF> workers(threads);
        for (std::size_t i = 0; i < workers.size(); i++)
        {
            std::size_t tag = vcfbwt::pfp::ParserVCF::WORKER | vcfbwt::pfp::ParserVCF::UNCOMPRESSED;
            if (i == workers.size() - 1) { tag = tag | vcfbwt::pfp::ParserVCF::LAST; }
            workers[i].init(params, "", reference_parse, tag);
            main_parser.register_worker(workers[i]);
        }

        #pragma omp parallel for schedule(static)
        for (std::size_t i = 0; i < vcf.size(); i++)
        {
            int this_thread = omp_get_thread_num();
            spdlog::info("Processing sample [{}/{}]: {}", i, vcf.size(), vcf[i].id());
            malloc_count_print_status();
        
            workers[this_thread](vcf[i]);
        }
    
        // close the main parser and exit
        main_parser.close();
    }
}
