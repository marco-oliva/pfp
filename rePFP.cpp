//
//  rePFP.cpp
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
    std::string out_prefix;
    std::string tmp_dir;
    std::size_t max_samples = 0;
    std::size_t threads = 1;
    std::size_t threshold = 0;

    vcfbwt::pfp::Params params;

    app.add_option("-v,--vcf", vcfs_file_names, "List of vcf files. Assuming in genome order!")->allow_extra_args(true)->expected(-1)->configurable();
    app.add_option("-r,--ref", refs_file_names, "List of reference files. Assuming in genome order!")->allow_extra_args(true)->expected(-1)->configurable();
    app.add_option("-o,--out-prefix", out_prefix, "Output prefix")->configurable();
    app.add_option("-m, --max", max_samples, "Max number of samples to analyze")->configurable();
    app.add_option("-w, --window-size", params.w, "Sliding window size")->check(CLI::Range(3, 200))->configurable();
    app.add_option("-p, --modulo", params.p, "Module used during parsing")->check(CLI::Range(50, 20000))->configurable();
    app.add_option("-t, --threads", threads, "Number of threads")->configurable();
    app.add_option("-T, --threshold", threshold, "Trigger Strings cost threshold")->configurable();
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

    // Set tmp file dir
    if (tmp_dir != "") { vcfbwt::TempFile::setDirectory(tmp_dir); }

    // ------------------------------------------------------------
    // First pass of pfp
    // ------------------------------------------------------------

    // Parse the VCF
    vcfbwt::VCF vcf(refs_file_names, vcfs_file_names, max_samples);

    // Set threads accordingly to configuration
    omp_set_num_threads(threads);

    vcfbwt::pfp::ReferenceParse reference_parse(vcf.get_reference(), params);

    vcfbwt::pfp::Parser main_parser(params, out_prefix, reference_parse);

    std::vector<vcfbwt::pfp::Parser> workers(threads);
    for (std::size_t i = 0; i < workers.size(); i++)
    {
        std::size_t tag = vcfbwt::pfp::Parser::WORKER | vcfbwt::pfp::Parser::UNCOMPRESSED;
        if (i == workers.size() - 1) { tag = tag | vcfbwt::pfp::Parser::LAST; }
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

    // ------------------------------------------------------------
    // AuPair
    // ------------------------------------------------------------
    vcfbwt::pfp::AuPair au_pair_algo(out_prefix, params.w);

    std::set<std::string_view> removed_trigger_strings;
    int removed_bytes = au_pair_algo.compress(removed_trigger_strings, threshold);

    spdlog::info("Removed {} bytes", removed_bytes);

    std::ofstream out_ts(out_prefix + ".removed_ts");
    for (auto& rts : removed_trigger_strings)
    {
        for (auto& c : rts) { out_ts.put(c); }
        out_ts.put(vcfbwt::pfp::ENDOFWORD);
    }
    out_ts.put(vcfbwt::pfp::ENDOFDICT);

    vcfbwt::DiskWrites::update(out_ts.tellp());
    out_ts.close();

    // ------------------------------------------------------------
    // Second pass of PFP
    // ------------------------------------------------------------

    params.ignore_ts_file = out_prefix + ".removed_ts";
    vcfbwt::pfp::ReferenceParse reference_parse_aupair(vcf.get_reference(), params);

    vcfbwt::pfp::Parser main_parser_aupair(params, out_prefix + "_compressed", reference_parse_aupair);

    std::vector<vcfbwt::pfp::Parser> workers_aupair(threads);
    for (std::size_t i = 0; i < workers.size(); i++)
    {
        std::size_t tag = vcfbwt::pfp::Parser::WORKER | vcfbwt::pfp::Parser::UNCOMPRESSED;
        if (i == workers_aupair.size() - 1) { tag = tag | vcfbwt::pfp::Parser::LAST; }
        workers_aupair[i].init(params, "", reference_parse_aupair, tag);
        main_parser_aupair.register_worker(workers_aupair[i]);
    }

    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < vcf.size(); i++)
    {
        int this_thread = omp_get_thread_num();
        spdlog::info("rePFP Processing sample [{}/{}]: {}", i, vcf.size(), vcf[i].id());
        malloc_count_print_status();

        workers_aupair[this_thread](vcf[i]);
    }

    // close the main parser and exit
    main_parser_aupair.close();

}
