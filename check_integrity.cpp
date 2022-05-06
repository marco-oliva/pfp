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


template <typename data_type>
void check(const std::string& input_dict_path, const std::string& input_parse_path, const std::string& input_occurrences_path, std::size_t window_size)
{
    spdlog::info("Reading dictionary");
    std::vector<std::vector<data_type>> dict;
    vcfbwt::pfp::ParserUtils<data_type>::read_dictionary(input_dict_path, dict);
    spdlog::info("Dictionary phrases: {}", dict.size());

    spdlog::info("Reading parse");
    std::vector<vcfbwt::size_type> parse;
    vcfbwt::pfp::ParserUtils<data_type>::read_parse(input_parse_path, parse);
    spdlog::info("Parse length: {}", parse.size());

    if (parse[0] != data_type(1)) { spdlog::error("parse[0] != 1"); exit(EXIT_FAILURE); }

    spdlog::info("Checking PFP");
    std::size_t errors = 0;
    for (std::size_t i = 1; i < parse.size(); i++)
    {
        if (parse[i] > dict.size()) { spdlog::error("parse[{}] > dict.size()", i); exit(EXIT_FAILURE); }
        if (parse[i] == 0) { spdlog::error("parse[{}] == 0", i); exit(EXIT_FAILURE); }

        // Check Trigger Strings
        std::vector<data_type> ts_prev(dict[parse[i-1] - 1].end() - window_size, dict[parse[i-1] - 1].end());
        std::vector<data_type> ts_curr(dict[parse[i] - 1].begin(), dict[parse[i] - 1].begin() + window_size);

        if (ts_prev != ts_curr)
        {
            errors += 1;
            std::stringstream previous, current;
            std::copy(dict[parse[i-1]-1].begin(), dict[parse[i-1]-1].end(), std::ostream_iterator<data_type>(previous, " "));
            std::copy(dict[parse[i]-1].begin(), dict[parse[i]-1].end(), std::ostream_iterator<data_type>(current, " "));
            spdlog::error("Missmatching trigger strings:\nP[{}]: {}\nP[{}] {}", i-1, previous.str(), i, current.str());
        }
    }
    spdlog::info("Reached end of parse. {} errors found.", errors);

    if (not input_occurrences_path.empty())
    {
        spdlog::info("Check occurrences");
        std::vector<vcfbwt::long_type> occ_computed(dict.size(), 0);
        for (auto& p : parse) { occ_computed[p - 1] += 1; }

        // Check the occ file
        bool occ_good = true;

        std::ifstream occ_stream(input_occurrences_path);
        if (parse.size() < std::numeric_limits<vcfbwt::short_type>::max())
        {
            spdlog::info("Occurrences type: vcfbwt::short_type");
            std::vector<vcfbwt::short_type> occ(dict.size(), 0);
            occ_stream.read((char*) occ.data(), sizeof(vcfbwt::short_type) * occ.size());

            for (std::size_t i = 0; i < occ.size(); i++)
            {
                if (occ[i] != occ_computed[i]) { spdlog::error("OccM[{}] = {}\tOccD[{}] = {}", i, occ_computed[i], i, occ[i]); }
                occ_good = occ_good and (occ[i] == occ_computed[i]);
            }
        }
        else
        {
            spdlog::info("Occurrences type: vcfbwt::long_type");
            std::vector<vcfbwt::long_type> occ(dict.size(), 0);
            occ_stream.read((char*) occ.data(), sizeof(vcfbwt::long_type) * occ.size());

            for (std::size_t i = 0; i < occ.size(); i++)
            {
                if (occ[i] != occ_computed[i]) { spdlog::error("OccM[{}] = {}\tOccD[{}] = {}", i, occ_computed[i], i, occ[i]); }
                occ_good = occ_good and (occ[i] == occ_computed[i]);
            }
        }
        if (occ_good) { spdlog::info("Occurrences file ok"); }
        else { spdlog::error("Error in occurrences file"); }
    }
}


int main(int argc, char **argv)
{
    CLI::App app("Check Integrity of PFP");

    std::string input_dict_path;
    std::string input_parse_path;
    std::string input_occurrences_path;
    std::size_t window_size;
    bool integers_pfp = false;

    app.add_option("-d,--dictionary", input_dict_path, "Dictionary")->required();
    app.add_option("-p,--parse", input_parse_path, "Parse")->required();
    app.add_option("-o,--occurences", input_occurrences_path, "Occurrences");
    app.add_flag("--integers", integers_pfp, "Integer (uint32_t) PFP");
    app.add_option("-w, --window", window_size, "Window size")->required();
    app.add_flag_callback("--version",vcfbwt::Version::print,"Version");
    app.allow_windows_style_options();

    CLI11_PARSE(app, argc, argv);

    // Print out configurations
    spdlog::info("Current Configuration:\n{}", app.config_to_str(true,true));

    // Unparse
    if (not integers_pfp)
    {
        check<vcfbwt::char_type>(input_dict_path, input_parse_path, input_occurrences_path, window_size);
    }
    else
    {
        check<uint32_t>(input_dict_path, input_parse_path, input_occurrences_path, window_size);
    }
}