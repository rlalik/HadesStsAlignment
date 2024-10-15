#include "lib.hpp"

// Hydra
#include <hloop.h>

// System
#include <filesystem>
#include <iostream>
#include <string>

#include <getopt.h>

namespace fs = std::filesystem;

bool replace(std::string& str, const std::string& from, const std::string& to)
{
    size_t start_pos = str.find(from);
    if (start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

auto main(int argc, char* argv[]) -> int
{
    int verbose {0};
    int events {0};
    int project {0};
    int beam_tilt {0};
    int all_tracks {0};
    std::string root_par_file;
    std::string ascii_par_file = "feb22_dst_params.txt";
    std::string output_file = "sts_alignment.bin";
    std::string log_file_name;
    std::string qa_file;
    float sigma {0.3};

    int c;
    while (1) {
        static struct option long_options[] = {/* These options set a flag. */
                                               // {"verbose", no_argument, &verbose, 1},
                                               // {"brief", no_argument, &verbose, 0},
                                               {"project", no_argument, &project, 1},
                                               {"tilt", no_argument, &beam_tilt, 1},
                                               {"all", no_argument, &all_tracks, 1},

                                               /* These options don’t set a flag.
                                                *              We distinguish them by their indices. */
                                               {"events", required_argument, 0, 'e'},
                                               {"ascii", required_argument, 0, 'a'},
                                               {"log", required_argument, 0, 'l'},
                                               {"root", required_argument, 0, 'r'},
                                               {"sigma", required_argument, 0, 's'},
                                               {"output", required_argument, 0, 'o'},
                                               {"qafile", required_argument, 0, 'q'},
                                               {"dump-param", required_argument, 0, 'p'},
                                               {"verbose", no_argument, 0, 'v'},
                                               {0, 0, 0, 0}};
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long(argc, argv, "e:a:l::r:s:vo:q:", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c) {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
                printf("option %s", long_options[option_index].name);
                if (optarg)
                    printf(" with arg %s", optarg);
                printf("\n");
                break;

            case 'e':
                events = atol(optarg);
                break;

            case 'a':
                ascii_par_file = optarg;
                break;

            case 'l':
                log_file_name = (optarg == nullptr) ? "sts_alignment.log" : optarg;
                break;

            case 'r':
                root_par_file = optarg;
                break;

            case 's':
                sigma = atof(optarg);
                break;

            case 'o':
                output_file = optarg;
                break;

            case 'q':
                qa_file = optarg;
                break;

            case 'v':
                verbose = true;  // atoi(optarg);
                printf("SET VERBOSE LEVEL %d\n", verbose);
                break;

            case '?':
                /* getopt_long already printed an error message. */
                break;

            default:
                abort();
        }
    }

    /* Instead of reporting ‘--verbose’
     *      and ‘--brief’ as they are encountered,
     *           we report the final status resulting from them. */
    // if (anapars.verbose)
    // puts("verbose flag is set");

    HLoop* loop = new HLoop(kTRUE);

    std::string first_file;

    /* Print any remaining command line arguments (not options). */
    if (optind < argc) {
        Bool_t ret;
        //              printf ("non-option ARGV-elements: ");
        while (optind < argc) {
            TString infile = argv[optind++];
            printf("INPUT: %s\n", infile.Data());
            if (infile.Contains(","))
                ret = loop->addMultFiles(infile);
            else if (infile.Contains(".root")) {
                if (!first_file.length())
                    first_file = infile;
                ret = loop->addFiles(infile);
            } else
                ret = loop->addFilesList(infile);

            if (!ret) {
                std::cerr << "READBACK: ERROR : cannot find inputfiles : " << infile.Data() << endl;
                std::exit(EXIT_FAILURE);
            }
        }
    }

    library::verbose = verbose;
    library lib(loop, output_file, root_par_file, ascii_par_file, log_file_name);

    lib.qa_file = qa_file;
    replace(lib.qa_file, "%f", fs::path(first_file).filename());

    lib.sigma = sigma;
    lib.project = project;
    lib.beam_tilt = beam_tilt;
    lib.all_tracks = all_tracks;

    lib.execute(events);

    return 0;
}
