// internal
#include "../source/straw_residual_model.hpp"

// external
#include <promille/promille.hpp>

// ROOT
#include <Math/Point3D.h>
#include <Math/Vector3D.h>

// system
#include <string>
#include <getopt.h>

using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;

auto main(int argc, char* argv[]) -> int
{
    /* Flag set by ‘--verbose’. */
    int verbose_flag;

    int c;

    while (1) {
        static struct option long_options[] = {/* These options set a flag. */
                                               {"verbose", no_argument, &verbose_flag, 1},
                                               {"brief", no_argument, &verbose_flag, 0},
                                               /* These options don’t set a flag.
                                                  We distinguish them by their indices. */
                                               // {"add", no_argument, 0, 'a'},
                                               // {"append", no_argument, 0, 'b'},
                                               // {"delete", required_argument, 0, 'd'},
                                               // {"create", required_argument, 0, 'c'},
                                               // {"file", required_argument, 0, 'f'},
                                               {0, 0, 0, 0}};
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long(argc, argv, "abc:d:f:", long_options, &option_index);

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

            case 'a':
                puts("option -a\n");
                break;

            case 'b':
                puts("option -b\n");
                break;

            case 'c':
                printf("option -c with value `%s'\n", optarg);
                break;

            case 'd':
                printf("option -d with value `%s'\n", optarg);
                break;

            case 'f':
                printf("option -f with value `%s'\n", optarg);
                break;

            case '?':
                /* getopt_long already printed an error message. */
                break;

            default:
                abort();
        }
    }

    /* Instead of reporting ‘--verbose’
       and ‘--brief’ as they are encountered,
       we report the final status resulting from them. */
    if (verbose_flag)
        puts("verbose flag is set");

    std::string input_data_file;

    /* Print any remaining command line arguments (not options). */
    if (optind < argc) {
        printf("non-option ARGV-elements: ");
        while (optind < argc) {
            printf("%s ", argv[optind]);
            input_data_file = std::string(argv[optind++]);
            break;
        }
        putchar('\n');
    }

    promille::promille<hsa::sts_residual<float, promille::euler::zyz>> mille("test_", "test");

    mille.set_local_parameter(0, "Bx");
    mille.set_local_parameter(1, "By");
    mille.set_local_parameter(2, "Tx");
    mille.set_local_parameter(3, "Ty");

    mille.add_global_parameter(1, -10, "Txg1");
    mille.add_global_parameter(2, -20, "Tyg1");
    mille.add_global_parameter(3, -30, "Tzg1");
    mille.add_global_parameter(5, -40, "Rag1");
    mille.add_global_parameter(6, -50, "Rbg1");
    mille.add_global_parameter(7, -60, "Rcg1");

    mille.add_global_parameter(16, -10, "Tz1");
    mille.add_global_parameter(15, -20, "Ty1");
    mille.add_global_parameter(14, -30, "Tx1");
    mille.add_global_parameter(11, -40, "Ra1");
    mille.add_global_parameter(12, -50, "Rb1");
    mille.add_global_parameter(13, -60, "Rc1");

    // mille.add_globals_set(0, 1, 2, 3, 4, 5, 6, 10, 20, 30, 40, 50, 60);
    mille.add_plane(0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16)
        .set_globals_configuration(promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FREE)
        .set_locals_configuration(promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FREE);

    mille.print(true);

    mille.add_measurement(0, 0, XYZPoint(0, 0, 1), XYZVector(2, 3, 4));

    /*
        auto mille = promille::mille_builder<promille::euler::zyz>("test_", "hades_basic_example_alignment_data.bin");
        mille.set_verbose(verbose_flag);

        if (input_data_file.length()) {
            mille.read_global_data(input_data_file);
        } else {
            mille.add_globals_set(0,
                                  {1, 0, promille::Kind::FREE},
                                  {2, 0, promille::Kind::FIXED},
                                  {3, 0, promille::Kind::FIXED},
                                  {4, float(TMath::DegToRad()) * 0, promille::Kind::FIXED},
                                  {5, 0, promille::Kind::FIXED},
                                  {6, 0, promille::Kind::FIXED},
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE);

            mille.add_globals_set(1,
                                  {1, 0, promille::Kind::FIXED},
                                  {2, 0, promille::Kind::FREE},
                                  {3, 0, promille::Kind::FIXED},
                                  {4, float(TMath::DegToRad()) * 90, promille::Kind::FIXED},
                                  {5, 0, promille::Kind::FIXED},
                                  {6, 0, promille::Kind::FIXED},
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE);

            mille.add_globals_set(2,
                                  {1, 0, promille::Kind::FIXED},
                                  {2, 0, promille::Kind::FREE},
                                  {3, 0, promille::Kind::FIXED},
                                  {4, float(TMath::DegToRad()) * 90, promille::Kind::FIXED},
                                  {5, 0, promille::Kind::FIXED},
                                  {6, 0, promille::Kind::FIXED},
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE);

            mille.add_globals_set(3,
                                  {1, 0, promille::Kind::FREE},
                                  {2, 0, promille::Kind::FIXED},
                                  {3, 0, promille::Kind::FIXED},
                                  {4, float(TMath::DegToRad()) * 0, promille::Kind::FIXED},
                                  {5, 0, promille::Kind::FIXED},
                                  {6, 0, promille::Kind::FIXED},
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE);

            mille.add_globals_set(4,
                                  {1, 0, promille::Kind::FREE},
                                  {2, 0, promille::Kind::FIXED},
                                  {3, 0, promille::Kind::FIXED},
                                  {4, float(TMath::DegToRad()) * 0, promille::Kind::FIXED},
                                  {5, 0, promille::Kind::FIXED},
                                  {6, 0, promille::Kind::FIXED},
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE);

            mille.add_globals_set(5,
                                  {1, 0, promille::Kind::FIXED},
                                  {2, 0, promille::Kind::FREE},
                                  {3, 0, promille::Kind::FIXED},
                                  {4, float(TMath::DegToRad()) * 90, promille::Kind::FIXED},
                                  {5, 0, promille::Kind::FIXED},
                                  {6, 0, promille::Kind::FIXED},
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE);

            mille.add_globals_set(6,
                                  {1, 0, promille::Kind::FREE},
                                  {2, 0, promille::Kind::FREE},
                                  {3, 0, promille::Kind::FIXED},
                                  {4, float(TMath::DegToRad()) * 45, promille::Kind::FIXED},
                                  {5, 0, promille::Kind::FIXED},
                                  {6, 0, promille::Kind::FIXED},
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE);

            mille.add_globals_set(7,
                                  {1, 0, promille::Kind::FREE},
                                  {2, 0, promille::Kind::FREE},
                                  {3, 0, promille::Kind::FIXED},
                                  {4, float(TMath::DegToRad()) * -45, promille::Kind::FIXED},
                                  {5, 0, promille::Kind::FIXED},
                                  {6, 0, promille::Kind::FIXED},
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE,
                                  promille::Kind::FREE);
        }

        mille.print(true);

        mille.add_measurement(0, 0, 1, 2, 3, 4, 0, 0, 0.1, 0.15);
        mille.add_measurement(1, 1, 2, 3, 4, 5, 0, 0, 0.1, 0.15);
        mille.add_measurement(2, 2, 3, 4, 5, 6, 0, 0, 0.1, 0.15);
        mille.add_measurement(3, 0, 1, 2, 3, 4, 0, 0, 0.1, 0.15);
        mille.add_measurement(4, 1, 2, 3, 4, 5, 0, 0, 0.1, 0.15);
        mille.add_measurement(5, 2, 3, 4, 5, 6, 0, 0, 0.1, 0.15);
        mille.end();

        mille.add_measurement(0, 1, 2, 3, 4, 5, 0, 0, 0.1, 0.15);
        mille.add_measurement(1, 2, 3, 4, 5, 6, 0, 0, 0.1, 0.15);
        mille.add_measurement(2, 3, 4, 5, 6, 7, 0, 0, 0.1, 0.15);
        mille.end();
    */
    return 0;
}
