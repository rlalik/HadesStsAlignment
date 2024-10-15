#include "lib.hpp"

// Hydra
#include <hades.h>
#include <hloop.h>
#include <hparasciifileio.h>
#include <hparrootfileio.h>
#include <hruntimedb.h>
#include <hspectrometer.h>
#include <hstsdetector.h>
#include <hstsgeompar.h>

// System
#include <iostream>
#include <string>

#include <getopt.h>

auto main(int argc, char* argv[]) -> int
{
    int verbose = 0;
    std::string root_par_file;
    std::string ascii_par_file;
    std::string output_file = "NewStsAlignment.txt";

    int c;
    while (1) {
        static struct option long_options[] = {/* These options set a flag. */
                                               // {"verbose", no_argument, &verbose, 1},
                                               // {"brief", no_argument, &verbose, 0},

                                               /* These options don’t set a flag.
                                                *              We distinguish them by their indices. */
                                               {"ascii", required_argument, 0, 'a'},
                                               {"root", required_argument, 0, 'r'},
                                               {"output", required_argument, 0, 'o'},
                                               {"verbose", no_argument, 0, 'v'},
                                               {0, 0, 0, 0}};
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long(argc, argv, "a:r:o:v:", long_options, &option_index);

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
                ascii_par_file = optarg;
                break;

            case 'r':
                root_par_file = optarg;
                break;

            case 'o':
                output_file = optarg;
                break;

            case 'v':
                verbose = atoi(optarg);
                printf("SET VERBOSE LEVEL %d\n", verbose);
                break;

            case '?':
                /* getopt_long already printed an error message. */
                break;

            default:
                abort();
        }
    }

    HLoop* loop = new HLoop(kTRUE);

    std::string input_results;

    /* Print any remaining command line arguments (not options). */
    if (optind < argc) {
        Bool_t ret;
        //              printf ("non-option ARGV-elements: ");
        while (optind < argc) {
            input_results = argv[optind++];
            printf("INPUT: %s\n", input_results.c_str());
            break;
        }
    } else {
        std::cerr << "ERROR : no input results given\n";
        std::exit(EXIT_FAILURE);
    }

    HGeomVector stsCellsLab_correction[STS_MAX_MODULES][STS_MAX_LAYERS];

    ifstream infile;
    infile.open(input_results);
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    if (!infile.is_open()) {
        cout << "Error reading new params" << endl;
        exit(1);
    }

    int param_id;
    float value, error, dummy1, dummy2;

    while (infile >> param_id >> value >> error >> dummy1 >> dummy2) {
        int vlay = param_id / 10;
        int mod = vlay / 4;
        int lay = vlay % 4;
        int param = param_id % 10;

        switch (param) {
            case 1:
                stsCellsLab_correction[mod][lay].X() = value;
                break;
            case 2:
                stsCellsLab_correction[mod][lay].Y() = value;
                break;
            case 3:
                stsCellsLab_correction[mod][lay].Z() = value;
                break;
            case 4:
                stsCellsLab_correction[mod][lay];
                break;
            case 5:
                stsCellsLab_correction[mod][lay];
                break;
            case 6:
                stsCellsLab_correction[mod][lay];
                break;
            default:
                abort();
        }
    }

    Hades* myHades = new Hades;
    HRuntimeDb* rtdb = gHades->getRuntimeDb();

    HSpectrometer* spec = gHades->getSetup();

    // create the detector and its setup
    // and add it in the spectrometer
    HStsDetector* det = new HStsDetector;

    Int_t stsMods[STS_MAX_LAYERS][STS_MAX_MODULES] = {{1, 1}, {1, 1}, {1, 1}, {1, 1}};

    spec->addDetector(new HStsDetector);
    for (int i = 0; i < STS_MAX_LAYERS; ++i)
        spec->getDetector("Sts")->setModules(i, stsMods[i]);

    auto input_cnt = 0;
    // Checking if parameters file was opened properly
    if (ascii_par_file.length()) {
        HParAsciiFileIo* input1 = new HParAsciiFileIo;
        input1->open((Text_t*)ascii_par_file.c_str(), "in");
        if (!input1->check()) {
            std::cerr << "Param file " << ascii_par_file << " not open!" << std::endl;
            // abort();
        } else {
            rtdb->setFirstInput(input1);
            input_cnt++;
        }
    }

    // Checking if parameters file was opened properly
    if (root_par_file.length()) {
        HParRootFileIo* input2 = new HParRootFileIo;
        input2->open((Text_t*)root_par_file.c_str(), "READ");

        if (!input2->check()) {
            std::cerr << "Param file " << root_par_file << " not open!" << std::endl;
            // abort();
        } else {
            if (input_cnt)
                rtdb->setSecondInput(input2);
            else
                rtdb->setFirstInput(input2);
        }
    }

    // Oracle as output
    HParAsciiFileIo* out = new HParAsciiFileIo;
    out->open(output_file.c_str(), "out");
    rtdb->setOutput(out);

    // create the parameter containers
    HStsGeomPar* pStrawGeomPar = (HStsGeomPar*)gHades->getRuntimeDb()->getContainer("StsGeomPar");

    rtdb->initContainers(444571079);  // listed in test db (sep10test)

    for (Int_t m = 0; m < STS_MAX_MODULES; ++m) {
        for (Int_t l = 0; l < STS_MAX_LAYERS; ++l) {
            HModGeomPar* fmodgeom = pStrawGeomPar->getModule(l, m);
            HGeomTransform& labTrans = fmodgeom->getLabTransform();

            printf("\n*** MOD  %d  LAY  %d ***\n", m, l);
            printf("    INPUT ALIGNMENT\n");
            labTrans.print();
            printf("    CORRECTIONS\n");
            stsCellsLab_correction[m][l].print();
            printf("    FINAL\n");
            auto old_trans = labTrans.getTransVector();
            old_trans += stsCellsLab_correction[m][l];
            labTrans.setTransVector(old_trans);
            labTrans.print();
        }
    }

    rtdb->saveOutput();
    // rtdb->print();

    delete myHades;

    return 0;
}
