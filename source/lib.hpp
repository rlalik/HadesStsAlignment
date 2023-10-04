#pragma once

#include <string>

#include <Rtypes.h>
#include <TH1.h>
#include <forwarddef.h>
#include <frpcdef.h>
#include <stsdef.h>

#include "mille_builder/mille_builder.hpp"

class HLoop;
class HCategory;
class HGeomVector;
class HStsGeomPar;
class HStsCalPar;

/**
 * @brief The core implementation of the executable
 *
 * This class makes up the library part of the executable, which means that the
 * main logic is implemented here. This kind of separation makes it easy to
 * test the implementation for the executable, because the logic is nicely
 * separated from the command-line logic implemented in the main function.
 */
struct library
{
    enum ElasticsErrorCode
    {
        OK,
        NOCAT,
        NOTUSED,
        BADMULT,
        NOPHIACC,
        NOTHETAACC,
        NOACC
    };

    /**
     * @brief Simply initializes the name member to the name of the project
     */
    library(HLoop* loop, const std::string& output_file, const std::string& root_par, const std::string& ascii_par);

    auto check_elastics_hf(Float_t phi_diff_min, Float_t phi_diff_max, Float_t thetap_diff_min, Float_t thetap_diff_max)
        -> std::tuple<ElasticsErrorCode, float, float>;

    auto check_elastics_hf(Float_t phi1,
                           Float_t theta1,
                           Float_t phi2,
                           Float_t theta2,
                           Float_t phi_diff_min,
                           Float_t phi_diff_max,
                           Float_t thetap_diff_min,
                           Float_t thetap_diff_max) -> std::tuple<ElasticsErrorCode, float, float>;

    auto find_beam_avgs(long long nevents) -> HGeomVector;
    auto execute(long long events = 0) -> void;

    // run paramaters
    std::string qa_file;
    float sigma {0};
    bool hack {false};
    bool beam_tilt {false};

    // members
    std::string output_file;
    HLoop* loop;

    HCategory* fGeantKine;
    HCategory* fStsCal;
    HCategory* fFRpcClus;
    HCategory* fForwardCand;
    HCategory* fParticleCand;

    HStsGeomPar* pStrawGeomPar;
    HStsCalPar* pCalPar;  // cal run par

    Int_t nLayers[STS_MAX_MODULES];
    Float_t cosa[STS_MAX_MODULES][STS_MAX_LAYERS];
    Float_t sina[STS_MAX_MODULES][STS_MAX_LAYERS];
    Float_t roty[STS_MAX_MODULES][STS_MAX_LAYERS];

    HGeomVector* stsCellsLab[STS_MAX_MODULES][STS_MAX_LAYERS][STS_MAX_CELLS];  // Centre of the strip [mm]
    HGeomVector* stsCellsLoc[STS_MAX_MODULES][STS_MAX_LAYERS][STS_MAX_CELLS];  // Centre of the strip [mm]

    TH1I BeamX;
    TH1I BeamY;
    TH1I BeamZ;

    mb::mille_builder<mb::euler::zyz> mille;

    static int verbose;
};
