#pragma once

#include "straw_residual_model.hpp"

// promille
#include <promille/euler_angles.hpp>
#include <promille/promille.hpp>

// Hydra
#include <forwarddef.h>
#include <frpcdef.h>

// ROOT
#include <TH1.h>
#include <TH2.h>

// system
#include <string>

#include <stsdef.h>

class HLoop;
class HCategory;
class HGeomVector;
class HStsGeomPar;
class HStsCalPar;

namespace hsa
{
// clang-format off
    /*
     * Normal parameters with a single dimension straw marked as free
     */
    static std::array<promille::Kind, 12> free_mask_globals[2][4] = {
        {
            {
                promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
                promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED
            }, {
                promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
                promille::Kind::FIXED, promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED
            }, {
                promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
                promille::Kind::FIXED, promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED
            }, {
                promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
                promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED
            },
        },{
            {
                promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
                promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED
            }, {
                promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
                promille::Kind::FIXED, promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,promille::Kind::FIXED
            }, {
                promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
                promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED
            }, {
                promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
                promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED
            },
        }};
// clang-format on

// // clang-format off
// std::array<promille::Kind, 12> free_mask_globals[2][4] = {
//     {
//         {
//             promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
//             promille::Kind::FIXED, promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
//             promille::Kind::FIXED, promille::Kind::FIXED
//         }, {
//             promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
//             promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
//             promille::Kind::FIXED, promille::Kind::FIXED
//         }, {
//             promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
//             promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
//             promille::Kind::FIXED, promille::Kind::FIXED
//         }, {
//             promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
//             promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
//             promille::Kind::FIXED, promille::Kind::FIXED
//         },
//     },{
//         {
//             promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
//             promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
//             promille::Kind::FIXED, promille::Kind::FIXED
//         }, {
//             promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
//             promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
//             promille::Kind::FIXED, promille::Kind::FIXED
//         }, {
//             promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
//             promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
//             promille::Kind::FIXED, promille::Kind::FIXED
//         }, {
//             promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
//             promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED,
//             promille::Kind::FIXED, promille::Kind::FIXED
//         },
//     }};
// // clang-format on

static std::array<promille::Kind, 4> free_mask_locals[2][4] = {
    {
        {promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FREE, promille::Kind::FIXED},
        {promille::Kind::FIXED, promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FREE},
        {promille::Kind::FIXED, promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FREE},
        {promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FREE, promille::Kind::FIXED},
    },
    {
        {promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FREE, promille::Kind::FIXED},
        {promille::Kind::FIXED, promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FREE},
        {promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE},
        {promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE},
    }};

// std::array<promille::Kind, 4> free_mask_locals[2][4] = {
//     {
//         {promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED},
//         {promille::Kind::FIXED, promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED},
//         {promille::Kind::FIXED, promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED},
//         {promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED},
//     },
//     {
//         {promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED},
//         {promille::Kind::FIXED, promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED},
//         {promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED},
//         {promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FIXED},
//     }};

// std::array<promille::Kind, 4> free_mask_locals[2][4] = {
//     {
//         {promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE},
//         {promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE},
//         {promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE},
//         {promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE},
//     },
//     {
//         {promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE},
//         {promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE},
//         {promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE},
//         {promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE},
//     }};

struct aligner_qa;

/**
 * @brief The core implementation of the executable
 *
 * This class makes up the library part of the executable, which means that the
 * main logic is implemented here. This kind of separation makes it easy to
 * test the implementation for the executable, because the logic is nicely
 * separated from the command-line logic implemented in the main function.
 */
struct forward_aligner_library
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
    forward_aligner_library(HLoop* loop, const std::string& output_file, const std::string& log_file_name);
    ~forward_aligner_library();

    auto set_qa_file(const char* qafile) -> void;

    auto standalone_init(const std::string& root_par, const std::string& ascii_par) -> bool;
    auto library_init() -> bool;

    auto model_init(std::array<promille::Kind, 12> free_mask_globals[2][4], std::array<promille::Kind, 4> free_mask_locals[2][4]) -> bool;

    auto find_beam_avgs(long long nevents) -> HGeomVector;
    auto execute(long long events = 0, long long first = 0) -> void;

    auto process_event(unsigned long event) -> int;

    // run paramaters
    float sigma {0};
    bool project {false};
    bool beam_tilt {false};
    bool all_tracks {false};

    // members
    std::string output_file;
    std::ofstream log_file;

    HLoop* loop;

    HCategory* fGeantKine;
    HCategory* fStsCal;
    HCategory* fFRpcClus;
    HCategory* fForwardCand;

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
    TH2I BeamXY;
    TH1I BeamZ;

    promille::mille<float> pro_mille;
    decltype(pro_mille.make_model_planes<hsa::sts_residual<float, promille::euler::zyz>>()) straw_planes;

    static int verbose;
    aligner_qa* qa {nullptr};
    std::array<long, 10> rejected = {0};
    bool sim {false};
};

};  // namespace hsa
