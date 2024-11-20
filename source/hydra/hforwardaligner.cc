//*-- AUTHOR : A.Zinchenko <Alexander.Zinchenko@jinr.ru>
//*-- Created : 08/08/2016
//*-- Modified: R. Lalik <Rafal.Lalik@ph.tum.de> 2016/10/10

//_HADES_CLASS_DESCRIPTION
/////////////////////////////////////////////////////////////
//
//  HForwardCandFinder
//
//  Tracking code for Sts
//
/////////////////////////////////////////////////////////////


#include "hforwardaligner.h"
#include "frpcdef.h"
#include "hades.h"
#include "hcategory.h"
#include "hevent.h"
#include "heventheader.h"
#include "hfilter.h"
#include "hforwardcand.h"
#include "hforwardcandfinderpar.h"
#include "hforwardcandsim.h"
#include "hfrpcdetector.h"
#include "hfrpchit.h"
#include "hgeomcompositevolume.h"
#include "hgeomvolume.h"
#include "hphysicsconstants.h"
#include "hruntimedb.h"
#include "hspectrometer.h"
#include "hstscal.h"
#include "hstscalsim.h"
#include "hstsdetector.h"
#include "hstsgeompar.h"

#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TMatrixFLazy.h>

#include <TKDTree.h>

#include <Math/Functor.h>

#include <cassert>
#include <iterator>

#include <promille/promille.hpp>
#include <promille/euler_angles.hpp>

#include "straw_residual_model.hpp"

#include <fmt/core.h>

HForwardAligner::HForwardAligner(const Text_t* name, const Text_t* title, const Text_t* output_bin_file)
: HReconstructor(name, title), fStsCal(nullptr), pFRpcHits(nullptr),
pForwardCand(nullptr), pStrawVFPar(nullptr), isSimulation(kFALSE)
{
    pro_mille = new promille::mille<float>("sts_", output_bin_file);
    straw_planes = new promille::model_planes<float, hsa::sts_residual<float, promille::euler::zyz, 4, 12>>(pro_mille);
}

HForwardAligner::~HForwardAligner()
{
    for (Int_t m = 0; m < STS_MAX_MODULES; ++m)
        for (Int_t l = 0; l < STS_MAX_LAYERS; ++l)
            for (Int_t c = 0; c < STS_MAX_CELLS; ++c)
            {
                delete stsCellsLoc[m][l][c];
                delete stsCellsLab[m][l][c];
            }

    delete pro_mille;
}

Bool_t HForwardAligner::init()
{
    HStsDetector* pSts = (HStsDetector*)(gHades->getSetup()->getDetector("Sts"));
    if (!pSts)
    {
        Error("HForwardAligner::init", "No STS detector found");
        return kFALSE;
    }

    HFRpcDetector* pFRpc = (HFRpcDetector*)(gHades->getSetup()->getDetector("FRpc"));
    if (!pFRpc)
    {
        Error("HForwardAligner::init", "No Forward RPC found");
        return kFALSE;
    }

    // GEANT input data
    pKine = gHades->getCurrentEvent()->getCategory(catGeantKine);
    if (pKine) { isSimulation = kTRUE; }

    // straw hits
    fStsCal = gHades->getCurrentEvent()->getCategory(catStsCal);
    if (!fStsCal)
    {
        Error("HForwardAligner::init()", "HStsCal input missing");
        return kFALSE;
    }

    pForwardCand = gHades->getCurrentEvent()->getCategory(catForwardCand);
    if (!pForwardCand)
    {
        Error("HForwardAligner::init()", "ForwardCand input missing");
        return kFALSE;
    }
    // create the parameter container
    pStsGeomPar = (HStsGeomPar*)gHades->getRuntimeDb()->getContainer("StsGeomPar");
    if (!pStsGeomPar)
    {
        Error("HForwardAligner::init()", "Parameter container for STS geometry not found");
        return kFALSE;
    }

    pCalPar = (HStsCalPar*)gHades->getRuntimeDb()->getContainer("StsCalPar");
    if (!pCalPar) {
        Error("HForwardCandFinderCosy::init()", "Parameter container StsCalPar not found");
        return kFALSE;
    }

    for (Int_t m = 0; m < STS_MAX_MODULES; ++m)
        for (Int_t l = 0; l < STS_MAX_LAYERS; ++l)
            for (Int_t c = 0; c < STS_MAX_CELLS; ++c)
            {
                stsCellsLoc[m][l][c] = nullptr;
                stsCellsLab[m][l][c] = nullptr;
            }

    return kTRUE;
}

Bool_t HForwardAligner::reinit()
{
    std::array<promille::Kind, 12> free_mask_globals[2][4] = {
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

    std::array<promille::Kind, 4> free_mask_locals[2][4] = {
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

    HGeomVector noTrans(0, 0, 0);
    for (Int_t m = 0; m < STS_MAX_MODULES; ++m) {
        for (Int_t l = 0; l < STS_MAX_LAYERS; ++l) {
            HModGeomPar* fmodgeom = pStsGeomPar->getModule(l, m);
            const HGeomTransform& labTrans = fmodgeom->getLabTransform();
            HGeomTransform labRot = labTrans;
            labRot.setTransVector(noTrans);

            std::array<Float_t, 9> rxy;
            for (int iy = 0; iy < 3; ++iy) {
                for (int ix = 0; ix < 3; ++ix) {
                    rxy[iy * 3 + ix] = labTrans.getRotMatrix().getElement(ix, iy);
                }
            }

            auto alpha = atan2(-labTrans.getRotMatrix().getElement(0, 1), labTrans.getRotMatrix().getElement(1, 1));
            auto beta = atan2(labTrans.getRotMatrix().getElement(2, 1), labTrans.getRotMatrix().getElement(2, 2));
            auto gamma = 0.0;

            straw_planes->add_plane((m + 1) * 100 + (l + 1) * 10,
                       1,
                       2,
                       3,
                       5,
                       6,
                       7,
                       pro_mille->add_global_parameter((m + 1) * 100 + (l + 1) * 10 + 1, labTrans.getTransVector().getX(), fmt::format("Tx_c_{:d}{:d}", m+1, l+1)),
                       pro_mille->add_global_parameter((m + 1) * 100 + (l + 1) * 10 + 2, labTrans.getTransVector().getY(), fmt::format("Ty_c_{:d}{:d}", m+1, l+1)),
                       pro_mille->add_global_parameter((m + 1) * 100 + (l + 1) * 10 + 3, labTrans.getTransVector().getZ(), fmt::format("Tz_c_{:d}{:d}", m+1, l+1)),
                       pro_mille->add_global_parameter((m + 1) * 100 + (l + 1) * 10 + 5, alpha, fmt::format("Ra_c_{:d}{:d}", m+1, l+1)),
                       pro_mille->add_global_parameter((m + 1) * 100 + (l + 1) * 10 + 6, beta, fmt::format("Rb_c_{:d}{:d}", m+1, l+1)),
                       pro_mille->add_global_parameter((m + 1) * 100 + (l + 1) * 10 + 7, gamma, fmt::format("Rc_c_{:d}{:d}", m+1, l+1)))
            .set_globals_configuration(free_mask_globals[m][l])
            .set_locals_configuration(free_mask_locals[m][l])
            .set_verbose(false)
            .model()
            .set_local_params(labTrans.getTransVector().getX(),
                              labTrans.getTransVector().getY(),
                              labTrans.getTransVector().getZ(),
                              alpha,
                              beta,
                              gamma);

            if (cosa[m][l] < 0) {
                cosa[m][l] = -cosa[m][l];
                sina[m][l] = -sina[m][l];
            }

            HGeomCompositeVolume* fMod = fmodgeom->getRefVolume();

            for (Int_t c = 0; c < STS_MAX_CELLS; ++c) {
                HGeomVolume* fVol = fMod->getComponent(c);

                if (!fVol)
                    break;

                if (stsCellsLab[m][l][c] == nullptr)
                    stsCellsLab[m][l][c] = new HGeomVector;

                if (stsCellsLoc[m][l][c] == nullptr)
                    stsCellsLoc[m][l][c] = new HGeomVector;

                HGeomVector* p = stsCellsLab[m][l][c];
                *p = fVol->getTransform().getTransVector();
                *stsCellsLoc[m][l][c] = *p;
                *p = labTrans.transFrom(*p);
            }
        }
    }

    return kTRUE;
}

Int_t HForwardAligner::execute()
{
    #ifdef VERBOSE_MODE
    static Int_t event = 0;
    printf("########### EVENT %d\n", event++);
    #endif
    // gErrorIgnoreLevel = kInfo; //debug level
    gErrorIgnoreLevel = kWarning; // debug level

    // // get Event vertex
    // HVertex vertex = gHades->getCurrentEvent()->getHeader()->getVertexReco();
    // if (vertex.getChi2() < 0)
    // { // vertex reco has not been reconstructed: fall back solution : Vertex cluster
    //     vertex = gHades->getCurrentEvent()->getHeader()->getVertexCluster();
    //     with_vertex_cluster = kTRUE;
    // }
    // else { with_vertex_cluster = kFALSE; }
    //
    // primary_vertex = vertex.getPos();
    //
    // // filter best tracks
    // filterTracks();
    //
    // // Select final tracks - remove ghosts
    // selectTracks();
    //
    // // Store tracks
    // storeVectors();


    return 0;
}
