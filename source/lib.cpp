#include "lib.hpp"

// Hydra
#include <hcategorymanager.h>
// #include <hdetgeompar.h>
#include <hdst.h>
#include <hforwardcand.h>
#include <hforwardtools.h>
// #include <hfrpccalpar.h>
// #include <hfrpccluster.h>
#include <hfrpcdetector.h>
// #include <hfrpchit.h>
#include <hgeomcompositevolume.h>
// #include <hgeomvector.h>
#include <hgeomvolume.h>
#include <hloop.h>
#include <hparasciifileio.h>
#include <hparrootfileio.h>
#include <hparticlebeamtilt.h>
#include <hparticlecand.h>
#include <hparticletool.h>
#include <hphysicsconstants.h>
#include <hruntimedb.h>
#include <hspectrometer.h>
#include <hstscal.h>
#include <hstscalpar.h>
#include <hstsdetector.h>
#include <hstsgeompar.h>

// ROOT
#include <Math/EulerAngles.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TFile.h>
#include <TGButton.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TRootEmbeddedCanvas.h>
#include <TStyle.h>
#include <TText.h>

// System
#include <iostream>
#include <vector>

using std::cout;
using std::endl;

int library::verbose = 0;

library::library(HLoop* loop,
                 const std::string& output_file,
                 const std::string& root_par,
                 const std::string& ascii_par,
                 const std::string& log_file_name)
    : loop {loop}
    , pro_mille("sts_", output_file.data())
    , straw_planes(pro_mille.make_model_planes<hsa::sts_residual<float, promille::euler::zyz>>())
{
    pro_mille.set_verbose(verbose);

    BeamX = TH1I("Beam_X", "beam (X)", 100, -20, 20);
    BeamY = TH1I("Beam_Y", "beam (Y)", 100, -20, 20);
    BeamXY = TH2I("Beam_XY", "beam (XY)", 100, -20, 20, 100, -20, 20);
    BeamZ = TH1I("Beam_Z", "beam (Z)", 200, -200, 00);

    if (!loop->setInput("-*,+HParticleCand,+HForwardCand,+HStsCal,+HFRpcHit,+HStart2Hit")) {  // reading file structure
        std::cerr << "READBACK: ERROR : cannot read input !" << std::endl;
        abort();
    }

    // clang-format: off
    Int_t mdcMods[6][4] = {{1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}};
    // clang-format: on
    HDst::setupSpectrometer("FEB22", mdcMods, "sts,frpc");
    HSpectrometer* spec = gHades->getSetup();
    Int_t stsMods[STS_MAX_LAYERS][STS_MAX_MODULES] = {{1, 1}, {1, 1}, {1, 1}, {1, 1}};
    Int_t frpcMods[FRPC_MAX_SECTORS][FRPC_MAX_MODULES] = {{1}, {1}, {1}, {1}};

    spec->addDetector(new HStsDetector);
    for (int i = 0; i < STS_MAX_LAYERS; ++i)
        spec->getDetector("Sts")->setModules(i, stsMods[i]);

    spec->addDetector(new HFRpcDetector);

    for (int i = 0; i < FRPC_MAX_SECTORS; ++i)
        spec->getDetector("FRpc")->setModules(i, frpcMods[i]);

    // TString asciiParFile = "feb21_dst_params.txt"; // File containing all parameters of the detector
    HRuntimeDb* rtdb = HRuntimeDb::instance();  // myHades -> getRuntimeDb();

    // Checking if parameters file was opened properly
    if (ascii_par.length()) {
        HParAsciiFileIo* input1 = new HParAsciiFileIo;
        input1->open((Text_t*)ascii_par.c_str(), "in");
        if (!input1->check()) {
            std::cerr << "Param file " << ascii_par << " not open!" << std::endl;
            // abort();
        } else {
            rtdb->setFirstInput(input1);
        }
    }

    // Checking if parameters file was opened properly
    if (root_par.length()) {
        HParRootFileIo* input2 = new HParRootFileIo;
        input2->open((Text_t*)root_par.c_str(), "READ");

        if (!input2->check()) {
            std::cerr << "Param file " << root_par << " not open!" << std::endl;
            // abort();
        } else {
            rtdb->setSecondInput(input2);
        }
    }

    // Initializing the parameter container for geometry
    pStrawGeomPar = (HStsGeomPar*)gHades->getRuntimeDb()->getContainer("StsGeomPar");
    // Checking if parameters container for geometry was created properly
    if (!pStrawGeomPar) {
        Error("HStsDigitizer::init()", "Parameter container for geometry not created");
        abort();
    }

    pCalPar = (HStsCalPar*)gHades->getRuntimeDb()->getContainer("StsCalPar");
    if (!pCalPar) {
        Error("HForwardCandFinderCosy::init()", "Parameter container StsCalPar not found");
        abort();
    }

    if (!rtdb->initContainers(444571079))
        if (!rtdb->initContainers(17000))
            rtdb->initContainers(1);
    // rtdb->print();

    for (Int_t m = 0; m < STS_MAX_MODULES; ++m)
        for (Int_t l = 0; l < STS_MAX_LAYERS; ++l)
            for (Int_t c = 0; c < STS_MAX_CELLS; ++c) {
                stsCellsLab[m][l][c] = nullptr;
                stsCellsLoc[m][l][c] = nullptr;
            }

    // clang-format off
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
    //             promille::Kind::FIXED,promille::Kind::FIXED
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

    // std::array<promille::Kind, 4> free_mask_locals[2][4] = {
    //     {
    //         {promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FREE, promille::Kind::FIXED},
    //         {promille::Kind::FIXED, promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FREE},
    //         {promille::Kind::FIXED, promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FREE},
    //         {promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FREE, promille::Kind::FIXED},
    //     },
    //     {
    //         {promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FREE, promille::Kind::FIXED},
    //         {promille::Kind::FIXED, promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FREE},
    //         {promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE},
    //         {promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE},
    //     }};

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

    std::array<promille::Kind, 4> free_mask_locals[2][4] = {
        {
            {promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE},
            {promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE},
            {promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE},
            {promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE},
        },
        {
            {promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE},
            {promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE},
            {promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE},
            {promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE},
        }};

    pro_mille.add_global_parameter(1, 0, "Txg1");
    pro_mille.add_global_parameter(2, 0, "Tyg1");
    pro_mille.add_global_parameter(3, 0, "Tzg1");
    pro_mille.add_global_parameter(5, 0, "Rag1");
    pro_mille.add_global_parameter(6, 0, "Rbg1");
    pro_mille.add_global_parameter(7, 0, "Rcg1");

    // elastics scattering plane
    // mille.add_plane<hsa::sts_residual<float, promille::euler::zyz>>(0, 1, 2, 3, 5, 6, 7, 11, 12, 13, 15, 16, 17)
    //     .set_globals_configuration(promille::Kind::FREE,
    //                                promille::Kind::FREE,
    //                                promille::Kind::FIXED,
    //                                promille::Kind::FIXED,
    //                                promille::Kind::FIXED,
    //                                promille::Kind::FIXED,
    //                                promille::Kind::FIXED,
    //                                promille::Kind::FIXED,
    //                                promille::Kind::FIXED,
    //                                promille::Kind::FIXED,
    //                                promille::Kind::FIXED,
    //                                promille::Kind::FIXED)
    //     .set_locals_configuration(promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE)
    //     .get_model()->set_local_params(0, 0, 0, 0, 0, 0);

    // straws planes
    straw_planes.set_verbose(verbose);

    HGeomVector noTrans(0, 0, 0);
    for (Int_t m = 0; m < STS_MAX_MODULES; ++m) {
        for (Int_t l = 0; l < STS_MAX_LAYERS; ++l) {
            HModGeomPar* fmodgeom = pStrawGeomPar->getModule(l, m);
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

            if (verbose) {
                printf("\n");
                printf("%f   %f   %f\n",
                       labTrans.getRotMatrix().getElement(0, 0),
                       labTrans.getRotMatrix().getElement(0, 1),
                       labTrans.getRotMatrix().getElement(0, 2));
                printf("%f   %f   %f\n",
                       labTrans.getRotMatrix().getElement(1, 0),
                       labTrans.getRotMatrix().getElement(1, 1),
                       labTrans.getRotMatrix().getElement(1, 2));
                printf("%f   %f   %f\n",
                       labTrans.getRotMatrix().getElement(2, 0),
                       labTrans.getRotMatrix().getElement(2, 1),
                       labTrans.getRotMatrix().getElement(2, 2));
                printf("\nANGLES= %f %f %f\n", alpha * TMath::RadToDeg(), beta * TMath::RadToDeg(), gamma * TMath::RadToDeg());
                printf("TRANS = %f %f %f\n",
                       labTrans.getTransVector().getX(),
                       labTrans.getTransVector().getY(),
                       labTrans.getTransVector().getZ());
            }

            straw_planes
                .add_plane((m + 1) * 100 + (l + 1) * 10,
                           1,
                           2,
                           3,
                           5,
                           6,
                           7,
                           pro_mille.add_global_parameter((m + 1) * 100 + (l + 1) * 10 + 1, 0, "Txc"),
                           pro_mille.add_global_parameter((m + 1) * 100 + (l + 1) * 10 + 2, 0, "Tyc"),
                           pro_mille.add_global_parameter((m + 1) * 100 + (l + 1) * 10 + 3, 0, "Tzc"),
                           pro_mille.add_global_parameter((m + 1) * 100 + (l + 1) * 10 + 5, 0, "Rac"),
                           pro_mille.add_global_parameter((m + 1) * 100 + (l + 1) * 10 + 6, 0, "Rbc"),
                           pro_mille.add_global_parameter((m + 1) * 100 + (l + 1) * 10 + 7, 0, "Rcc"))
                .set_globals_configuration(free_mask_globals[m][l])
                .set_locals_configuration(free_mask_locals[m][l])
                .set_verbose(verbose)
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

    pro_mille.write_param_file();

    fGeantKine = HCategoryManager::getCategory(catGeantKine, kTRUE, "catGeantKine");

    fStsCal = HCategoryManager::getCategory(catStsCal, kTRUE, "catStsCal");
    if (!fStsCal) {
        cout << "No catStsCal!" << endl;
    }

    // Forward Detector RPC Hits
    fFRpcClus = HCategoryManager::getCategory(catFRpcClus, kTRUE, "catFRpcCluster");
    if (!fFRpcClus) {
        cout << "No catFRpcClus!" << endl;
    }

    // Forward Detector Candidates
    fForwardCand = HCategoryManager::getCategory(catForwardCand, kTRUE, "catForwardCand");
    if (!fForwardCand) {
        cout << "No catForwardCand!" << endl;
    }

    // Forward Detector Candidates
    fParticleCand = HCategoryManager::getCategory(catParticleCand, kTRUE, "catParticleCand");
    if (!fParticleCand) {
        cout << "No catParticleCand!" << endl;
    }

    // Checking if ParticleCand and ForwardCand are open
    if (!fParticleCand or !fForwardCand)
        abort();

    if (log_file_name.length()) {
        log_file = std::ofstream(log_file_name, std::ios::out);
    }
}

auto library::check_elastics_hf(Float_t phi_diff_min, Float_t phi_diff_max, Float_t thetap_diff_min, Float_t thetap_diff_max)
    -> std::tuple<ElasticsErrorCode, float, float>
{
    auto particle_cand_cnt = fParticleCand->getEntries();
    auto forward_cand_cnt = fForwardCand->getEntries();

    // Checking for 1 track in HADES and 1 track in FwDet
    if (particle_cand_cnt < 1 or forward_cand_cnt < 1)
        return {BADMULT, 0, 0};

    // Particle 1 - HADES Particle
    HParticleCand* cand1 = dynamic_cast<HParticleCand*>(fParticleCand->getObject(0));
    // Particle 2 - FwDet Particle
    HForwardCand* cand2 = dynamic_cast<HForwardCand*>(fForwardCand->getObject(0));

    if (!cand2->isFlagBit(HForward::kIsUsed) or !cand1->isFlagBit(Particle::kIsUsed)) {
        return {NOTUSED, 0, 0};
    }

    const auto t_phi_P1 = cand1->getPhi();  // [deg]
    const auto t_theta_P1 = cand1->getTheta();  // [deg]

    // Correction for track start position (Z_START)
    //     auto t_z_start_P2 = cand1->getZ();          // [mm]
    //     FIXME not valid any more
    //     cand2->setStartXYZ(0.0, 0.0, t_z_start_P2); // We take the track start position from the
    // reconstructed HADES track paired with FwDet track
    //     cand2->calcPoints(); // Recalculating ToF and Distance to account for correction
    // cand2->calc4vectorProperties(HPhysicsConstants::mass(14));  // Assume proton mass

    const auto t_phi_P2 = cand2->getPhi();  // [deg]
    const auto t_theta_P2 = cand2->getTheta();  // [deg]

    // Two Particles Variables
    const auto min_phi = min(t_phi_P1, t_phi_P2);
    const auto max_phi = max(t_phi_P1, t_phi_P2);
    auto t_phi_diff = max_phi - min_phi;  // [deg]
    // auto t_phi_diff = abs(t_phi_P1 - t_phi_P2);  // [deg]
    if (t_phi_diff > 180.)
        t_phi_diff = 360. - t_phi_diff;

    const auto t_tan_theta_product = TMath::Tan(TMath::DegToRad() * t_theta_P1) * TMath::Tan(TMath::DegToRad() * t_theta_P2);  // [ ]

    const auto phi_acc = ((t_phi_diff >= phi_diff_min) && (t_phi_diff <= phi_diff_max));
    const auto theta_acc = ((t_tan_theta_product >= thetap_diff_min) && (t_tan_theta_product <= thetap_diff_max));

    if (phi_acc && theta_acc)
        return {OK, t_phi_diff, t_tan_theta_product};
    if (phi_acc)
        return {NOTHETAACC, t_phi_diff, t_tan_theta_product};
    if (theta_acc)
        return {NOPHIACC, t_phi_diff, t_tan_theta_product};
    return {NOACC, t_phi_diff, t_tan_theta_product};
}

auto library::check_elastics_hf(Float_t phi1,
                                Float_t theta1,
                                Float_t phi2,
                                Float_t theta2,
                                Float_t phi_diff_min,
                                Float_t phi_diff_max,
                                Float_t thetap_diff_min,
                                Float_t thetap_diff_max) -> std::tuple<ElasticsErrorCode, float, float>
{
    // Correction for track start position (Z_START)
    //     auto t_z_start_P2 = cand1->getZ();          // [mm]
    //     FIXME not valid any more
    //     cand2->setStartXYZ(0.0, 0.0, t_z_start_P2); // We take the track start position from the
    // reconstructed HADES track paired with FwDet track
    //     cand2->calcPoints(); // Recalculating ToF and Distance to account for correction
    // cand2->calc4vectorProperties(HPhysicsConstants::mass(14));  // Assume proton mass

    // Two Particles Variables
    const auto min_phi = min(phi1, phi2);
    const auto max_phi = max(phi1, phi2);
    auto t_phi_diff = max_phi - min_phi;  // [deg]
    // auto t_phi_diff = abs(t_phi_P1 - t_phi_P2);  // [deg]
    if (t_phi_diff > 180.)
        t_phi_diff = 360. - t_phi_diff;

    const auto t_tan_theta_product = TMath::Tan(TMath::DegToRad() * theta1) * TMath::Tan(TMath::DegToRad() * theta2);  // [ ]

    const auto phi_acc = ((t_phi_diff >= phi_diff_min) && (t_phi_diff <= phi_diff_max));
    const auto theta_acc = ((t_tan_theta_product >= thetap_diff_min) && (t_tan_theta_product <= thetap_diff_max));

    if (phi_acc && theta_acc)
        return {OK, t_phi_diff, t_tan_theta_product};
    if (phi_acc)
        return {NOTHETAACC, t_phi_diff, t_tan_theta_product};
    if (theta_acc)
        return {NOPHIACC, t_phi_diff, t_tan_theta_product};
    return {NOACC, t_phi_diff, t_tan_theta_product};
}

auto library::execute(long long events) -> void
{
    // Popup the GUI...
    gStyle->SetOptStat(1);

    TH2I* h_vert_reco_xy {nullptr};
    TH1I* h_vert_reco_z {nullptr};

    TH1I* h_fcand_chi2_ndf {nullptr};

    TH1I* h_delta_phi {nullptr};
    TH1I* h_tan_prod {nullptr};
    TH1I* h_elastic_error {nullptr};
    TH2I* h_elastic_mult {nullptr};

    TH2I* h_fcand_txty {nullptr};
    TH2I* h_fcand_theta_phi {nullptr};
    TH2I* h_fproj_txty {nullptr};
    TH2I* h_fproj_vert {nullptr};
    TH2I* h_el_theta_phi {nullptr};
    TH2I* h_el_theta_phi_rejected {nullptr};
    TH2I* h_el_theta_phi_all {nullptr};
    TH2I* h_proj_diff_theta_phi {nullptr};
    TH1I* h_p1_p2_dist {nullptr};
    TH1I* h_p1_beam_dist {nullptr};

    TH1I* h_qa_sign {nullptr};
    TH1I* h_qa_t_d_dist {nullptr};
    TH1I* h_qa_drift_dist {nullptr};
    TH1I* h_qa_residuals {nullptr};
    TH1I* h_qa_uresiduals {nullptr};

    TH2I* h_qa_sign_plane {nullptr};
    TH2I* h_qa_t_d_dist_plane {nullptr};
    TH2I* h_qa_drift_dist_plane {nullptr};
    TH2I* h_qa_residuals_plane {nullptr};
    TH2I* h_qa_uresiduals_plane {nullptr};

    if (qa_file.length()) {
        h_vert_reco_xy = new TH2I("h_vert_reco_xy", "h_vert_reco_xy", 100, -20, 20, 400, -20, 20);
        h_vert_reco_z = new TH1I("h_vert_reco_z", "h_vert_reco_z", 200, -200, 0);

        h_fcand_chi2_ndf = new TH1I("h_fcand_chi2_ndf", "h_fcand_chi2_ndf", 50, 0, 5);

        h_delta_phi = new TH1I("h_delta_phi", "h_delta_phi", 200, 170, 190);
        h_tan_prod = new TH1I("h_tan_prod", "h_tan_prod", 400, 0, 1);
        h_elastic_error = new TH1I("h_elastic_error", "h_elastic_error", 7, 0, 7);
        h_elastic_mult = new TH2I("h_elastic_mult", "h_elastic_mult", 10, 0, 10, 10, 0, 10);
        // TH2I* h_fcand_delta = new TH2I("h_fcand_delta", "h_fcand_delta", 200, -50, 50, 200, -50, 50);
        h_fcand_txty = new TH2I("h_fcand_txty", "h_fcand_txty;tx;ty", 200, -0.2, 0.2, 200, -0.2, 0.2);
        h_fcand_theta_phi = new TH2I("h_fcand_theta_phi", "h_fcand_theta_phi;phi [deg];theta [deg]", 360, 0, 360, 200, 0, 10);
        h_fproj_txty = new TH2I("h_fproj_txty", "h_fproj_txty;tx;ty", 200, -0.2, 0.2, 200, -0.2, 0.2);
        h_fproj_vert = new TH2I("h_fproj_vert", "h_fproj_vert;x [mm];y [mm]", 200, -0.2, 0.2, 200, -0.2, 0.2);
        h_el_theta_phi = new TH2I("h_el_theta_phi", "h_el_theta_phi;phi [deg];theta [deg]", 100, 175, 180, 100, 0.2, 0.4);
        h_el_theta_phi_rejected =
            new TH2I("h_el_theta_phi_rejected", "h_el_theta_phi_rejected;phi [deg];theta [deg]", 180, 0, 180, 100, 0, 1);
        h_el_theta_phi_all = new TH2I("h_el_theta_phi_all", "h_el_theta_phi_all;phi [deg];theta [deg]", 180, 0, 180, 100, 0, 1);
        h_proj_diff_theta_phi =
            new TH2I("h_proj_diff_theta_phi", "h_proj_diff_theta_phi;#Delta phi [deg];#Delta theta [deg]", 100, -4, 4, 100, -2, 2);
        h_p1_p2_dist = new TH1I("h_p1_p2_dist", "h_p1_p2_dist;dist [mm];counts", 100, 0, 15);
        h_p1_beam_dist = new TH1I("h_p1_beam_dist", "h_p1_beam_dist;dist [mm];counts", 100, 0, 15);

        h_qa_sign = new TH1I("h_qa_sign", "h_qa_sign;sign;counts", 2, -1, 2);
        h_qa_t_d_dist = new TH1I("h_qa_t_d_dist", "h_qa_t_d_dist;track-wire distance [mm];counts", 100, 0, 10);
        h_qa_drift_dist = new TH1I("h_qa_drift_dist", "h_qa_drift_dist;track-wire distance [mm];counts", 100, 0, 10);
        h_qa_residuals = new TH1I("h_qa_residuals", "h_qa_residuals;residuals [mm];counts", 100, -1, 1);
        h_qa_uresiduals = new TH1I("h_qa_uresiduals", "h_qa_uresiduals;residuals [mm];counts", 100, -1, 1);

        h_qa_sign_plane = new TH2I("h_qa_sign_plane", "h_qa_sign_plane;sign;plane;counts", 3, -1, 2, 8, 1, 9);
        h_qa_t_d_dist_plane =
            new TH2I("h_qa_t_d_dist_plane", "h_qa_t_d_dist_plane;track-wire distance [mm];plane;counts", 100, 0, 10, 8, 1, 9);
        h_qa_drift_dist_plane =
            new TH2I("h_qa_drift_dist_plane", "h_qa_drift_dist_plane;track-wire distance [mm];plane;counts", 100, 0, 10, 8, 1, 9);
        h_qa_residuals_plane = new TH2I("h_qa_residuals_plane", "h_qa_residuals_plane;residuals [mm];plane;counts", 100, -1, 1, 8, 1, 9);
        h_qa_uresiduals_plane = new TH2I("h_qa_uresiduals_plane", "h_qa_uresiduals_plane;residuals [mm];plane;counts", 100, -1, 1, 8, 1, 9);
    }

    auto nevts = events == 0 ? loop->getEntries() : min(events, loop->getEntries());

    bool sim = (fGeantKine != nullptr);

    if (sim) {
        printf("*** SIM MODE ***\n");
    } else {
        printf("*** EXP MODE ***\n");
    }

    auto beam_avg = find_beam_avgs(10000);
    beam_avg.print();
    if (sim) {
        beam_avg = HGeomVector(0, 0, 0);
        beam_avg.print();
    }

    printf("   Beam AVG = %f  %f    %f\n", beam_avg.X(), beam_avg.Y(), beam_avg.Z());

    const Double_t pi_value = 3.14159265359;  // Pi constant value
    const Double_t c_value = 299.792458;  // Speed of light (In mm/ns)
    const Double_t proton_mass = 938.272081358;  // Proton mass [MeV/c^2]

    const Double_t z_target = -115.0;  // Position of target (middle point) [mm]
    const Double_t beam_Ekin = 4500;  // Kinetic energy of beam particle [MeV]
    const Double_t beam_E01 = beam_Ekin + proton_mass;  // Total energy of beam particle [MeV]
    const Double_t beam_gamma_cm = sqrt((beam_E01 + proton_mass) / (2 * proton_mass));  // Gamma of Center of Mass system
    const Double_t beam_p01 = sqrt(pow(beam_E01, 2) - pow(proton_mass, 2));  // Momentum of beam particle [MeV/c]
    const Double_t inverse_gamma2 = 1.0 / (beam_gamma_cm * beam_gamma_cm);  // 0.29429;

    pro_mille.print();

    long good_vectors_count = 0;

    std::array<long, 10> rejected = {0};

    if (!fParticleCand)
        abort();
    if (!fForwardCand)
        abort();

    const bool use_beam_tilt = beam_tilt && !sim;

    std::cout << "============= CONFIG =============\n"
              << " Beam tilt: " << use_beam_tilt;
    if (beam_tilt && sim) {
        std::cout << " (tilt enabled but disabled by sim)";
    }
    std::cout << '\n';
    std::cout << " Project: " << project << "\n Events to analyze: " << nevts << "\n============= CONFIG =============\n";

    for (int event = 0l; event < nevts; ++event) {
        loop->nextEvent(event);
        if (event % 10000 == 0)
            cout << event << "  " << 100.0 * event / nevts << "%" << endl;

        auto particle_cand_cnt = fParticleCand->getEntries();
        auto forward_cand_cnt = fForwardCand->getEntries();

        // Checking for 1 track in HADES and 1 track in FwDet
        if (particle_cand_cnt != 1 or forward_cand_cnt != 1) {
            // BADMULT
            rejected[2]++;
            continue;
        }

        // Particle 1 - HADES Particle
        // HParticleCand* cand1 = dynamic_cast<HParticleCand*>(fParticleCand->getObject(0));
        // Particle 2 - FwDet Particle
        // HForwardCand* cand2 = dynamic_cast<HForwardCand*>(fForwardCand->getObject(0));

        auto p_cand = dynamic_cast<HParticleCand*>(fParticleCand->getObject(0));
        auto f_cand = dynamic_cast<HForwardCand*>(fForwardCand->getObject(0));

        if (!f_cand->isFlagBit(HForward::kIsUsed) or !p_cand->isFlagBit(Particle::kIsUsed)) {
            // NOTUSED
            rejected[1]++;
            continue;
        }

        /* BEAM TILT */
        HEventHeader* evHeader = loop->getEventHeader();
        HParticleBeamTilt beamTilt;
        beamTilt.setEventHeader(evHeader);

        Double_t thDeg, phiDeg;  // angles after beam tilt correction
        beamTilt.correctAngles(p_cand->getTheta(), p_cand->getPhi(), thDeg, phiDeg);  // "cand" is HParticleCand

        const TLorentzVector* beam = beamTilt.beamLVector();
        TVector3 tilted_dir = beam->Vect();  //(beam->X(), beam->Y(), beam->Z());
        tilted_dir *= 1. / tilted_dir.Z();

        float el_phi_min, el_phi_max, el_theta_min, el_theta_max;
        if (sim) {
            el_phi_min = 179.;
            el_phi_max = 180.;
            el_theta_min = inverse_gamma2 - 0.05;
            el_theta_max = inverse_gamma2 + 0.05;
        } else {
            el_phi_min = 177.;
            el_phi_max = 180.;
            el_theta_min = inverse_gamma2 - 0.05;
            el_theta_max = inverse_gamma2 + 0.05;
        }

        // check conditions for elastics scattering
        auto [res, phi_diff, theta_tan_prod] =
            check_elastics_hf(phiDeg, thDeg, f_cand->getPhi(), f_cand->getTheta(), el_phi_min, el_phi_max, el_theta_min, el_theta_max);

        if (qa_file.length()) {
            h_elastic_error->Fill(res);
            h_elastic_mult->Fill(fParticleCand->getEntries(), fForwardCand->getEntries());
            h_el_theta_phi_all->Fill(phi_diff, theta_tan_prod);
        }

        if (res != OK) {
            if (res == NOCAT)
                rejected[0]++;
            // else if (res == NOTUSED) rejected[1]++;
            // else if (res == BADMULT) rejected[2]++;
            else {
                if (res == NOPHIACC)
                    rejected[3]++;
                if (res == NOTHETAACC)
                    rejected[4]++;
                if (res == NOACC)
                    rejected[5]++;
                if (qa_file.length()) {
                    h_el_theta_phi_rejected->Fill(phi_diff, theta_tan_prod);
                }
            }
            continue;
        }

        // printf("Angles tilt correction: phi  %f -> %f   theta %f -> %f\n", p_cand->getPhi(), phiDeg, p_cand->getTheta(), thDeg);
        // beam->Print();
        // tilted_dir.Print();

        // if (sim and fcand->getChi2()/fcand->getNDF() > 0.8) continue;
        // if (!sim and fcand->getChi2()/fcand->getNDF() > 2.0) continue;

        p_cand->calc4vectorProperties(HPhysicsConstants::mass(14));
        f_cand->calc4vectorProperties(HPhysicsConstants::mass(14));

        /* Alignment steps:
         * 1. Find p_H, p_F
         * 2. Check the distance between these two tracks < CUT
         * 3. Check the distance |p_H - beam_avg| < CUT
         * 4. b_F = POCA(p_H, p_F)
         *
         * 5. if (project):
         *      p_F' = p_beam - p_H:
         *    else
         *      p_F' = p_F,
         *
         * 6. a. beam_tilt==1 - use tilted beam
         *    b. beam_tilt==0 - use (0,0,1) beam vector
         *
         */

        /* We will use beam vector and hades-p1 to find the assumed reaction origin.
         * But if the distance between tracks is to large, we want to drop this track.
         * We will count only vectors which passed the selection.
         */

        // STEP 1
        HGeomVector p1_base, p1_dir;
        HParticleTool::calcSegVector(
            p_cand->getZ(), p_cand->getR(), p_cand->getPhi() * TMath::DegToRad(), p_cand->getTheta() * TMath::DegToRad(), p1_base, p1_dir);

        HGeomVector p2_base, p2_dir;
        HParticleTool::calcSegVector(
            f_cand->getZ(), f_cand->getR(), f_cand->getPhi() * TMath::DegToRad(), f_cand->getTheta() * TMath::DegToRad(), p2_base, p2_dir);

        auto the_vertex = HParticleTool::calcVertexAnalytical(p1_base, p1_dir, p2_base, p2_dir);

        HGeomVector be_base(the_vertex);
        HGeomVector be_dir(0, 0, 1);

        if (use_beam_tilt) {
            be_dir = HGeomVector(beam->X() / beam->Z(), beam->Y() / beam->X(), 1.0);
        }

        // STEP 2+3
        auto dist_beam_p1 = HParticleTool::calculateMinimumDistance(p1_base, p1_dir, be_base, be_dir);
        auto dist_p1_p2 = HParticleTool::calculateMinimumDistance(p1_base, p1_dir, p2_base, p2_dir);

        // STEP 5
        TVector3 the_track;
        HGeomVector the_base;

        if (!project) {
            the_track = TVector3(p2_dir.X(), p2_dir.Y(), p2_dir.Z());
            the_track *= 1. / the_track.Z();
            the_base = HGeomVector(p2_base.X() - the_track.X() * p2_base.Z(), p2_base.Y() - the_track.Y() * p2_base.Z(), 0.);
        } else {
            // Define reference beam vector
            TLorentzVector beam_lvec;
            beam_lvec.SetPxPyPzE(0, 0, beam_p01, beam_E01);

            // Find expected elastics vector
            auto miss_lvec = (beam_tilt ? *beam : beam_lvec) - *p_cand;

            // Extract direction and normalize to Z() == 1.0
            the_track = miss_lvec.Vect();
            the_track *= 1. / the_track.Z();

            // Calculate BASE for F-track, project it to Z=0
            the_base = HGeomVector(the_vertex.X() - the_track.X() * the_vertex.Z(), the_vertex.Y() - the_track.Y() * the_vertex.Z(), 0.);
        }

        // mille.add_measurement(0,
        //                       1.0,
        //                       hsa::XYZPoint(the_base.X() - the_track.X() * the_base.Z(), the_base.Y() - the_track.Y() * the_base.Z(), 0),
        //                       hsa::XYZVector(the_track.X() / the_track.Z(), the_track.Y() / the_track.Z(), 1),
        //                       0,
        //                       0,
        //                       0,
        //                       dist_p1_p2,
        //                       1.0);

        auto proj_phi_diff = f_cand->getPhi() - the_track.Phi() * TMath::RadToDeg();
        if (proj_phi_diff > 180)
            proj_phi_diff += -180;
        auto proj_theta_diff = f_cand->getTheta() - the_track.Theta() * TMath::RadToDeg();
        if (proj_theta_diff > 180)
            proj_theta_diff += -180;

        // printf("ANG DIFF: %f %f\n", proj_phi_diff, proj_theta_diff);
        if (qa_file.length()) {
            h_vert_reco_xy->Fill(be_base.X(), be_base.Y());
            h_vert_reco_z->Fill(be_base.Z());

            h_fcand_chi2_ndf->Fill(f_cand->getChi2() / f_cand->getNDF());

            h_p1_beam_dist->Fill(dist_beam_p1);
            h_p1_p2_dist->Fill(dist_p1_p2);

            h_delta_phi->Fill(phi_diff);
            h_tan_prod->Fill(theta_tan_prod);
            h_el_theta_phi->Fill(phi_diff, theta_tan_prod);

            h_fcand_txty->Fill(p2_dir.X() / p2_dir.Z(), p2_dir.Y() / p2_dir.Z());
            h_fproj_txty->Fill(the_track.X() / the_track.Z(), the_track.Y() / the_track.Z());
            h_fproj_vert->Fill(the_base.X(), the_base.Y());
            h_fcand_theta_phi->Fill(f_cand->getPhi(), f_cand->getTheta());

            h_proj_diff_theta_phi->Fill(proj_phi_diff, proj_theta_diff);
        }
        // auto p2_p1_diff = p2_lvec - *f_cand;

        if (sim and f_cand->getChi2() / f_cand->getNDF() > 0.3) {
            rejected[6]++;
            continue;
        }
        if (!sim and f_cand->getChi2() / f_cand->getNDF() > 1.0) {
            rejected[6]++;
            continue;
        }

        if (sim and dist_beam_p1 > 0.5) {
            rejected[7]++;
            // continue;
        }

        if (!sim and dist_beam_p1 > 0.5) {
            rejected[7]++;
            // continue;
        }

        if (sim and dist_p1_p2 > 0.5) {
            rejected[8]++;
            // continue;
        }

        if (!sim and dist_p1_p2 > 0.5) {
            rejected[8]++;
            // continue;
        }

        // f_cand->print();
        // printf(" BEAM  "); beam_lvec.Print();
        if (verbose) {
            printf(" P1-H    ");
            (p_cand->Vect() * (1.0 / p_cand->Vect().Z())).Print();
            printf(" P1-H B  ");
            p1_base.print();
            printf(" P1-H T  ");
            p1_dir.print();
            printf(" P2-F  ");
            (f_cand->Vect() * (1.0 / f_cand->Vect().Z())).Print();
            printf(" P2-F B  ");
            p2_base.print();
            printf(" P2-F T  ");
            p2_dir.print();
            // printf(" P2    ");
            // (miss_lvec.Vect() * (1.0 / miss_lvec.Vect().Z())).Print();

            printf("  Add: %ld   BASE   %f  %f  %f   DIR   %f  %f  %f   ALG   %f  %f  %f   RES   %f   %f\n",
                   good_vectors_count,
                   the_base.X(),
                   the_base.Y(),
                   the_base.Z(),
                   the_track.X(),
                   the_track.Y(),
                   the_track.Z(),
                   0.0,
                   0.0,
                   0.0,
                   0.0,
                   0.1);
        }

        // straw_planes1.plane(0).add_measurement(
        //                       sigma,
        //                       hsa::XYZPoint(the_base.X() - the_track.X() * the_base.Z(), the_base.Y() - the_track.Y() * the_base.Z(), 0),
        //                       hsa::XYZVector(the_track.X() / the_track.Z(), the_track.Y() / the_track.Z(), 1),
        //                       hsa::XYZPoint(straw_loc.getX(), straw_loc.getY(), straw_loc.getZ()),
        //                       hsa::XYZVector(0, 1, 0),
        //                       dr);

        // straw_planes.plane(0).add_measurement(
        //     sigma,
        //     hsa::XYZPoint(the_base.X() - the_track.X() * the_base.Z(), the_base.Y() - the_track.Y() * the_base.Z(), 0),
        //     hsa::XYZVector(the_track.X() / the_track.Z(), the_track.Y() / the_track.Z(), 1),
        //     hsa::XYZPoint(straw_loc.getX(), straw_loc.getY(), straw_loc.getZ()),
        //     hsa::XYZVector(0, 1, 0),
        //     dr);

        auto track_base = hsa::XYZPoint(the_base.X() - the_track.X() * the_base.Z(), the_base.Y() - the_track.Y() * the_base.Z(), 0);
        auto track_dir = hsa::XYZVector(the_track.X() / the_track.Z(), the_track.Y() / the_track.Z(), 1);

        if (log_file.is_open()) {
            log_file << "===NEW_EVENT===\n";
            log_file << "-TRACK " << track_base.x() << ' ' << track_base.y() << ' ' << track_base.z() << ' ' << track_dir.x() << ' '
                     << track_dir.y() << ' ' << track_dir.z() << '\n';
        }

        const int Nev = f_cand->getNofHits();
        for (int k = 0; k < Nev; ++k) {
            HStsCal* fstscal = HCategoryManager::getObject(fstscal, fStsCal, f_cand->getStsHitIndex(k));
            Char_t mod, lay, ud;
            Int_t straw;

            fstscal->getAddress(mod, lay, straw, ud);

            const auto& straw_loc = *(stsCellsLoc[mod][lay][straw]);

            if (verbose) {
            }
            auto dr = f_cand->getStsDriftRadius(k);
            // printf(" Add: %d BASE %f %f %f DIR %f %f %f ALG %f %f %f RES %f %f\n",
            //        mod * 4 + lay,
            //        p2_base.X(),
            //        p2_base.Y(),
            //        p2_base.Z(),
            //        p2_dir.X(),
            //        p2_dir.Y(),
            //        p2_dir.Z(),
            //        u,
            //        0.,
            //        0.,
            //        dr,
            //        0.1);

            auto straw_base = hsa::XYZPoint(straw_loc.getX(), straw_loc.getY(), straw_loc.getZ());
            auto straw_dir = hsa::XYZVector(0, 1, 0);

            auto plane = straw_planes.plane((mod + 1) * 100 + (lay + 1) * 10)
                             .add_measurement(sigma, track_base, track_dir, straw_base, straw_dir, dr);

            auto model = plane.model();

            if (log_file.is_open()) {
                log_file << "-PLANE " << (int)mod << ' ' << (int)lay << ' ' << straw << ' ' << (int)ud << '\n';
                log_file << "-STRAW " << straw_base.x() << ' ' << straw_base.y() << ' ' << straw_base.z() << ' ' << straw_dir.x() << ' '
                         << straw_dir.y() << ' ' << straw_dir.z() << '\n';
                log_file << model;
            }

            if (qa_file.length()) {
                h_qa_sign->Fill(model.res_sign);
                h_qa_t_d_dist->Fill(model.track_wire_distance);
                h_qa_drift_dist->Fill(dr);
                h_qa_residuals->Fill(model.residual());
                h_qa_uresiduals->Fill(model.residual() * model.res_sign);

                auto plane_id = mod * 4 + lay + 1;
                h_qa_sign_plane->Fill(model.res_sign, plane_id);
                h_qa_t_d_dist_plane->Fill(model.track_wire_distance, plane_id);
                h_qa_drift_dist_plane->Fill(dr, plane_id);
                h_qa_residuals_plane->Fill(model.residual(), plane_id);
                h_qa_uresiduals_plane->Fill(model.residual() * model.res_sign, plane_id);
            }

            if (verbose) {
                plane.model().print_params();
            }
        }

        good_vectors_count++;
        pro_mille.end();
    }

    printf("Got %ld good vectors, rejected: ", good_vectors_count);
    for (int i = 0; i < rejected.size(); ++i)
        printf("  %ld", rejected[i]);
    printf("\n");

    TF1* f_res = new TF1("f_res", "gaus(0)+pol0(3)", -300, 300);
    f_res->SetParameters(100, 0, 10, 10);

    TText* tt = nullptr;
    if (sim)
        tt = new TText(0.5, 0.5, "Simulations");
    else
        tt = new TText(0.5, 0.5, "Data");
    tt->SetTextAlign(22);
    tt->SetTextColor(kRed + 2);
    tt->SetTextFont(43);
    tt->SetTextSize(20);
    tt->Draw();

    TFile* qafile {nullptr};

    if (qa_file.length()) {
        qafile = TFile::Open(qa_file.c_str(), "RECREATE");

        if (qafile) {
            auto can_vertex = new TCanvas("can_vertex", "can_vertex", 1200, 900);
            can_vertex->Divide(4, 3);

            can_vertex->cd(1);
            BeamX.Draw();
            BeamX.Write();

            can_vertex->cd(2);
            BeamY.Draw();
            BeamY.Write();

            can_vertex->cd(3);
            BeamXY.Draw("colz");
            BeamXY.Write();

            can_vertex->cd(4);
            BeamZ.Draw();
            BeamZ.Write();

            can_vertex->cd(5);
            h_vert_reco_xy->Draw("colz");
            h_vert_reco_xy->Write();

            can_vertex->cd(6);
            h_vert_reco_z->Draw();
            h_vert_reco_z->Write();

            can_vertex->cd(7);
            h_p1_p2_dist->Draw();
            h_p1_p2_dist->Write();

            can_vertex->cd(8);
            h_p1_beam_dist->Draw();
            h_p1_beam_dist->Write();

            can_vertex->cd(9);
            h_fcand_chi2_ndf->Draw();
            h_fcand_chi2_ndf->Write();

            can_vertex->Write();

            auto can_elastics = new TCanvas("can_elastics", "can_elastics", 1200, 900);
            can_elastics->Divide(4, 3);
            can_elastics->cd(1);
            h_delta_phi->Draw();
            h_delta_phi->Write();

            can_elastics->cd(2);
            h_tan_prod->Draw();
            h_tan_prod->Write();

            // can_elastics->cd(3);
            // h_fcand_delta->Draw("colz");
            // h_vert_reco_xy->Write();

            can_elastics->cd(3);
            h_elastic_error->Draw("text,h");
            h_elastic_error->Write();
            gPad->SetLogy();

            can_elastics->cd(4);
            h_elastic_mult->Draw("colz");
            h_elastic_mult->Write();

            can_elastics->cd(5);
            h_fcand_txty->Draw("colz");
            h_fcand_txty->Write();

            can_elastics->cd(6);
            h_fproj_txty->Draw("colz");
            h_fproj_txty->Write();

            can_elastics->cd(7);
            h_fproj_vert->Draw("colz");
            h_fproj_vert->Write();

            can_elastics->cd(8);
            h_fcand_theta_phi->Draw("colz");
            h_fcand_theta_phi->Write();

            can_elastics->cd(9);
            // h_el_theta_phi_all->Draw("colz");
            h_el_theta_phi->Draw("colz");
            h_el_theta_phi->Write();
            // h_el_theta_phi_rejected->Draw("colz,same");

            can_elastics->cd(12);
            h_proj_diff_theta_phi->Draw("colz");
            h_proj_diff_theta_phi->Write();

            can_elastics->Write();

            auto can_qa = new TCanvas("can_qa", "can_qa", 1500, 600);
            can_qa->Divide(5, 2);

            can_qa->cd(1);
            h_qa_sign->Draw("h,text");
            h_qa_sign->Write();

            can_qa->cd(2);
            h_qa_t_d_dist->Draw();
            h_qa_t_d_dist->Write();

            can_qa->cd(3);
            h_qa_drift_dist->Draw();
            h_qa_drift_dist->Write();

            can_qa->cd(4);
            h_qa_residuals->Draw();
            h_qa_residuals->Write();

            can_qa->cd(5);
            h_qa_uresiduals->Draw();
            h_qa_uresiduals->Write();

            can_qa->cd(6);
            h_qa_sign_plane->Draw("colz,text");
            h_qa_sign_plane->Write();

            can_qa->cd(7);
            h_qa_t_d_dist_plane->Draw("colz");
            h_qa_t_d_dist_plane->Write();

            can_qa->cd(8);
            h_qa_drift_dist_plane->Draw("colz");
            h_qa_drift_dist_plane->Write();

            can_qa->cd(9);
            h_qa_residuals_plane->Draw("colz");
            h_qa_residuals_plane->Write();

            can_qa->cd(10);
            h_qa_uresiduals_plane->Draw("colz");
            h_qa_uresiduals_plane->Write();

            can_qa->Write();

            qafile->Close();
        }
    }
}

auto library::find_beam_avgs(long long nevents) -> HGeomVector
{
    decltype(nevents) cnt = 0;
    for (int event = 0l; event < nevents; ++event) {
        loop->nextEvent(event);

        HEventHeader* header = gHades->getCurrentEvent()->getHeader();
        HVertex vertex = header->getVertexReco();
        if (vertex.getChi2() < 0)
            continue;

        cnt++;
        BeamX.Fill(vertex.getX());
        BeamY.Fill(vertex.getY());
        BeamZ.Fill(vertex.getZ());
    }

    BeamX.Print();
    BeamY.Print();
    BeamZ.Print();

    HGeomVector res;

    TF1 f_avg("f_avg", "gaus(0)", -5, 5);
    BeamX.Fit(&f_avg);
    res.X() = f_avg.GetParameter(1);
    BeamY.Fit(&f_avg);
    res.Y() = f_avg.GetParameter(1);
    BeamZ.Fit(&f_avg);
    res.Z() = f_avg.GetParameter(1);

    printf("== BEAM AVG: Events to analyze: %lld   Analyzed: %lld\n", nevents, cnt);

    return res;
}
