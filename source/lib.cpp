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

// look & feel libs
#include <fmt/core.h>
#include <tabulate/table.hpp>

// System
#include <iostream>
#include <vector>

using std::cout;
using std::endl;

int forward_aligner_library::verbose = 0;

forward_aligner_library::forward_aligner_library(HLoop* loop,
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

    if (!loop->setInput("-*,+HForwardCand,+HStsCal,+HFRpcHit,+HStart2Hit")) {  // reading file structure
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

    pro_mille.add_global_parameter(1, 0, "Tx_g");
    pro_mille.add_global_parameter(2, 0, "Ty_g");
    pro_mille.add_global_parameter(3, 0, "Tz_g");
    pro_mille.add_global_parameter(5, 0, "Ra_g");
    pro_mille.add_global_parameter(6, 0, "Rb_g");
    pro_mille.add_global_parameter(7, 0, "Rc_g");

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
                           pro_mille.add_global_parameter((m + 1) * 100 + (l + 1) * 10 + 1, 0, fmt::format("Tx_c_{:d}{:d}", m+1, l+1)),
                           pro_mille.add_global_parameter((m + 1) * 100 + (l + 1) * 10 + 2, 0, fmt::format("Ty_c_{:d}{:d}", m+1, l+1)),
                           pro_mille.add_global_parameter((m + 1) * 100 + (l + 1) * 10 + 3, 0, fmt::format("Tz_c_{:d}{:d}", m+1, l+1)),
                           pro_mille.add_global_parameter((m + 1) * 100 + (l + 1) * 10 + 5, 0, fmt::format("Ra_c_{:d}{:d}", m+1, l+1)),
                           pro_mille.add_global_parameter((m + 1) * 100 + (l + 1) * 10 + 6, 0, fmt::format("Rb_c_{:d}{:d}", m+1, l+1)),
                           pro_mille.add_global_parameter((m + 1) * 100 + (l + 1) * 10 + 7, 0, fmt::format("Rc_c_{:d}{:d}", m+1, l+1)))
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

    // Checking if ParticleCand and ForwardCand are open
    if (!fForwardCand)
        abort();

    if (log_file_name.length()) {
        log_file = std::ofstream(log_file_name, std::ios::out);
    }
}

auto forward_aligner_library::execute(long long events, long long first) -> void
{
    // Popup the GUI...
    gStyle->SetOptStat(1);

    TH1I* h_fcand_chi2_ndf {nullptr};
    TH2I* h_fcand_txty {nullptr};
    TH2I* h_fcand_theta_phi {nullptr};
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
        h_fcand_chi2_ndf = new TH1I("h_fcand_chi2_ndf", "h_fcand_chi2_ndf", 50, 0, 5);

        h_fcand_txty = new TH2I("h_fcand_txty", "h_fcand_txty;tx;ty", 200, -0.2, 0.2, 200, -0.2, 0.2);
        h_fcand_theta_phi = new TH2I("h_fcand_theta_phi", "h_fcand_theta_phi;phi [deg];theta [deg]", 360, 0, 360, 200, 0, 10);

        h_qa_sign = new TH1I("h_qa_sign", "h_qa_sign;sign;counts", 2, -1, 2);
        h_qa_t_d_dist = new TH1I("h_qa_t_d_dist", "h_qa_t_d_dist;track-wire distance [mm];counts", 140, 0, 7);
        h_qa_drift_dist = new TH1I("h_qa_drift_dist", "h_qa_drift_dist;track-wire distance [mm];counts", 140, 0, 7);
        h_qa_residuals = new TH1I("h_qa_residuals", "h_qa_residuals;residuals [mm];counts", 100, -1, 1);
        h_qa_uresiduals = new TH1I("h_qa_uresiduals", "h_qa_uresiduals;residuals [mm];counts", 100, -1, 1);

        h_qa_sign_plane = new TH2I("h_qa_sign_plane", "h_qa_sign_plane;sign;plane;counts", 2, -1, 2, 8, 1, 9);
        h_qa_t_d_dist_plane =
            new TH2I("h_qa_t_d_dist_plane", "h_qa_t_d_dist_plane;track-wire distance [mm];plane;counts", 140, 0, 7, 8, 1, 9);
        h_qa_drift_dist_plane =
            new TH2I("h_qa_drift_dist_plane", "h_qa_drift_dist_plane;track-wire distance [mm];plane;counts", 140, 0, 7, 8, 1, 9);
        h_qa_residuals_plane = new TH2I("h_qa_residuals_plane", "h_qa_residuals_plane;residuals [mm];plane;counts", 200, -2, 2, 8, 1, 9);
        h_qa_uresiduals_plane = new TH2I("h_qa_uresiduals_plane", "h_qa_uresiduals_plane;residuals [mm];plane;counts", 200, -2, 2, 8, 1, 9);
    }

    auto nevts = events == 0 ? loop->getEntries() : min(first + events, loop->getEntries());

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

    for (auto event = first; event < nevts; ++event) {
        loop->nextEvent(event);
        if (event % 10000 == 0)
            fmt::print("{:d}/{:d}  {:.2f}%\n", event, nevts, 100.0 * event / nevts);

        auto forward_cand_cnt = fForwardCand->getEntries();

        if (forward_cand_cnt == 0) {
            // BADMULT
            rejected[1]++;
            continue;
        }

        for (auto i = 0u; i < forward_cand_cnt; ++i)
        {
            auto f_cand = dynamic_cast<HForwardCand*>(fForwardCand->getObject(i));

            if (!f_cand->isFlagBit(HForward::kIsUsed)) {
                // NOTUSED
                rejected[2]++;
                continue;
            }

            f_cand->calc4vectorProperties(HPhysicsConstants::mass(14));

            HGeomVector f_base, f_dir;
            HParticleTool::calcSegVector(
                f_cand->getZ(), f_cand->getR(), f_cand->getPhi() * TMath::DegToRad(), f_cand->getTheta() * TMath::DegToRad(), f_base, f_dir);

            if (qa_file.length()) {
                h_fcand_chi2_ndf->Fill(f_cand->getChi2() / f_cand->getNDF());

                h_fcand_txty->Fill(f_dir.X() / f_dir.Z(), f_dir.Y() / f_dir.Z());
                h_fcand_theta_phi->Fill(f_cand->getPhi(), f_cand->getTheta());
            }

            if (sim and f_cand->getChi2() / f_cand->getNDF() > 2.0) {
                rejected[6]++;
                continue;
            }
            if (!sim and f_cand->getChi2() / f_cand->getNDF() > 1.0) {
                rejected[6]++;
                continue;
            }

            TVector3 the_track = TVector3(f_dir.X(), f_dir.Y(), f_dir.Z());
            the_track *= 1. / the_track.Z();
            HGeomVector the_base = HGeomVector(f_base.X() - the_track.X() * f_base.Z(), f_base.Y() - the_track.Y() * f_base.Z(), 0.);

            if (verbose) {
                fmt::print("============================ EVENT: {} ============================\n", event);
                printf(" P2-F  ");
                (f_cand->Vect() * (1.0 / f_cand->Vect().Z())).Print();
                printf(" P2-F B  ");
                f_base.print();
                printf(" P2-F T  ");
                f_dir.print();
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
                       0.16);
            }

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

                const auto& straw_loc = *(stsCellsLoc[(int)mod][(int)lay][straw]);

                auto dr = f_cand->getStsDriftRadius(k);
                if (dr != dr) continue;


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
        }

        pro_mille.end();
    }

    fmt::print("End of evnt loop.\n");

    tabulate::Table rejected_stats;
    rejected_stats.add_row(tabulate::RowStream{} << "Good vectors" << "-- not used--" << "FCand cnt == 0"
    << "kIsUsed != 0" << "-- not used--" << "-- not used--"<< "-- not used--" << "Chi2 cut" );

    tabulate::RowStream rs{};
    rs << good_vectors_count;

    // fmt::print("Got {} good vectors, rejected: ", good_vectors_count);
    for (const auto & r : rejected)
        rs << r;
    // printf("  %ld", rejected[i]);

    // // fmt::print("Got {} good vectors, rejected: ", good_vectors_count);
    // for (int i = 0; i < rejected.size(); ++i)
    //     rs << rejected[i];
    //     // printf("  %ld", rejected[i]);
    printf("\n");

    // rejected_stats.add_row(rs);
    // std::cout << rejected_stats << '\n';

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
            auto can_vertex = new TCanvas("can_vertex", "can_vertex", 900, 600);
            can_vertex->Divide(3, 2);

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

            h_fcand_chi2_ndf->Draw();
            h_fcand_chi2_ndf->Write();

            can_vertex->Write();

            auto can_elastics = new TCanvas("can_elastics", "can_elastics", 900, 900);
            can_elastics->Divide(2, 2);
            can_elastics->cd(1);
            h_fcand_txty->Draw("colz");
            h_fcand_txty->Write();

            can_elastics->cd(2);
            h_fcand_theta_phi->Draw("colz");
            h_fcand_theta_phi->Write();

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

auto forward_aligner_library::find_beam_avgs(long long nevents) -> HGeomVector
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
        BeamXY.Fill(vertex.getX(), vertex.getY());
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

    printf("== BEAM AVG: Events to analyze: %lld   Analyzed: %lld\n", nevents, cnt);

    return res;
}
