#ifndef HFORWARDALIGNER_H
#define HFORWARDALIGNER_H

#include "hforwardtools.h"
#include "hgeomvector.h"
#include "hreconstructor.h"
#include "stsdef.h"

class HCategory;
class HForwardCand;
class HForwardCandFinderPar;
class HFRpcHit;
class HStsCalPar;
class HStsGeomPar;

namespace promille {
    template<typename T>
    class mille;

    template<typename T, typename ResidualModel>
    struct model_planes;

    template<typename T, size_t Nlocals, size_t Nglobals>
    struct residual_model_base;

    namespace euler
    {
        template<typename T>
        struct euler_base;

        template<typename T>
        struct zyz;
    }
}

namespace hsa
{
    template<typename T, template<class> class R, size_t Nlocals, size_t Nglobals>
    struct sts_residual;
}

class HForwardAligner : public HReconstructor
{
public:
    HForwardAligner(const Text_t* name, const Text_t* title, const Text_t* output_bin_file);
    virtual ~HForwardAligner();

    virtual Bool_t init();
    virtual Bool_t reinit();
    virtual Int_t execute();

protected:

private:
    HCategory* pKine;                   // Kine category for sim
    HCategory* fStsCal;                 // Input array of straw hits
    HCategory* pFRpcHits;               // Input array of rpc hits
    HCategory* pForwardCand;            // Output array of candidates
    HStsGeomPar* pStsGeomPar;           // strips geometry
    HStsCalPar* pCalPar;                // cal run par
    HForwardCandFinderPar* pStrawVFPar; // vector finder parameters
    Bool_t isSimulation;                // flag to mark simulation run
    HGeomVector primary_vertex;         // the primary vector cluster

    HGeomVector* stsCellsLoc[STS_MAX_MODULES][STS_MAX_LAYERS][STS_MAX_CELLS]; // Centre of the strip [mm]
    HGeomVector* stsCellsLab[STS_MAX_MODULES][STS_MAX_LAYERS][STS_MAX_CELLS]; // Centre of the strip [mm]

    Int_t nLayers[STS_MAX_MODULES];
    Float_t cosa[STS_MAX_MODULES][STS_MAX_LAYERS];
    Float_t sina[STS_MAX_MODULES][STS_MAX_LAYERS];
    Float_t roty[STS_MAX_MODULES][STS_MAX_LAYERS];

    promille::mille<float> * pro_mille;
    promille::model_planes<float, hsa::sts_residual<float, promille::euler::zyz, 4, 12>> * straw_planes;

};

#endif

