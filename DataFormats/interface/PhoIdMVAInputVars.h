#ifndef FLASHgg_PhoIdMVAInputVars_h
#define FLASHgg_PhoIdMVAInputVars_h

#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"

namespace flashgg {

    /** class to keep the actual values of the input variables to the photon MVA id.
        Similar to class DiPhotonMVAResult but without the MVA output value. */
    class PhoIdMVAInputVars
    {

    public:
        PhoIdMVAInputVars();
        virtual ~PhoIdMVAInputVars() {}

        // Input variables
        float scRawE;
        float r9;
        float covIEtaIEta;
        float phiWidth;
        float etaWidth;
        float covIEtaIPhi;
        float s4;
        float pfPhoIso03;
        float pfPhoIso03Corr;
        float pfChgIso03;
        float pfChgIso03worst;
        float scEta;
        float rho;
        float esEffSigmaRR;
        float esEnovSCRawEn;
    };

}

#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

