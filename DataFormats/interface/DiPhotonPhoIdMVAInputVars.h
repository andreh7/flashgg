#ifndef FLASHgg_DiPhotonPhoIdMVAInputVars_h
#define FLASHgg_DiPhotonPhoIdMVAInputVars_h

#include "DataFormats/Common/interface/Ptr.h"
#include "flashgg/DataFormats/interface/PhoIdMVAInputVars.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"

namespace flashgg {

    /** class to keep the input values for the photon ID MVAs for a diphoton pair.
        We store these input variables as a pair because at least one input variable
        depends on the selected vertex which is taken from the diphoton candidate.
    */
    class DiPhotonPhoIdMVAInputVars
    {

    public:
        DiPhotonPhoIdMVAInputVars();
        virtual ~DiPhotonPhoIdMVAInputVars() {}

        PhoIdMVAInputVars &getInputsLeading() { return inputsLeading_; }
        PhoIdMVAInputVars &getInputsSubLeading() { return inputsSubLeading_; }

        const PhoIdMVAInputVars &getInputsLeading() const { return inputsLeading_; }
        const PhoIdMVAInputVars &getInputsSubLeading() const { return inputsSubLeading_; }
        
    private:
        /** the actual input variables */
        PhoIdMVAInputVars inputsLeading_, inputsSubLeading_;
   };

    /** maps from DiPhotonCandidate to DiPhotonPhoIdMVAInputVars. This can be used to get
        DiPhotonPhoIdMVAInputVars objects from DiPhotonCandidate objects without having
        to modify the structure of the latter. */
    typedef std::map<edm::Ptr<DiPhotonCandidate>, DiPhotonPhoIdMVAInputVars> DiPhotonPhoIdMVAInputVarsAssociation;
    // typedef edm::ValueMap<DiPhotonPhoIdMVAInputVars> DiPhotonPhoIdMVAInputVarsAssociation;
}

#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

