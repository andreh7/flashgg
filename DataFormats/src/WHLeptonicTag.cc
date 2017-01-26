#include "flashgg/DataFormats/interface/WHLeptonicTag.h"
#include <algorithm>

#include "DataFormats/Math/interface/LorentzVector.h"


using namespace flashgg;

WHLeptonicTag::WHLeptonicTag() : DiPhotonTagBase::DiPhotonTagBase()
{}

WHLeptonicTag::~WHLeptonicTag()
{}


WHLeptonicTag::WHLeptonicTag( edm::Ptr<DiPhotonCandidate> diPho, edm::Ptr<DiPhotonMVAResult> mvares ) : DiPhotonTagBase::DiPhotonTagBase( diPho, *mvares ) {}
WHLeptonicTag::WHLeptonicTag( edm::Ptr<DiPhotonCandidate> dipho, DiPhotonMVAResult mvares ) : DiPhotonTagBase::DiPhotonTagBase( dipho, mvares ) {}

double WHLeptonicTag::transverseMass() const { return transverseMass(met(), electrons(), muons()); }

double 
WHLeptonicTag::transverseMass(const edm::Ptr<flashgg::Met> &met, 
                              const std::vector<edm::Ptr<flashgg::Electron> > &electrons,
                              const std::vector<edm::Ptr<Muon> > &muons) {

    // momentum of the highest pt particle
    const math::XYZTLorentzVector *lepton = NULL; 

    if (electrons.size() == 0) {

        if (muons.size() == 0)
            // no leptons at all
            return -1;
        else
            lepton = &(muons[0]->p4());

    } else {

        // we have electrons
        if (muons.size() == 0) {

            // we have only electrons, no muons
            lepton = &(electrons[0]->p4());

        } else {

            // must decide between electron and muon
            if (muons[0]->pt() > electrons[0]->pt()) {

                // take the muon
                lepton = &(muons[0]->p4());

            } else {
                // take the electron
                lepton = &(electrons[0]->p4());
            }
        }
    }

    const double leptonET = lepton->Et();
    const double leptonPhi = lepton->Phi();

    // MT^2 = 2 * ET1 * ET2 * (1 - cos (delta phi))

    const double metMag = met->getCorPt();
    const double metPhi = met->getCorPhi();

    const double dphi = deltaPhi(leptonPhi, metPhi);

    double mt2 = 2 * metMag * leptonET * ( 1 - cos(dphi));

    return sqrt(max(0., mt2));
}


// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

