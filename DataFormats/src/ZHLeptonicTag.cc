#include "flashgg/DataFormats/interface/ZHLeptonicTag.h"
#include "flashgg/DataFormats/interface/WHLeptonicTag.h"
#include <algorithm>

using namespace flashgg;

ZHLeptonicTag::ZHLeptonicTag() : DiPhotonTagBase::DiPhotonTagBase()
{}

ZHLeptonicTag::~ZHLeptonicTag()
{}


ZHLeptonicTag::ZHLeptonicTag( edm::Ptr<DiPhotonCandidate> diPho, edm::Ptr<DiPhotonMVAResult> mvares ) : DiPhotonTagBase::DiPhotonTagBase( diPho, *mvares ) {}
ZHLeptonicTag::ZHLeptonicTag( edm::Ptr<DiPhotonCandidate> dipho, DiPhotonMVAResult mvares ) : DiPhotonTagBase::DiPhotonTagBase( dipho, mvares ) {}

double ZHLeptonicTag::transverseMass() const { return WHLeptonicTag::transverseMass(met(), electrons(), muons()); }

double ZHLeptonicTag::zcandMass(const std::vector<edm::Ptr<flashgg::Electron> > &electrons,
                                const std::vector<edm::Ptr<Muon> > &muons) {
    
    vector<const reco::Candidate*> leptons;

    // try electrons
    for (const edm::Ptr<flashgg::Electron> &electron : electrons) {
        leptons.push_back(&(*electron));
    }

    double eleMassDiff;
    double eleMass = zcandMassHelper(leptons, eleMassDiff);

    leptons.clear();

    // try muons electrons
    for (const edm::Ptr<Muon> &muon : muons) {
        leptons.push_back(&(*muon));
    }
    
    double muMassDiff;
    double muMass = zcandMassHelper(leptons, muMassDiff);

    // return the one closer to mZ = 91.2 GeV
    if (muMassDiff < eleMassDiff)
        return muMass;
    else
        return eleMass;

}

double 
ZHLeptonicTag::zcandMassHelper(const std::vector<const reco::Candidate*> &leptons, double &massDiffToZ) {

    double mass = -1;

    // assume leptons are all the same flavour, look 
    // for opposite sign pairs
    vector<const reco::Candidate *> posLeptons, negLeptons;
    
    for (const reco::Candidate* lepton : leptons) {
        if (lepton->charge() > 0)
            posLeptons.push_back(lepton);
        else if (lepton->charge() < 0)
            negLeptons.push_back(lepton);
    }

    // assume leptons are ordered in decreasing pt,
    // prefer higher pt leptons first
    if (posLeptons.size() > 0 && negLeptons.size() > 0) {
        mass = (posLeptons[0]->p4() + negLeptons[0]->p4()).M();
    }

    massDiffToZ = fabs(mass - 91.2);

    return mass;
}

double ZHLeptonicTag::zcandMass() const { return zcandMass(electrons(), muons()); }




// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

