#ifndef FLASHgg_ZHLeptonicTag_h
#define FLASHgg_ZHLeptonicTag_h

#include "flashgg/DataFormats/interface/DiPhotonTagBase.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"

namespace flashgg {

    class ZHLeptonicTag: public DiPhotonTagBase
    {
    public:
        ZHLeptonicTag();
        ZHLeptonicTag( edm::Ptr<DiPhotonCandidate>, edm::Ptr<DiPhotonMVAResult> );
        ZHLeptonicTag( edm::Ptr<DiPhotonCandidate>, DiPhotonMVAResult );
        ~ZHLeptonicTag();

        ZHLeptonicTag *clone() const override { return ( new ZHLeptonicTag( *this ) ); }
        
        const std::vector<edm::Ptr<Muon> > muons() const { return Muons_;}
        const std::vector<edm::Ptr<flashgg::Electron> > electrons() const {return Electrons_;}
        const std::vector<edm::Ptr<Jet> > jets() const { return Jets_;}
        const edm::Ptr<flashgg::Met>  met() const { return MET_;}
        
        void setJets( std::vector<edm::Ptr<Jet> > Jets ) { Jets_ = Jets; }
        void setMuons( std::vector<edm::Ptr<Muon> > Muons ) {Muons_ = Muons;}
        void setMET( edm::Ptr<flashgg::Met> MET ) {MET_ = MET;}
        void setElectrons( std::vector<edm::Ptr<Electron> > Electrons ) {Electrons_ = Electrons;}

        /** @return the transverse mass calculated from the highest pt lepton 
            and MET or -1 if there is no lepton */
        double transverseMass() const;

        /** checks if there are two same flavour opposite sign leptons and returns
            their mass or -1 if there isn't any */
        static double zcandMass(const std::vector<edm::Ptr<flashgg::Electron> > &electrons,
                                const std::vector<edm::Ptr<Muon> > &muons);


        double zcandMass() const;

        DiPhotonTagBase::tag_t tagEnum() const override {return DiPhotonTagBase::kZHLeptonic; }
    private:
        /** @param massDiffToZ will be set to the difference of the found lepton mass w.r.t to the 
            nominal Z mass */
        static double zcandMassHelper(const std::vector<const reco::Candidate*> &leptons, double &massDiffToZ);

    private:
        std::vector<edm::Ptr<Muon> > Muons_;
        std::vector<edm::Ptr<Electron> > Electrons_;
        std::vector<edm::Ptr<Jet> > Jets_;
        edm::Ptr<flashgg::Met> MET_;
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

