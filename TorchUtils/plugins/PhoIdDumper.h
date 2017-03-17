#ifndef TorchUtils_PhoIdDumper_h
#define TorchUtils_PhoIdDumper_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"


#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/TorchUtils/interface/PhoIdWriter.h"
#include "flashgg/DataFormats/interface/DiPhotonPhoIdMVAInputVars.h"
#include "flashgg/TorchUtils/interface/TrackWriter.h"

namespace flashgg {
    class PhoIdDumper : public edm::EDAnalyzer
    {
    public:
        explicit PhoIdDumper( const edm::ParameterSet & );
        virtual ~PhoIdDumper();

        static void fillDescriptions( edm::ConfigurationDescriptions &descriptions );


    protected:

        virtual void beginJob() override;
        virtual void analyze( const edm::Event &, const edm::EventSetup & ) override;
        virtual void endJob() override;

        edm::EDGetTokenT<edm::View<flashgg::DiPhotonCandidate> > diphotonToken_;

        /** output file names for the Torch tensors with rechits */
        const std::string outputFname;
        
        /** whether to write the rechits in sparse or dense format */
        const bool writeRecHitsSparseFlag;

        /** window sizes */
        const unsigned windowHalfWidth;
        const unsigned windowHalfHeight;
        
        /** run/ls/event number for comparison purposes */
        std::vector<unsigned> runNumber, lsNumber;
        std::vector<unsigned long long> eventNumber;

        /** rechit information to be written out at the end */
        std::vector<std::vector<PhoIdWriter::RecHitData> > recHits;
        std::vector<float> weights;
        std::vector<float> labels; // 1 = prompt photon, 0 = fake or non-prompt photon
        std::vector<float> mvaID; // official MVA id
        std::vector<float> genDeltaR; // deltaR to matched gen photon

        std::vector<float> chgIsoWrtChosenVtx;
        std::vector<float> chgIsoWrtWorstVtx;

        /** whether to also write the official photon ID mva input variables or not */
        const bool writePhotonIdInputVarsFlag;

        /** value map for diphoton input variables association */
        edm::EDGetTokenT<flashgg::DiPhotonPhoIdMVAInputVarsAssociation> phoIdInputVarsToken;

        //----------
        // photon BDT id input variables
        //----------

        std::vector<float> scRawE,
            r9,
            covIEtaIEta,
            phiWidth,
            etaWidth,
            covIEtaIPhi,
            s4,
            pfPhoIso03,
            pfChgIso03,       // duplicate of above
            pfChgIso03worst,  // duplicate of above
            scEta,
            rho,
            esEffSigmaRR;

        //----------
        // other photon variables
        //----------
        std::vector<float> photonEt;

        //----------

        /** boundary between barrel and endcap */
        const float etaMaxBarrel = 1.5;

        TrackWriter *trackWriter;

        /** to write out the photon id variables */
        PhoIdWriter *phoIdWriter;
            
        
        //----------------------------------------

        virtual void wrapCoordinates(PhoIdWriter::RecHitData &rechit) = 0;

        /** finds the crystal with the maximum energy, normalizes to that and applies the given window around
            the maximum */
        void applyWindowAndNormalizeEnergy(std::vector<PhoIdWriter::RecHitData> &rechits, int windowHalfWidth, int windowHalfHeight);

        //----------------------------------------

        /** subdet specific: extract rechit coordinates */
        virtual void fillRecHits(const flashgg::Photon &photon, std::vector<PhoIdWriter::RecHitData> &rechits) = 0;

        /** called to check if a photon is in barrel/endcap */
        virtual bool isPhotonInSubdet(const flashgg::Photon &photon) = 0;

        void addPhoton(const edm::EventID &eventId, 
                       const flashgg::Photon &photon, 
                       const edm::Ptr<reco::Vertex> &photonVertex,
                       float weight, float mvaID,
                       float chosenVertexChargedIso,
                       float worstVertexChargedIso,
                       const PhoIdMVAInputVars *phoIdInputVars
                       );

        friend class PhoIdWriterTorch;
        friend class PhoIdWriterNumpy;
    };

    //----------------------------------------------------------------------

    /** dumps endcap photons */
    class PhoIdDumperBarrel : public PhoIdDumper
    {
        //----------------------------------------
    public:
        explicit PhoIdDumperBarrel( const edm::ParameterSet &iParams);

    protected:
        // wrap around in phi for the barrel
        virtual void wrapCoordinates(PhoIdWriter::RecHitData &rechit);

        virtual void fillRecHits(const flashgg::Photon &photon, std::vector<PhoIdWriter::RecHitData> &rechits);

        virtual bool isPhotonInSubdet(const flashgg::Photon &photon);
    };

    //----------------------------------------------------------------------

    /** dumps barrel photons */
    class PhoIdDumperEndcap : public PhoIdDumper
    {
    public:

        explicit PhoIdDumperEndcap( const edm::ParameterSet &iParams);

    protected:
        virtual void wrapCoordinates(PhoIdWriter::RecHitData &rechit);

        void fillRecHits(const flashgg::Photon &photon, std::vector<PhoIdWriter::RecHitData> &rechits);

        virtual bool isPhotonInSubdet(const flashgg::Photon &photon);
    };


} // namespace flashgg

#endif

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
