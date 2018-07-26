#ifndef TorchUtils_TrackWriter_h
#define TorchUtils_TrackWriter_h

#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "flashgg/DataFormats/interface/Photon.h"


#include "flashgg/TorchUtils/interface/TorchWriter.h"
#include "flashgg/DataFormats/interface/VertexCandidateMap.h"


namespace flashgg
{

  /** class to find and write out track information 
      of tracks close to photons */
  class TrackWriter
  {
  protected:

    edm::EDGetTokenT<edm::View<pat::PackedCandidate> > packedCandidatesToken;

    /** maximum delta R for a track to be added */
    const double maxDeltaR = 0.4;

    /** keep this during one event, assuming we do not run in multi-threaded
	mode */
    edm::Handle<edm::View<pat::PackedCandidate> > patCandidates;

    const edm::Event *event;

    math::XYZVector magneticField;

    /** things to write out.

	First index is the photon index,
	second index is the index within the given photon.
    */
    

    /** pt divided by photon et */
    std::vector<std::vector<float> > trackpt;
    std::vector<std::vector<float> > etaAtVertex;
    std::vector<std::vector<float> > phiAtVertex;
    std::vector<std::vector<int> >   charge, pdgId;

    /** coordinates of vertex associated to track */
    std::vector<std::vector<float> > vtxX, vtxY, vtxZ;
    std::vector<std::vector<unsigned> > vtxIndex;

    //----------------------------------------

  public:

    /** see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEDMGetDataFromEvent#Consumes_and_Helpers */
    TrackWriter(edm::ParameterSet const& iPS, edm::ConsumesCollector && iC);

    void beginRun(const edm::Run& run, const edm::EventSetup& iSetup);

    /** called at the beginning of each event */
    void newEvent(const edm::Event &event);
    
    /** finds tracks close to this photon and stores the corresponding data */
    void addPhoton(const edm::Ptr<pat::Photon> &photon, const edm::Ptr<reco::Vertex> &photonVertex,
		   const flashgg::VertexCandidateMap &vtxcandmap);

    void propagateToECAL(const pat::PackedCandidate &cand, 
			 const reco::Vertex &vtx,
			 double &etaAtEcal, double &phiAtEcal);

    virtual ~TrackWriter();

    static bool vetoPackedCand( const pat::Photon &photon, const edm::Ptr<pat::PackedCandidate> &pfcand );

  };
}

#endif
