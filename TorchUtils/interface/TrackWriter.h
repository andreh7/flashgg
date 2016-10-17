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

    /** things to write out.

	First index is the photon index,
	second index is the index within the given photon.
    */
    

    /** pt divided by photon et */
    std::vector<std::vector<float> > relpt;
    std::vector<std::vector<float> > detaAtVertex;
    std::vector<std::vector<float> > dphiAtVertex;
    std::vector<std::vector<int> >   charge;

    /** dz(track vertex - photon vertex) */
    std::vector<std::vector<float> > vtxDz;

    template<typename DataType>
    void writeFlattenedVector(TorchWriter &tw, 
			      const std::vector<std::vector<DataType> > &values,
			      unsigned totNumItems
			      );

  public:

    /** see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEDMGetDataFromEvent#Consumes_and_Helpers */
    TrackWriter(edm::ParameterSet const& iPS, edm::ConsumesCollector && iC);

    /** called at the beginning of each event */
    void newEvent(const edm::Event &event);
    
    /** finds tracks close to this photon and stores the corresponding data */
    void addPhoton(const flashgg::Photon &photon, const edm::Ptr<reco::Vertex> &photonVertex);

    /** called at the end to write out the collected data */
    void writeOut(TorchWriter &tw);

  };
}

#endif
