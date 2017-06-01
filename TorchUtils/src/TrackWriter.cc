#include "flashgg/TorchUtils/interface/TrackWriter.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include <iostream>

using namespace std;

namespace flashgg
{
  TrackWriter::TrackWriter(edm::ParameterSet const& iPS, edm::ConsumesCollector && iC)
  {
    packedCandidatesToken = iC.consumes<edm::View<pat::PackedCandidate> >(edm::InputTag("packedPFCandidates"));
  }

  //----------------------------------------------------------------------
  void TrackWriter::newEvent(const edm::Event &event)
  {
    // event.getByToken(packedCandidatesToken, patCandidates);
    this->event = &event;
  }
  //----------------------------------------------------------------------
  void TrackWriter::addPhoton(const flashgg::Photon &photon, const edm::Ptr<reco::Vertex> &photonVertex)
  {
    edm::Handle<edm::View<pat::PackedCandidate> > patCandidates;
    this->event->getByToken(packedCandidatesToken, patCandidates);

    // find tracks in a cone of dR (e.g. = 0.3) around this photon
    // vector<pat::PackedCandidate>          "packedPFCandidates"        ""                "PAT"
    
    double photonEta = photon.eta(); 
    double photonPhi = photon.phi();
    double photonEt  = photon.et();

    vector<float> trackpt, detaAtVertex, dphiAtVertex;
    vector<int> charge;

    /** coordinates of vertex associated to track  */
    vector<float> vtxX, vtxY, vtxZ;

    // see https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/DataFormats/PatCandidates/interface/PackedCandidate.h
    // for PackedCandidate
    for (const pat::PackedCandidate &cand : *patCandidates.product())
      {
	// get the track
	// see https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/DataFormats/TrackReco/interface/Track.h
	const reco::Track *track = cand.bestTrack();

	// cout << "track=" << track << endl;

	if (track == NULL)
	  continue;

	// deltaPhi(A,B) calculates phi(A) - phi(B) (modulo wrapping around)
	// our 'reference' is the photon
	double dphi = deltaPhi(track->phi(), photonPhi);
	double deta = track->eta() - photonEta;
	double dr2 = dphi * dphi + deta * deta; 

	// photonVertex.id() and cand.vertexRef().id() seem to be the same all the time but 
	// the vertex locations are not -> decide by dz
	// (we actually store the dz into the output file
	// so we can decide later on whether this is the same vertex or not)

	if (dr2 < maxDeltaR * maxDeltaR)
        {
	  // keep this track for this photon
	  trackpt.push_back(track->pt());
	  detaAtVertex.push_back(deta);
	  dphiAtVertex.push_back(dphi);
	  charge.push_back(track->charge());

	  // track minus photon vertex
	  vtxX.push_back(cand.vertexRef()->x());
	  vtxY.push_back(cand.vertexRef()->y());
	  vtxZ.push_back(cand.vertexRef()->z());
        }

      } // loop over tracks of this event

    // add this photon
    this->trackpt.push_back(trackpt);
    this->detaAtVertex.push_back(detaAtVertex);
    this->dphiAtVertex.push_back(dphiAtVertex);
    this->charge.push_back(charge);
    this->vtxX.push_back(vtxX);
    this->vtxY.push_back(vtxY);
    this->vtxZ.push_back(vtxZ);
  }

  //----------------------------------------------------------------------

  TrackWriter::~TrackWriter() 
  {
  }

  //----------------------------------------------------------------------

}

