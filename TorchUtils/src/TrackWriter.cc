#include "flashgg/TorchUtils/interface/TrackWriter.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "DataFormats/Common/interface/RefToPtr.h"


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

  /** copied from class PhotonIdUtils */
  bool TrackWriter::vetoPackedCand( const pat::Photon &photon, const edm::Ptr<pat::PackedCandidate> &pfcand ) {
        edm::RefVector<pat::PackedCandidateCollection> associated =  photon.associatedPackedPFCandidates();

        for( unsigned int ipc = 0; ipc < associated.size(); ipc++ ) {
            edm::Ptr<pat::PackedCandidate> associatedPtr = edm::refToPtr( associated[ipc] );
            if( associatedPtr == pfcand )  { return true; }
        }

        return false;
    }

  //----------------------------------------

  void TrackWriter::addPhoton(const flashgg::Photon &photon, const edm::Ptr<reco::Vertex> &photonVertex,
			      const flashgg::VertexCandidateMap &vtxcandmap
			      )
  {
    edm::Handle<edm::View<pat::PackedCandidate> > patCandidates;
    this->event->getByToken(packedCandidatesToken, patCandidates);

    // find tracks in a cone of dR (e.g. = 0.3) around this photon
    // vector<pat::PackedCandidate>          "packedPFCandidates"        ""                "PAT"
    
    double photonEta = photon.eta(); 
    double photonPhi = photon.phi();
    double photonEt  = photon.et();

    vector<float> trackpt, etaAtVertex, phiAtVertex;
    vector<int> charge, pdgId;

    /** coordinates of vertex associated to track  */
    vector<float> vtxX, vtxY, vtxZ;

    // see https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/DataFormats/PatCandidates/interface/PackedCandidate.h
    // for PackedCandidate


    for (auto vtxToTrackIter = vtxcandmap.begin();
	 vtxToTrackIter != vtxcandmap.end();
	 ++vtxToTrackIter)
      {
	// note that here we do NOT require the vertex of 
	// the track to be the same as the (selected) one
	// of the photon since later we may also want
	// to calculate the charged isolation with respect
	// to other vertices
	const edm::Ptr<reco::Vertex> vtx = vtxToTrackIter->first;
	const edm::Ptr<pat::PackedCandidate> candPtr = vtxToTrackIter->second;

	const pat::PackedCandidate &cand = *candPtr;

	// skip candidate if it is part of the photon
	if (vetoPackedCand(photon, candPtr))
	  continue;

	// deltaPhi(A,B) calculates phi(A) - phi(B) (modulo wrapping around)
	// our 'reference' is the photon
	double dphi = deltaPhi(cand.phi(), photonPhi);
	double deta = cand.eta() - photonEta;
	double dr2 = dphi * dphi + deta * deta; 

	// photonVertex.id() and cand.vertexRef().id() seem to be the same all the time but 
	// the vertex locations are not -> decide by dz
	// (we actually store the dz into the output file
	// so we can decide later on whether this is the same vertex or not)

	if (dr2 < maxDeltaR * maxDeltaR)
        {
	  // keep this track for this photon
	  trackpt.push_back(cand.pt());
	  etaAtVertex.push_back(cand.eta());
	  phiAtVertex.push_back(cand.phi());
	  charge.push_back(cand.charge());
	  pdgId.push_back(cand.pdgId());

	  // track minus photon vertex
	  vtxX.push_back(cand.vertexRef()->x());
	  vtxY.push_back(cand.vertexRef()->y());
	  vtxZ.push_back(cand.vertexRef()->z());
        }

      } // loop over tracks of this event

    // add this photon
    this->trackpt.push_back(trackpt);
    this->etaAtVertex.push_back(etaAtVertex);
    this->phiAtVertex.push_back(phiAtVertex);
    this->charge.push_back(charge);
    this->pdgId.push_back(pdgId);
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

