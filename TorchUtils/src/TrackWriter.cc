#include "flashgg/TorchUtils/interface/TrackWriter.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "DataFormats/Common/interface/RefToPtr.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "FastSimulation/Particle/interface/RawParticle.h"
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"


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

  void TrackWriter::addPhoton(const edm::Ptr<pat::Photon> &photon, const edm::Ptr<reco::Vertex> &photonVertex,
			      const flashgg::VertexCandidateMap &vtxcandmap
			      )
  {
    edm::Handle<edm::View<pat::PackedCandidate> > patCandidates;
    this->event->getByToken(packedCandidatesToken, patCandidates);

    // find tracks in a cone of dR (e.g. = 0.3) around this photon
    // vector<pat::PackedCandidate>          "packedPFCandidates"        ""                "PAT"
    
    double photonEta = photon->eta(); 
    double photonPhi = photon->phi();
    double photonEt  = photon->et();

    vector<float> trackpt, etaAtVertex, phiAtVertex;
    vector<int> charge, pdgId;

    /** coordinates of vertex associated to track  */
    vector<float> vtxX, vtxY, vtxZ;
    vector<unsigned> vtxIndex;

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
	if (vetoPackedCand(*photon, candPtr))
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

	  // note that we take the coordinates of the vertex to which the
	  // track is associated, NOT the vertex returned by the track 
	  vtxX.push_back(vtx->x());
	  vtxY.push_back(vtx->y());
	  vtxZ.push_back(vtx->z());

	  vtxIndex.push_back(vtx.key());

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
    this->vtxIndex.push_back(vtxIndex);
  }

  //----------------------------------------------------------------------

  void TrackWriter::propagateToECAL(const pat::PackedCandidate &cand, 
				    const reco::Vertex &vtx,
				    double &etaAtEcal, double &phiAtEcal) {

    // code adapted from RecoParticleFlow/PFTracking/src/PFTrackTransformer.cc
    BaseParticlePropagator theParticle =
      BaseParticlePropagator(
			     RawParticle(XYZTLorentzVector(cand.px(),
							   cand.py(),
							   cand.pz(),
							   cand.energy()),
					 XYZTLorentzVector(vtx.x(),
							   vtx.y(),
							   vtx.z(),
							   0.)),
			     0.,0.,magneticField.z());

    theParticle.setCharge(cand.charge());

    // propagate to ECAL
    theParticle.propagateToEcalEntrance(false);

    // get the new momentum coordinates

    if(theParticle.getSuccess() == 0) {
      // propagation failed for some reason (e.g. did not intersect with ECAL barrel)
      etaAtEcal = -9999;
      phiAtEcal = -9999;
    } else {
      etaAtEcal = theParticle.eta();
      phiAtEcal = theParticle.phi();
    }
  }

  //----------------------------------------------------------------------

  TrackWriter::~TrackWriter() 
  {
  }

  //----------------------------------------------------------------------

  void 
  TrackWriter::beginRun(const edm::Run& run, const edm::EventSetup& iSetup) {
     
    // from RecoParticleFlow/PFTracking/plugins/PFV0Producer.cc
    edm::ESHandle<MagneticField> magneticField;
    iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  
    this->magneticField = math::XYZVector(magneticField->inTesla(GlobalPoint(0,0,0)));
  }


  //----------------------------------------------------------------------

}

