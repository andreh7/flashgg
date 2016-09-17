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
  void TrackWriter::addPhoton(const flashgg::Photon &photon)
  {
    edm::Handle<edm::View<pat::PackedCandidate> > patCandidates;
    this->event->getByToken(packedCandidatesToken, patCandidates);

    // find tracks in a cone of dR (e.g. = 0.3) around this photon
    // vector<pat::PackedCandidate>          "packedPFCandidates"        ""                "PAT"
    
    double photonEta = photon.eta(); 
    double photonPhi = photon.phi();
    double photonEt  = photon.et();

    vector<float> relpt, detaAtVertex, dphiAtVertex;
    vector<int> charge;

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

	if (dr2 < maxDeltaR * maxDeltaR)
        {
	  // keep this track for this photon
	  relpt.push_back(track->pt() / photonEt);
	  detaAtVertex.push_back(deta);
	  dphiAtVertex.push_back(dphi);
	  charge.push_back(track->charge());
        }
      } // loop over tracks of this event

    // add this photon
    this->relpt.push_back(relpt);
    this->detaAtVertex.push_back(detaAtVertex);
    this->dphiAtVertex.push_back(dphiAtVertex);
    this->charge.push_back(charge);

  }

  //----------------------------------------------------------------------

  template<typename DataType>
  void TrackWriter::writeFlattenedVector(TorchWriter &tw, 
					 const std::vector<std::vector<DataType> > &values,
					 unsigned totNumItems
					 )
  {
    vector<DataType> flatValues(totNumItems);

    unsigned nextPos = 0;

    for (unsigned i = 0; i < values.size(); ++i)
    {
      for (unsigned j = 0; j < values[i].size(); ++j)
	flatValues[nextPos++] = values[i][j];
  
    } // loop over tracks

    assert(nextPos == totNumItems);

    tw.writeTypeVector<DataType>(flatValues);
  }

  //----------------------------------------------------------------------

  void TrackWriter::writeOut(TorchWriter &tw)
  {
    // fill tensors to be written out, similar to photon rechits sparse
    // format

    //----------
    // calculate first indices
    //----------
    vector<int32_t> firstIndex(relpt.size()), numTracks(relpt.size());

    // Torch's indexing is one based
    int32_t nextStartIndex = 1;

    for (unsigned i = 0; i < relpt.size(); ++i)
    {
      firstIndex[i] = nextStartIndex;
      numTracks[i] = relpt[i].size();
      
      nextStartIndex += relpt[i].size();
    } // loop over photons

    //----------
    
    int32_t totNumTracks = nextStartIndex - 1;

    // number of tensors to write out
    const unsigned tableSize = 6;

    tw.writeInt(tw.MAGIC_TABLE);
    tw.writeInt(tw.getNextObjectIndex());
    tw.writeInt(tableSize); // number of variables

    tw.writeInt(tw.MAGIC_STRING); tw.writeString("firstIndex");  tw.writeTypeVector(firstIndex);
    tw.writeInt(tw.MAGIC_STRING); tw.writeString("numTracks");  tw.writeTypeVector(numTracks);


    tw.writeInt(tw.MAGIC_STRING); tw.writeString("relpt"       ); writeFlattenedVector(tw,relpt, totNumTracks);
    tw.writeInt(tw.MAGIC_STRING); tw.writeString("detaAtVertex"); writeFlattenedVector(tw,detaAtVertex, totNumTracks);
    tw.writeInt(tw.MAGIC_STRING); tw.writeString("dphiAtVertex"); writeFlattenedVector(tw,dphiAtVertex, totNumTracks);
    tw.writeInt(tw.MAGIC_STRING); tw.writeString("charge"      ); writeFlattenedVector(tw,charge, totNumTracks);
  }

  //----------------------------------------------------------------------

}

