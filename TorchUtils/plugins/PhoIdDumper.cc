#include "flashgg/TorchUtils/plugins/PhoIdDumper.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "flashgg/TorchUtils/interface/PhoIdWriterTorch.h"
#include "flashgg/TorchUtils/interface/TrackWriterTorch.h"
#include "flashgg/TorchUtils/interface/PhoIdWriterNumpy.h"
#include "flashgg/TorchUtils/interface/TrackWriterNumpy.h"

#include <boost/algorithm/string/predicate.hpp>


#include "flashgg/MicroAOD/interface/PhotonIdUtils.h"

#include "flashgg/TorchUtils/interface/TrackWriter.h"

#define MYDEBUG

using namespace std;
using namespace edm;

// **********************************************************************

namespace flashgg {

    //----------------------------------------------------------------------

    void PhoIdDumper::applyWindowAndNormalizeEnergy(std::vector<PhoIdWriter::RecHitData> &rechits, int windowHalfWidth, int windowHalfHeight, bool normalizeRecHitsToMax)
    {
        if (rechits.size() < 1)
            return;

        // find maximum energy rechit
        unsigned maxIndex = 0;

        for (unsigned i = 1; i < rechits.size(); ++i)
            if (rechits[i].energy > rechits[maxIndex].energy)
                maxIndex = i;

        float maxEnergy = rechits[maxIndex].energy;
        int xmax = rechits[maxIndex].dx;
        int ymax = rechits[maxIndex].dy;

        // coordinates are zero based
        const int centerX = windowHalfWidth;
        const int centerY = windowHalfHeight;

        // note the reverse order (for deleting rechits outside the window)
        for (int i = rechits.size() - 1; i >= 0; --i)
            {
                // normalize energy
                if (normalizeRecHitsToMax && maxEnergy > 0)
                    rechits[i].energy /= maxEnergy;
                
                // center coordinate
                rechits[i].dx -= xmax; rechits[i].dx += centerX;
                rechits[i].dy -= ymax; rechits[i].dy += centerY;

                // wrap around in phi for the barrel
                wrapCoordinates(rechits[i]);

                // apply window 
                if (rechits[i].dx < 0 || rechits[i].dx >= 2 * windowHalfWidth + 1 || 
                    rechits[i].dy < 0 || rechits[i].dy >= 2 * windowHalfHeight + 1)
                    // rechit is outside window
                    // note that we already copied the values of the max rechit
                    // so we don't have to worry about it shifting its index
                    rechits.erase(rechits.begin() + i);

            } // loop over rechits            
    }

    //----------------------------------------
    

    /** copy of PhotonIdUtils::pfIsoChgWrtVtx(..) but with slightly
        different parameter types */
    double pfIsoChgWrtVtx(const pat::Photon &photon,
                          const edm::Ptr<reco::Vertex> vtx,
                          const flashgg::VertexCandidateMap &vtxcandmap,
                          float coneSize, float coneVetoBarrel, float coneVetoEndcap,
                          float ptMin
                          )
    {
        bool removeOverlappingCandidates_ = true;
        OverlapRemovalAlgo *overlapAlgo_ = NULL;
        double deltaPhiRotation_ = 0;

        float isovalue = 0;

        float coneVeto = 0;
        if( photon.isEB() )      { coneVeto = coneVetoBarrel; }
        else if( photon.isEE() ) { coneVeto = coneVetoEndcap; }


        math::XYZVector SCdirection( photon.superCluster()->x() - vtx->x(),
                                     photon.superCluster()->y() - vtx->y(),
                                     photon.superCluster()->z() - vtx->z()
                                     );

        std::cout << "SELECTED VERTEX: " << vtx->z() 
                  << " PHOTON ET=" << photon.et()
                  << std::endl;
        // std::cout << "VERTICES:" << std::endl;
        // for (auto pp = vtxcandmap.begin(); pp != vtxcandmap.end(); ++pp) {
        //     std::cout << "  vtx: " << pp->first->z() << std::endl;
        // }

        auto mapRange = std::equal_range( vtxcandmap.begin(), vtxcandmap.end(), vtx, flashgg::compare_with_vtx() );
        if( mapRange.first == mapRange.second ) { 
            std::cout << "SKIPPING ALL CANDIDATES BECAUSE NO VERTICES" << std::endl;
            return -1.; } // no entries for this vertex
        for( auto pair_iter = mapRange.first ; pair_iter != mapRange.second ; pair_iter++ ) {
            edm::Ptr<pat::PackedCandidate> pfcand = pair_iter->second;

            if( abs( pfcand->pdgId() ) == 11 || abs( pfcand->pdgId() ) == 13 ) { 
                std::cout << "SKIPPING CANDIDATE WITH pt=" << pfcand->pt() << " BECAUSE OF PDGID" << std::endl;
                continue; } //J. Tao not e/mu
            if( removeOverlappingCandidates_ &&
                ( ( overlapAlgo_ == 0 &&  TrackWriter::vetoPackedCand( photon, pfcand ) ) ||
                  ( overlapAlgo_ != 0 && ( *overlapAlgo_ )( photon, pfcand ) ) ) ) { 
                std::cout << "SKIPPING CANDIDATE WITH pt=" << pfcand->pt() << " BECAUSE OF OVERLAP" << std::endl;
                continue; }


            if( pfcand->pt() < ptMin )         { 
                std::cout << "SKIPPING CANDIDATE WITH pt=" << pfcand->pt() << " BECAUSE OF MIN PT" << std::endl;
                continue; 
            }
            float dRTkToVtx  = deltaR( pfcand->momentum().Eta(), pfcand->momentum().Phi(),
                                       SCdirection.Eta(), SCdirection.Phi() + deltaPhiRotation_ ); // rotate SC in phi if requested (random cone isolation)
            if( dRTkToVtx > coneSize) { 
                std::cout << "SKIPPING CANDIDATE WITH"
                          << " pt=" << pfcand->pt() 
                          << " eta=" << pfcand->eta() 
                          << " phi=" << pfcand->phi() 
                          << " BECAUSE DELTAR=" << dRTkToVtx 
                          << " OUTSIDE OUTER CONE OF " << coneSize 
                          << " (DETA=" << fabs(pfcand->momentum().Eta() - SCdirection.Eta())
                          << ", DPHI=" << deltaPhi(pfcand->momentum().Phi(), SCdirection.Phi() + deltaPhiRotation_)
                          << ")"
                          << std::endl;
                continue; 
            }

            if( dRTkToVtx < coneVeto ) { 
                std::cout << "SKIPPING CANDIDATE WITH pt=" << pfcand->pt() << " BECAUSE DELTAR=" << dRTkToVtx << " INSIDE INNER CONE OF " << coneVeto << std::endl;
                continue; 
            }

            std::cout << "ACCEPTING CANDIDATE WITH pt=" << pfcand->pt() 
                      << " VTXDZ=" << pfcand->vertex().z() - vtx->z() 
                      << std::endl;

            isovalue += pfcand->pt();
        }
        return isovalue;

    }

    void PhoIdDumper::addPhoton(const edm::EventID &eventId, const flashgg::Photon &photon, 
                                const edm::Ptr<reco::Vertex> &photonVertex,
                                float weight, float mvaID,
                                float chosenVertexChargedIso,
                                float worstVertexChargedIso,
                                const PhoIdMVAInputVars *phoIdInputVars,
                                const flashgg::DiPhotonCandidate &diphoton,
                                const flashgg::VertexCandidateMap &vtxcandmap
                                )
    {
        float label = photon.genMatchType() == flashgg::Photon::kPrompt ? 1 : 0;

        float genDeltaR = photon.genDeltaR();

        std::vector<PhoIdWriter::RecHitData> rechits;

        if (isPhotonInSubdet(photon))
            {
#ifdef MYDEBUG
                {
                    pfIsoChgWrtVtx(photon,
                                   photonVertex,
                                   vtxcandmap,
                                   0.3, // coneSize
                                   0.02, // coneVetoBarrel
                                   0.02, // coneVetoEndcap
                                   0.1 // ptMin
                                   );
                
                }
#endif




                // photon is in the detector region (barrel/endcap) this instance works with
                // TODO: do we need to check this ? We already check in the following function
                //       for each rechit
                fillRecHits(photon, rechits);
                applyWindowAndNormalizeEnergy(rechits, windowHalfWidth, windowHalfHeight, normalizeRecHitsToMax);

                // ignore 'empty' photons for the moment
                if (rechits.size() < 1)
                    return;

                // event identification
                this->runNumber.push_back(eventId.run());
                this->lsNumber.push_back(eventId.luminosityBlock());
                this->eventNumber.push_back(eventId.event());

                this->recHits.push_back(rechits);
                this->weights.push_back(weight);
                this->labels.push_back(label);
                this->mvaID.push_back(mvaID);
                this->genDeltaR.push_back(genDeltaR);

                // track isolation variables
                this->chgIsoWrtChosenVtx.push_back(chosenVertexChargedIso);
                this->chgIsoWrtWorstVtx.push_back(worstVertexChargedIso);

                // photon ID input variables
                if (phoIdInputVars != NULL)
                    {
                        scRawE          .push_back( phoIdInputVars->scRawE          );
                        r9              .push_back( phoIdInputVars->r9              );
                        covIEtaIEta     .push_back( phoIdInputVars->covIEtaIEta     );
                        phiWidth        .push_back( phoIdInputVars->phiWidth        );
                        etaWidth        .push_back( phoIdInputVars->etaWidth        );
                        covIEtaIPhi     .push_back( phoIdInputVars->covIEtaIPhi     );
                        s4              .push_back( phoIdInputVars->s4              );
                        pfPhoIso03      .push_back( phoIdInputVars->pfPhoIso03      );
                        pfChgIso03      .push_back( phoIdInputVars->pfChgIso03      );
                        pfChgIso03worst .push_back( phoIdInputVars->pfChgIso03worst );
                        scEta           .push_back( phoIdInputVars->scEta           );
                        rho             .push_back( phoIdInputVars->rho             );
                        esEffSigmaRR    .push_back( phoIdInputVars->esEffSigmaRR    );

#ifdef MYDEBUG
                        std::cout << "pfChgIso03=" << pfChgIso03.back() << std::endl;
#endif
                    }

                // other photon variables
                photonEt.push_back(photon.et());
                photonEta.push_back(photon.eta());
                photonPhi.push_back(photon.phi());

                photonVertexX.push_back(photon.vertex().x());
                photonVertexY.push_back(photon.vertex().y());
                photonVertexZ.push_back(photon.vertex().z());
                photonVertexIndex.push_back(photonVertex.key());

                photonSCx.push_back(photon.superCluster()->x());
                photonSCy.push_back(photon.superCluster()->y());
                photonSCz.push_back(photon.superCluster()->z());

                diphotonMass.push_back(diphoton.mass());

                // tracks
                trackWriter->addPhoton(photon, photonVertex, vtxcandmap);
                
            }
    }

    //----------------------------------------------------------------------

    PhoIdDumper::PhoIdDumper( const edm::ParameterSet &iConfig, bool isBarrel_):
        isBarrel(isBarrel_),
        diphotonToken_( consumes<edm::View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "diphotonsInput" ) ) ),

        outputFname( iConfig.getUntrackedParameter<std::string>("output") ),
        writeRecHitsSparseFlag( iConfig.getUntrackedParameter<bool>("writeSparse")),
        windowHalfWidth( iConfig.getUntrackedParameter<unsigned>("windowHalfWidth")),
        windowHalfHeight( iConfig.getUntrackedParameter<unsigned>("windowHalfHeight")),

        writePhotonIdInputVarsFlag ( iConfig.getUntrackedParameter<bool>("writePhotonIdInputVars")),

        normalizeRecHitsToMax ( iConfig.getUntrackedParameter<bool>("normalizeRecHitsToMax")),
        
        vertexCandidateMapToken_( consumes<VertexCandidateMap>(InputTag("flashggVertexMapNonUnique")))

    {
        if (writePhotonIdInputVarsFlag)
            {
                phoIdInputVarsToken = consumes<flashgg::DiPhotonPhoIdMVAInputVarsAssociation>( iConfig.getParameter<InputTag> ( "photonIdInputVarsInputTag" ));
            }

        if (boost::algorithm::ends_with(outputFname, ".npz")) {
            phoIdWriter = new PhoIdWriterNumpy();
            trackWriter = new TrackWriterNumpy(iConfig, consumesCollector());
        } else  {
            // assume Torch output format
            phoIdWriter = new PhoIdWriterTorch();
            trackWriter = new TrackWriterTorch(iConfig, consumesCollector());
        }

    }

    //----------------------------------------------------------------------

    PhoIdDumper::~PhoIdDumper()
    {

    }

    //----------------------------------------------------------------------

    void
    PhoIdDumper::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
    {
        Handle<edm::View<flashgg::DiPhotonCandidate> > diphotons;
        iEvent.getByToken( diphotonToken_, diphotons );

        Handle<flashgg::DiPhotonPhoIdMVAInputVarsAssociation> phoIdInputVarsHandle;

        if (writePhotonIdInputVarsFlag)
            iEvent.getByToken(phoIdInputVarsToken, phoIdInputVarsHandle);

        Handle<VertexCandidateMap> vertexCandidateMap;
        iEvent.getByToken( vertexCandidateMapToken_, vertexCandidateMap );
        const flashgg::VertexCandidateMap vtxToCandMap = *( vertexCandidateMap.product() );

        trackWriter->newEvent(iEvent);

        //----------------------------------------

        unsigned diphotonIndex = 0;
        for ( auto diphoton = diphotons.product()->begin(); diphoton != diphotons.product()->end(); ++diphoton, ++diphotonIndex)
            {
                // TODO: should we take the square root of the event weight for photons ?
                // TODO: check if this corresponds to the final event weight ?!
                float weight = diphoton->centralWeight();

                const PhoIdMVAInputVars *phoIdInputVarsLeading = NULL, *phoIdInputVarsSubLeading = NULL;

                // make sure the edm::Ref is not destroyed after the 'if' body
                std::vector<edm::Ptr<flashgg::DiPhotonCandidate>> diphotonPointers = diphotons->ptrs();

                if (writePhotonIdInputVarsFlag)
                    {
                        // do NOT make a copy, we use a pointer instead so that the object
                        // is not destroyed at the end of the 'if' body and the pointers
                        // to the leading and subleading photon variable values remain valid

                        const DiPhotonPhoIdMVAInputVars *diphoInputVars = &(phoIdInputVarsHandle->at(diphotonPointers[diphotonIndex]));
                        phoIdInputVarsLeading = &(diphoInputVars->getInputsLeading());
                        phoIdInputVarsSubLeading = &(diphoInputVars->getInputsSubLeading());
                    }

                addPhoton(iEvent.id(),
                          *diphoton->leadingPhoton(), 
                          diphoton->vtx(),
                          weight, 
                          diphoton->leadingView()->phoIdMvaWrtChosenVtx(),
                          diphoton->leadingView()->pfChIso03WrtChosenVtx(),
                          diphoton->leadingPhoton()->pfChgIsoWrtWorstVtx04(),
                          phoIdInputVarsLeading,
                          *diphoton,
                          vtxToCandMap
                          );
                addPhoton(iEvent.id(),
                          *diphoton->subLeadingPhoton(), 
                          diphoton->vtx(), 
                          weight, 
                          diphoton->subLeadingView()->phoIdMvaWrtChosenVtx(),
                          diphoton->subLeadingView()->pfChIso03WrtChosenVtx(),
                          diphoton->subLeadingPhoton()->pfChgIsoWrtWorstVtx04(),
                          phoIdInputVarsSubLeading,
                          *diphoton,
                          vtxToCandMap
                          );

                // only consider the first pair (how are they sorted ?)
                break;
            } // loop over diphotons
    } // analyze

    //----------------------------------------------------------------------

    void
    PhoIdDumper::beginJob()
    {
    }

    //----------------------------------------------------------------------

    void
    PhoIdDumper::endJob()
    {
        phoIdWriter->writeTo(*this, outputFname);
    }

    //----------------------------------------------------------------------

    void
    PhoIdDumper::fillDescriptions( edm::ConfigurationDescriptions &descriptions )
    {
        //The following says we do not know what parameters are allowed so do no validation
        // Please change this to state exactly what you do use, even if it is no parameters
        edm::ParameterSetDescription desc;
        desc.setUnknown();
        descriptions.addDefault( desc );
    }

    //----------------------------------------------------------------------
    // PhoIdDumperBarrel
    //----------------------------------------------------------------------

    PhoIdDumperBarrel::PhoIdDumperBarrel( const edm::ParameterSet &iParams) : 
        PhoIdDumper(iParams, true)
    {
    }

    //----------------------------------------------------------------------

    void PhoIdDumperBarrel::wrapCoordinates(PhoIdWriter::RecHitData &rechit)
    {
        while (rechit.dx < 0)
            rechit.dx += 360;
            
        while (rechit.dx >= 360)
            rechit.dx -= 360;
    }

    //----------------------------------------------------------------------

    void PhoIdDumperBarrel::fillRecHits(const flashgg::Photon &photon, std::vector<PhoIdWriter::RecHitData> &rechits)
    {
        for (auto rechit : *photon.recHits())
            {
                if (rechit.detid().det() != DetId::Ecal)
                    continue;

                if (rechit.detid().subdetId() != EcalBarrel)
                    continue;

                EBDetId dt = EBDetId(rechit.detid());

                PhoIdWriter::RecHitData data;
                data.energy = rechit.energy();
                data.dx = dt.ieta();
                data.dy = dt.iphi();

                // eta goes from -85 to -1 then jumps to +1 to +85
                // (i.e. there is no zero..)
                if (data.dx < 0)
                    data.dx += 1;

                // iphi goes from 1 to 360, move it to 0.. 359
                data.dy -= 1;
                
                rechits.push_back(data);

            } // loop over rechits
    }

    //----------------------------------------------------------------------

    bool PhoIdDumperBarrel::isPhotonInSubdet(const flashgg::Photon &photon)
    {
        return fabs(photon.eta()) < etaMaxBarrel;
    }

    //----------------------------------------------------------------------
    // PhoIdDumperEndcap
    //----------------------------------------------------------------------

    PhoIdDumperEndcap::PhoIdDumperEndcap( const edm::ParameterSet &iParams) : 
        PhoIdDumper(iParams, false)
    {
    }

    //----------------------------------------------------------------------

    void PhoIdDumperEndcap::wrapCoordinates(PhoIdWriter::RecHitData &rechit)
    {
        // nothing to do in the endcap
    }

    //----------------------------------------------------------------------

    void PhoIdDumperEndcap::fillRecHits(const flashgg::Photon &photon, std::vector<PhoIdWriter::RecHitData> &rechits)
    {
        for (auto rechit : *photon.recHits())
            {
                if (rechit.detid().det() != DetId::Ecal)
                    continue;

                if (rechit.detid().subdetId() != EcalEndcap)
                    continue;

                EEDetId dt = EEDetId(rechit.detid());

                PhoIdWriter::RecHitData data;
                data.energy = rechit.energy();
                data.dx = dt.ix();
                data.dy = dt.iy();

                rechits.push_back(data);

            } // loop over rechits
    }

    //----------------------------------------------------------------------

    bool PhoIdDumperEndcap::isPhotonInSubdet(const flashgg::Photon &photon)
    {
        return fabs(photon.eta()) >= etaMaxBarrel;
    }

    //----------------------------------------------------------------------

} // namespace flashgg

typedef flashgg::PhoIdDumperBarrel FlashggPhoIdDumperBarrel;
DEFINE_FWK_MODULE( FlashggPhoIdDumperBarrel );

typedef flashgg::PhoIdDumperEndcap FlashggPhoIdDumperEndcap;
DEFINE_FWK_MODULE( FlashggPhoIdDumperEndcap );


// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

