#include <memory>

#include <vector>
#include <fstream>
#include <cassert>

#include <cstdint>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include "flashgg/DataFormats/interface/VBFTag.h"
#include "flashgg/DataFormats/interface/UntaggedTag.h"
#include "flashgg/DataFormats/interface/DiPhotonTagBase.h"
#include "flashgg/DataFormats/interface/TTHHadronicTag.h"
#include "flashgg/DataFormats/interface/TTHLeptonicTag.h"
#include "flashgg/DataFormats/interface/VHTightTag.h"
#include "flashgg/DataFormats/interface/VHEtTag.h"
#include "flashgg/DataFormats/interface/VHLooseTag.h"
#include "flashgg/DataFormats/interface/VHHadronicTag.h"
#include "flashgg/DataFormats/interface/VBFTagTruth.h"
#include "flashgg/DataFormats/interface/ZPlusJetTag.h"
#include "flashgg/DataFormats/interface/DiPhotonPhoIdMVAInputVars.h"

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "flashgg/TorchUtils/interface/TorchWriter.h"
#include "flashgg/TorchUtils/interface/TrackWriter.h"


using namespace std;
using namespace edm;

// **********************************************************************

namespace flashgg {

    //----------------------------------------------------------------------

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

        struct RecHitData
        {
            float energy;
            int dx; // eta in the barrel
            int dy; // phi in the barrel
        };


        /** rechit information to be written out at the end */
        std::vector<std::vector<RecHitData> > recHits;
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

        TrackWriter trackWriter;

        //----------------------------------------

        virtual void wrapCoordinates(RecHitData &rechit) = 0;

        /** finds the crystal with the maximum energy, normalizes to that and applies the given window around
            the maximum */
        void applyWindowAndNormalizeEnergy(std::vector<RecHitData> &rechits, int windowHalfWidth, int windowHalfHeight)
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
                if (maxEnergy > 0)
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

        /** subdet specific: extract rechit coordinates */
        virtual void fillRecHits(const flashgg::Photon &photon, std::vector<RecHitData> &rechits) = 0;

        /** called to check if a photon is in barrel/endcap */
        virtual bool isPhotonInSubdet(const flashgg::Photon &photon) = 0;

        void addPhoton(const flashgg::Photon &photon, 
                       const edm::Ptr<reco::Vertex> &photonVertex,
                       float weight, float mvaID,
                       float chosenVertexChargedIso,
                       float worstVertexChargedIso,
                       const PhoIdMVAInputVars *phoIdInputVars
                       )
        {
            float label = photon.genMatchType() == flashgg::Photon::kPrompt ? 1 : 0;

            float genDeltaR = photon.genDeltaR();

            std::vector<RecHitData> rechits;

            if (isPhotonInSubdet(photon))
            {
                // photon is in the detector region (barrel/endcap) this instance works with
                // TODO: do we need to check this ? We already check in the following function
                //       for each rechit
                fillRecHits(photon, rechits);
                applyWindowAndNormalizeEnergy(rechits, windowHalfWidth, windowHalfHeight);

                // ignore 'empty' photons for the moment
                if (rechits.size() < 1)
                    return;

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
                    }

                // other photon variables
                photonEt.push_back(photon.et());

                // tracks
                trackWriter.addPhoton(photon, photonVertex);
                
            }
        }

        void writeRecHits(TorchWriter &tw, const std::vector<std::vector<RecHitData> > &rechits, 
                          int windowHalfWidth, int windowHalfHeight)
        {
            std::vector<unsigned> sizes;
            sizes.push_back(rechits.size());
            sizes.push_back(1);
            sizes.push_back(2 * windowHalfWidth + 1);
            sizes.push_back(2 * windowHalfHeight + 1);

            // fill the rechits (relative) energies
            std::vector<float> data(sizes[0] * sizes[1] * sizes[2] * sizes[3], 0.);
            for (unsigned event = 0; event < rechits.size(); ++event)
            {
                unsigned baseIndex = event * sizes[1] * sizes[2] * sizes[3];

                for (const RecHitData rechit : rechits[event])
                {
                    unsigned offset = rechit.dx * sizes[3] + rechit.dy;
                    data[baseIndex + offset] = rechit.energy;
                }

            } // loop over events

            tw.writeTypeTensor(sizes, data);
        }

        /** function to write out one element of each rechit as a table of tables */
        template<typename DataType, typename Func>
        void writeRecHitsValues(TorchWriter &tw, 
                                const std::vector<std::vector<RecHitData> > &rechits,
                                unsigned totNumRecHits, 
                                Func&& valueFunc)
        {
            std::vector<DataType> values(totNumRecHits);

            unsigned nextPos = 0;

            for (unsigned i = 0; i < rechits.size(); ++i)
            {
                for (unsigned j = 0; j < rechits[i].size(); ++j)
                    values[nextPos++] = valueFunc(rechits[i][j]);
            } // loop over photons

            cout << "nextPos=" << nextPos << " totNumRecHits=" << totNumRecHits << endl;
            assert(nextPos == totNumRecHits);

            tw.writeTypeVector<DataType>(values);
        }
                                

        /** writes the rechits as a list of list of rechits rather than a tensor which has many zeros in it.
            Assumes that the given window already has been applied. 

            writes a table with the following indices:
              energy     : 1d tensor of (normalized) rechit energies
              x          : 1d tensor of x coordinates (written out as one based)
              y          : 1d tensor of y coordinates (written out as one based)
              firstIndex : 1d tensor of first indices in the energy, x and y tensors (index of this
                           tensor the index of the photon)
              numRecHits : number of rechits for the given photon
        */
        void writeRecHitsSparse(TorchWriter &tw, const std::vector<std::vector<RecHitData> > &rechits)
        {
            const unsigned tableSize = 5;

            tw.writeInt(tw.MAGIC_TABLE);
            tw.writeInt(tw.getNextObjectIndex());
            tw.writeInt(tableSize);

            std::vector<int32_t> firstIndex(rechits.size()), numRecHits(rechits.size());
            unsigned nextStartIndex = 1;

            for (unsigned i = 0; i < rechits.size(); ++i)
            {
                firstIndex[i] = nextStartIndex;
                numRecHits[i] = rechits[i].size();
                
                nextStartIndex += rechits[i].size();

            } // loop over photons
            
            unsigned totNumRecHits = nextStartIndex - 1;
            
            tw.writeInt(tw.MAGIC_STRING); tw.writeString("firstIndex");  tw.writeTypeVector(firstIndex);
            tw.writeInt(tw.MAGIC_STRING); tw.writeString("numRecHits");  tw.writeTypeVector(numRecHits);

            // we add one to the x and y indices to stick to torch's one based indexing
            // (note that in the dense writing routine this is not needed because
            // we calculate the address directly)
            tw.writeInt(tw.MAGIC_STRING); tw.writeString("energy");  writeRecHitsValues<float>(tw, rechits, totNumRecHits, [](const RecHitData &rh)  { return rh.energy; });
            tw.writeInt(tw.MAGIC_STRING); tw.writeString("x");       writeRecHitsValues<int32_t>(tw, rechits, totNumRecHits, [](const RecHitData &rh)  { return rh.dx + 1; });
            tw.writeInt(tw.MAGIC_STRING); tw.writeString("y");       writeRecHitsValues<int32_t>(tw, rechits, totNumRecHits, [](const RecHitData &rh)  { return rh.dy + 1; });
        }


        void writeTorchData()
        {
            // write a dict/table with 
            //  X = rechits
            //  y = labels
            //  weight = weights
            //  mvaid (for comparison)
            //  charged isolation w.r.t. chosen vertex
            //  charged isolation w.r.t. worst vertex
            unsigned tableSize = 7;

            if (writePhotonIdInputVarsFlag)
                // add sub-recrod with photon id input variables
                tableSize += 1;
            
            // for writing out tracks
            tableSize += 1;

            // for writing other photon variables
            tableSize += 1;

            std::ofstream os(outputFname.c_str());
            TorchWriter tw(os);

            tw.writeInt(tw.MAGIC_TABLE);

            tw.writeInt(tw.getNextObjectIndex());
            
            tw.writeInt(tableSize);

            tw.writeInt(tw.MAGIC_STRING); tw.writeString("X");
            if (writeRecHitsSparseFlag)
                writeRecHitsSparse(tw, recHits);
            else
                writeRecHits(tw, recHits, windowHalfWidth, windowHalfHeight);

            tw.writeInt(tw.MAGIC_STRING); tw.writeString("y");      tw.writeTypeVector(labels);
            tw.writeInt(tw.MAGIC_STRING); tw.writeString("weight"); tw.writeTypeVector(weights);
            tw.writeInt(tw.MAGIC_STRING); tw.writeString("mvaid");  tw.writeTypeVector(mvaID);
            tw.writeInt(tw.MAGIC_STRING); tw.writeString("genDR");  tw.writeTypeVector(genDeltaR);

            tw.writeInt(tw.MAGIC_STRING); tw.writeString("chgIsoWrtChosenVtx");  tw.writeTypeVector(chgIsoWrtChosenVtx);
            tw.writeInt(tw.MAGIC_STRING); tw.writeString("chgIsoWrtWorstVtx");   tw.writeTypeVector(chgIsoWrtWorstVtx);

            if (writePhotonIdInputVarsFlag)
                {
                    tw.writeInt(tw.MAGIC_STRING); tw.writeString("phoIdInput");      

                    tw.writeInt(tw.MAGIC_TABLE);
                    tw.writeInt(tw.getNextObjectIndex());
                    tw.writeInt(13); // number of input variables

                    tw.writeInt(tw.MAGIC_STRING); tw.writeString("scRawE"          );  tw.writeTypeVector(scRawE          );
                    tw.writeInt(tw.MAGIC_STRING); tw.writeString("r9"              );  tw.writeTypeVector(r9              );
                    tw.writeInt(tw.MAGIC_STRING); tw.writeString("covIEtaIEta"     );  tw.writeTypeVector(covIEtaIEta     );
                    tw.writeInt(tw.MAGIC_STRING); tw.writeString("phiWidth"        );  tw.writeTypeVector(phiWidth        );
                    tw.writeInt(tw.MAGIC_STRING); tw.writeString("etaWidth"        );  tw.writeTypeVector(etaWidth        );
                    tw.writeInt(tw.MAGIC_STRING); tw.writeString("covIEtaIPhi"     );  tw.writeTypeVector(covIEtaIPhi     );
                    tw.writeInt(tw.MAGIC_STRING); tw.writeString("s4"              );  tw.writeTypeVector(s4              );
                    tw.writeInt(tw.MAGIC_STRING); tw.writeString("pfPhoIso03"      );  tw.writeTypeVector(pfPhoIso03      );
                    tw.writeInt(tw.MAGIC_STRING); tw.writeString("pfChgIso03"      );  tw.writeTypeVector(pfChgIso03      );
                    tw.writeInt(tw.MAGIC_STRING); tw.writeString("pfChgIso03worst" );  tw.writeTypeVector(pfChgIso03worst );
                    tw.writeInt(tw.MAGIC_STRING); tw.writeString("scEta"           );  tw.writeTypeVector(scEta           );
                    tw.writeInt(tw.MAGIC_STRING); tw.writeString("rho"             );  tw.writeTypeVector(rho             );
                    tw.writeInt(tw.MAGIC_STRING); tw.writeString("esEffSigmaRR"    );  tw.writeTypeVector(esEffSigmaRR    );
                }

            //----------
            // other photon variables (e.g. Et for Et flattening)
            //----------
            {
                tw.writeInt(tw.MAGIC_STRING); tw.writeString("phoVars");
                tw.writeInt(tw.MAGIC_TABLE);
                tw.writeInt(tw.getNextObjectIndex());
                tw.writeInt(1); // number of variables

                tw.writeInt(tw.MAGIC_STRING); tw.writeString("phoEt"          );  tw.writeTypeVector(photonEt);
            }

            //----------


            tw.writeInt(tw.MAGIC_STRING); tw.writeString("tracks");      
            trackWriter.writeOut(tw);
        }

    };

// ******************************************************************************************
// ******************************************************************************************

    PhoIdDumper::PhoIdDumper( const edm::ParameterSet &iConfig ):
        diphotonToken_( consumes<edm::View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "diphotonsInput" ) ) ),

        outputFname( iConfig.getUntrackedParameter<std::string>("output") ),
        writeRecHitsSparseFlag( iConfig.getUntrackedParameter<bool>("writeSparse")),
        windowHalfWidth( iConfig.getUntrackedParameter<unsigned>("windowHalfWidth")),
        windowHalfHeight( iConfig.getUntrackedParameter<unsigned>("windowHalfHeight")),

        writePhotonIdInputVarsFlag ( iConfig.getUntrackedParameter<bool>("writePhotonIdInputVars")),

        trackWriter(iConfig, consumesCollector())

    {
        if (writePhotonIdInputVarsFlag)
        {
            phoIdInputVarsToken = consumes<flashgg::DiPhotonPhoIdMVAInputVarsAssociation>( iConfig.getParameter<InputTag> ( "photonIdInputVarsInputTag" ));
        }

    }

    PhoIdDumper::~PhoIdDumper()
    {

    }

    void
    PhoIdDumper::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
    {
        Handle<edm::View<flashgg::DiPhotonCandidate> > diphotons;
        iEvent.getByToken( diphotonToken_, diphotons );

        Handle<flashgg::DiPhotonPhoIdMVAInputVarsAssociation> phoIdInputVarsHandle;

        if (writePhotonIdInputVarsFlag)
            iEvent.getByToken(phoIdInputVarsToken, phoIdInputVarsHandle);

        trackWriter.newEvent(iEvent);

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

            addPhoton(*diphoton->leadingPhoton(), 
                      diphoton->vtx(),
                      weight, 
                      diphoton->leadingView()->phoIdMvaWrtChosenVtx(),
                      diphoton->leadingView()->pfChIso03WrtChosenVtx(),
                      diphoton->leadingPhoton()->pfChgIsoWrtWorstVtx04(),
                      phoIdInputVarsLeading
                      );
            addPhoton(*diphoton->subLeadingPhoton(), 
                      diphoton->vtx(), 
                      weight, 
                      diphoton->subLeadingView()->phoIdMvaWrtChosenVtx(),
                      diphoton->subLeadingView()->pfChIso03WrtChosenVtx(),
                      diphoton->subLeadingPhoton()->pfChgIsoWrtWorstVtx04(),
                      phoIdInputVarsSubLeading
                      );

            // only consider the first pair (how are they sorted ?)
            break;
        } // loop over diphotons
    } // analyze

    void
    PhoIdDumper::beginJob()
    {
    }

    void
    PhoIdDumper::endJob()
    {
        writeTorchData();
    }

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

    /** dumps endcap photons */
    class PhoIdDumperBarrel : public PhoIdDumper
    {
        //----------------------------------------
    public:
        explicit PhoIdDumperBarrel( const edm::ParameterSet &iParams) : 
            PhoIdDumper(iParams)
        {
        }

        //----------------------------------------

    protected:
        // wrap around in phi for the barrel
        virtual void wrapCoordinates(RecHitData &rechit)
        {
            while (rechit.dx < 0)
                rechit.dx += 360;
            
            while (rechit.dx >= 360)
                rechit.dx -= 360;
        }

        //----------------------------------------

        virtual void fillRecHits(const flashgg::Photon &photon, std::vector<RecHitData> &rechits)
        {
            for (auto rechit : *photon.recHits())
            {
                if (rechit.detid().det() != DetId::Ecal)
                    continue;

                if (rechit.detid().subdetId() != EcalBarrel)
                    continue;

                EBDetId dt = EBDetId(rechit.detid());

                RecHitData data;
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

        //----------------------------------------

        virtual bool isPhotonInSubdet(const flashgg::Photon &photon)
        {
            return fabs(photon.eta()) < etaMaxBarrel;
        }

    };

    //----------------------------------------------------------------------

    /** dumps barrel photons */
    class PhoIdDumperEndcap : public PhoIdDumper
    {
        //----------------------------------------
    public:

        explicit PhoIdDumperEndcap( const edm::ParameterSet &iParams) : 
            PhoIdDumper(iParams)
        {
        }

        //----------------------------------------

    protected:
        virtual void wrapCoordinates(RecHitData &rechit)
        {
            // nothing to do in the endcap
        }

        //----------------------------------------

        void fillRecHits(const flashgg::Photon &photon, std::vector<RecHitData> &rechits)
        {
            for (auto rechit : *photon.recHits())
            {
                if (rechit.detid().det() != DetId::Ecal)
                    continue;

                if (rechit.detid().subdetId() != EcalEndcap)
                    continue;

                EEDetId dt = EEDetId(rechit.detid());

                RecHitData data;
                data.energy = rechit.energy();
                data.dx = dt.ix();
                data.dy = dt.iy();

                rechits.push_back(data);

            } // loop over rechits
        }

        //----------------------------------------

        virtual bool isPhotonInSubdet(const flashgg::Photon &photon)
        {
            return fabs(photon.eta()) >= etaMaxBarrel;
        }
        //----------------------------------------


    };

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

