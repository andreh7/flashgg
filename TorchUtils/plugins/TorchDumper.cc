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

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"


using namespace std;
using namespace edm;

// **********************************************************************

namespace flashgg {

    const unsigned MAGIC_NUMBER = 1;
    const unsigned MAGIC_STRING = 2;
    const unsigned MAGIC_TABLE  = 3;
    const unsigned MAGIC_TORCH  = 4;

    template<typename DataType>
    inline void writeType(std::ostream &os, const DataType &value)
    {
        os.write((const char*) &value, sizeof(value));
    }

    inline void writeInt(std::ostream &os, int32_t value)
    {
        os.write((const char*) &value, sizeof(value));
    }

    /** assumes that double is 8 bytes */
    inline void writeDouble(std::ostream &os, double value)
    {
        os.write((const char*) &value, sizeof(value));
    }

    inline void writeLong(std::ostream &os, int64_t value)
    {
        os.write((const char*)&value, sizeof(value));
    }

    inline void writeFloat(std::ostream &os, float value)
    {
        os.write((const char*)&value, sizeof(value));
    }


    inline void writeString(std::ostream &os, const std::string &str)
    {
        size_t len = str.size();
        writeInt(os, len);

        os.write((const char*) &str[0], len);
    }

    template<typename DataType>
    void writeTypeTensorHelper(std::ostream &os, unsigned &objectIndex, const std::vector<unsigned> &sizes, const std::vector<DataType> &data, const std::string &tensorTypeName,
                         const std::string &storageTypeName)
    {
        // data must be stored in the order such the an increase of the index in the last dimension
        // by one corresponds to an indcreas of the index into data by one etc.

        writeInt(os, MAGIC_TORCH);

        // write object index (we do NOT reuse objects here, we assume that we only write one
        // object per file)
        writeInt(os, objectIndex++);

        // version string
        writeString(os, "V 1");

        // class name
        writeString(os, tensorTypeName);

        //----------

        // write the number of coordinates
        writeInt(os, sizes.size());

        // write the sizes of each dimension
        for (unsigned dimsize : sizes)
            writeLong(os, dimsize);

        // calculate strides
        std::vector<uint64_t> strides(sizes.size());
        {
            uint64_t product = 1;
            for (unsigned i = sizes.size(); i > 0; --i)
                {
                    strides[i-1] = product;
                    product *= sizes[i-1];
                }

            assert(product == data.size());
        }
        for (uint64_t stride : strides)
            writeLong(os, stride);

        // write storage offset
        writeLong(os, 1);

        //----------
        // write the FloatStorage object: the actual data
        //----------
        writeInt(os, MAGIC_TORCH);

        // write object index (we do NOT reuse objects here, we assume that we only write one
        // object per file)
        writeInt(os, objectIndex++);

        // version string
        writeString(os, "V 1");

        // class name
        writeString(os, storageTypeName);

        // size
        writeLong(os, data.size());

        // the actual content
        for (unsigned i = 0 ; i < data.size(); ++i)
            writeType<DataType>(os, data[i]);

    }

    template<typename DataType>
    void writeTypeTensor(std::ostream &os, unsigned &objectIndex, const std::vector<unsigned> &sizes, const std::vector<DataType> &data);

    // template specializations
    template<>
    inline void writeTypeTensor<float>(std::ostream &os, unsigned &objectIndex, const std::vector<unsigned> &sizes, const std::vector<float> &data) 
    {
        writeTypeTensorHelper<float>(os, objectIndex, sizes, data, "torch.FloatTensor", "torch.FloatStorage");
    }

    template<>
    inline void writeTypeTensor<int32_t>(std::ostream &os, unsigned &objectIndex, const std::vector<unsigned> &sizes, const std::vector<int32_t> &data) 
    {
        writeTypeTensorHelper<int32_t>(os, objectIndex, sizes, data, "torch.IntTensor", "torch.IntStorage");
    }

    //----------------------------------------------------------------------

    class TorchDumper : public edm::EDAnalyzer
    {
    public:
        explicit TorchDumper( const edm::ParameterSet & );
        virtual ~TorchDumper();

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

        /** boundary between barrel and endcap */
        const float etaMaxBarrel = 1.5;

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

        void addPhoton(const flashgg::Photon &photon, float weight, float mvaID,
                       float chosenVertexChargedIso,
                       float worstVertexChargedIso
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
            }
        }

        void writeRecHits(std::ostream &os, unsigned &objectIndex, const std::vector<std::vector<RecHitData> > &rechits, 
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

            writeTypeTensor(os, objectIndex, sizes, data);
        }

        /** function to write out one element of each rechit as a table of tables */
        template<typename DataType, typename Func>
        void writeRecHitsValues(std::ostream &os, unsigned &objectIndex, 
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

            writeTypeVector<DataType>(os, objectIndex, values);
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
        void writeRecHitsSparse(std::ostream &os, unsigned &objectIndex, const std::vector<std::vector<RecHitData> > &rechits)
        {
            const unsigned tableSize = 5;

            writeInt(os, MAGIC_TABLE);
            writeInt(os, objectIndex++);
            writeInt(os, tableSize);

            std::vector<int32_t> firstIndex(rechits.size()), numRecHits(rechits.size());
            unsigned nextStartIndex = 1;

            for (unsigned i = 0; i < rechits.size(); ++i)
            {
                firstIndex[i] = nextStartIndex;
                numRecHits[i] = rechits[i].size();
                
                nextStartIndex += rechits[i].size();

            } // loop over photons
            
            unsigned totNumRecHits = nextStartIndex - 1;
            
            writeInt(os, MAGIC_STRING); writeString(os, "firstIndex");  writeTypeVector(os, objectIndex, firstIndex);
            writeInt(os, MAGIC_STRING); writeString(os, "numRecHits");  writeTypeVector(os, objectIndex, numRecHits);

            // we add one to the x and y indices to stick to torch's one based indexing
            // (note that in the dense writing routine this is not needed because
            // we calculate the address directly)
            writeInt(os, MAGIC_STRING); writeString(os, "energy");  writeRecHitsValues<float>(os, objectIndex, rechits, totNumRecHits, [](const RecHitData &rh)  { return rh.energy; });
            writeInt(os, MAGIC_STRING); writeString(os, "x");       writeRecHitsValues<int32_t>(os, objectIndex, rechits, totNumRecHits, [](const RecHitData &rh)  { return rh.dx + 1; });
            writeInt(os, MAGIC_STRING); writeString(os, "y");       writeRecHitsValues<int32_t>(os, objectIndex, rechits, totNumRecHits, [](const RecHitData &rh)  { return rh.dy + 1; });
        }


        template<typename DataType>
        void writeTypeVector(std::ostream &os, unsigned &objectIndex, const std::vector<DataType> &data)
        {
            std::vector<unsigned> sizes;
            sizes.push_back(data.size());

            writeTypeTensor<DataType>(os, objectIndex, sizes, data);
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
            const unsigned tableSize = 7;

            std::ofstream os(outputFname.c_str());

            writeInt(os, MAGIC_TABLE);

            unsigned objectIndex = 1;
            writeInt(os, objectIndex++);
            
            writeInt(os, tableSize);

            writeInt(os, MAGIC_STRING); writeString(os, "X");      
            if (writeRecHitsSparseFlag)
                writeRecHitsSparse(os, objectIndex, recHits);
            else
                writeRecHits(os, objectIndex, recHits, windowHalfWidth, windowHalfHeight);

            writeInt(os, MAGIC_STRING); writeString(os, "y");      writeTypeVector(os, objectIndex, labels);
            writeInt(os, MAGIC_STRING); writeString(os, "weight"); writeTypeVector(os, objectIndex, weights);
            writeInt(os, MAGIC_STRING); writeString(os, "mvaid");  writeTypeVector(os, objectIndex, mvaID);
            writeInt(os, MAGIC_STRING); writeString(os, "genDR");  writeTypeVector(os, objectIndex, genDeltaR);

            writeInt(os, MAGIC_STRING); writeString(os, "chgIsoWrtChosenVtx");  writeTypeVector(os, objectIndex, chgIsoWrtChosenVtx);
            writeInt(os, MAGIC_STRING); writeString(os, "chgIsoWrtWorstVtx");   writeTypeVector(os, objectIndex, chgIsoWrtWorstVtx);
        }

    };

// ******************************************************************************************
// ******************************************************************************************

    TorchDumper::TorchDumper( const edm::ParameterSet &iConfig ):
        diphotonToken_( consumes<edm::View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "diphotonsInput" ) ) ),

        outputFname( iConfig.getUntrackedParameter<std::string>("output") ),
        writeRecHitsSparseFlag( iConfig.getUntrackedParameter<bool>("writeSparse")),
        windowHalfWidth( iConfig.getUntrackedParameter<unsigned>("windowHalfWidth")),
        windowHalfHeight( iConfig.getUntrackedParameter<unsigned>("windowHalfHeight"))
    {
    }

    TorchDumper::~TorchDumper()
    {

    }

    void
    TorchDumper::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
    {
        Handle<edm::View<flashgg::DiPhotonCandidate> > diphotons;
        iEvent.getByToken( diphotonToken_, diphotons );

        //for ( auto photon = photons.product()->begin(); photon != photons.product()->end(); ++photon)
        for ( auto diphoton : *diphotons.product()) {
            // TODO: should we take the square root of the event weight for photons ?
            // TODO: check if this corresponds to the final event weight ?!
            float weight = diphoton.centralWeight();

            addPhoton(*diphoton.leadingPhoton(), weight, 
                      diphoton.leadingView()->phoIdMvaWrtChosenVtx(),
                      diphoton.leadingView()->pfChIso03WrtChosenVtx(),
                      diphoton.leadingPhoton()->pfChgIsoWrtWorstVtx04()
                      );
            addPhoton(*diphoton.subLeadingPhoton(), weight, 
                      diphoton.subLeadingView()->phoIdMvaWrtChosenVtx(),
                      diphoton.subLeadingView()->pfChIso03WrtChosenVtx(),
                      diphoton.subLeadingPhoton()->pfChgIsoWrtWorstVtx04()

                      );

            // only consider the first pair (how are they sorted ?)
            break;
        }
    } // analyze

    void
    TorchDumper::beginJob()
    {
    }

    void
    TorchDumper::endJob()
    {
        writeTorchData();
    }

    void
    TorchDumper::fillDescriptions( edm::ConfigurationDescriptions &descriptions )
    {
        //The following says we do not know what parameters are allowed so do no validation
        // Please change this to state exactly what you do use, even if it is no parameters
        edm::ParameterSetDescription desc;
        desc.setUnknown();
        descriptions.addDefault( desc );
    }

    //----------------------------------------------------------------------

    /** dumps endcap photons */
    class TorchDumperBarrel : public TorchDumper
    {
        //----------------------------------------
    public:
        explicit TorchDumperBarrel( const edm::ParameterSet &iParams) : 
            TorchDumper(iParams)
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
    class TorchDumperEndcap : public TorchDumper
    {
        //----------------------------------------
    public:

        explicit TorchDumperEndcap( const edm::ParameterSet &iParams) : 
            TorchDumper(iParams)
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

typedef flashgg::TorchDumperBarrel FlashggTorchDumperBarrel;
DEFINE_FWK_MODULE( FlashggTorchDumperBarrel );

typedef flashgg::TorchDumperEndcap FlashggTorchDumperEndcap;
DEFINE_FWK_MODULE( FlashggTorchDumperEndcap );


// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

