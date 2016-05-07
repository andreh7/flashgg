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

    const unsigned MAGIC_TORCH = 4;


    inline void writeInt(std::ostream &os, int32_t value)
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

    void writeFloatTensor(std::ostream &os, const std::vector<unsigned> &sizes, const std::vector<float> &data)
    {
        // data must be stored in the order such the an increase of the index in the last dimension
        // by one corresponds to an indcreas of the index into data by one etc.

        unsigned objectIndex = 1;

        writeInt(os, MAGIC_TORCH);

        // write object index (we do NOT reuse objects here, we assume that we only write one
        // object per file)
        writeInt(os, objectIndex++);

        // version string
        writeString(os, "V 1");

        // class name
        writeString(os, "torch.FloatTensor");

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
        writeString(os, "torch.FloatStorage");

        // size
        writeLong(os, data.size());

        // the actual content
        for (unsigned i = 0 ; i < data.size(); ++i)
            writeFloat(os, data[i]);

    }

    //----------------------------------------------------------------------

    class TorchDumper : public edm::EDAnalyzer
    {
    public:
        explicit TorchDumper( const edm::ParameterSet & );
        ~TorchDumper();

        static void fillDescriptions( edm::ConfigurationDescriptions &descriptions );


    private:

        virtual void beginJob() override;
        virtual void analyze( const edm::Event &, const edm::EventSetup & ) override;
        virtual void endJob() override;

        edm::EDGetTokenT<edm::View<flashgg::Photon> > photonToken_;
        
        /** window sizes */
        const unsigned barrelWindowHalfWidth = 3;
        const unsigned barrelWindowHalfHeight = 11;

        const unsigned endcapWindowHalfWidth = 3;
        const unsigned endcapWindowHalfHeight = 11;

        struct RecHitData
        {
            float energy;
            int dx; // eta in the barrel
            int dy; // phi in the barrel
        };


        /** rechit information to be written out at the end */
        std::vector<std::vector<RecHitData> > barrelRecHits;
        std::vector<float> barrelWeights;
        std::vector<float> barrelLabels; // 1 = prompt photon, 0 = fake or non-prompt photon

        std::vector<std::vector<RecHitData> > endcapRecHits;
        std::vector<float> endcapWeights;
        std::vector<float> endcapLabels; // 1 = prompt photon, 0 = fake or non-prompt photon

        //----------------------------------------

        /** finds the crystal with the maximum energy, normalizes to that and applies the given window around
            the maximum */
        void applyWindowAndNormalizeEnergy(std::vector<RecHitData> &rechits, int windowHalfWidth, int windowHalfHeight, bool isBarrel)
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
                if (isBarrel)
                {
                    while (rechits[i].dx < 0)
                        rechits[i].dx += 360;

                    while (rechits[i].dx >= 360)
                        rechits[i].dx -= 360;
                }

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

        void fillBarrelRecHits(const flashgg::Photon &photon, std::vector<RecHitData> &rechits)
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

        void fillEndcapRecHits(const flashgg::Photon &photon, std::vector<RecHitData> &rechits)
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


        void addPhoton(const flashgg::Photon &photon, float weight)
        {
            float label = photon.genMatchType() == flashgg::Photon::kPrompt ? 1 : 0;

            std::vector<RecHitData> rechits;

            const float etaMaxBarrel = 1.5;

            if (fabs(photon.eta()) < etaMaxBarrel)
            {
                // barrel
                fillBarrelRecHits(photon, rechits);
                applyWindowAndNormalizeEnergy(rechits, barrelWindowHalfWidth, barrelWindowHalfHeight, true);

                // ignore 'empty' photons for the moment
                if (rechits.size() < 1)
                    return;

                barrelRecHits.push_back(rechits);
                barrelWeights.push_back(weight);
                barrelLabels.push_back(label);
            }
            else
            {
                // endcap
                fillEndcapRecHits(photon, rechits);
                applyWindowAndNormalizeEnergy(rechits, endcapWindowHalfWidth, endcapWindowHalfHeight, false);

                // ignore 'empty' photons for the moment
                if (rechits.size() < 1)
                    return;

                endcapRecHits.push_back(rechits);
                endcapWeights.push_back(weight);
                endcapLabels.push_back(label);

            }
        }

        void writeRecHits(const std::string &fname, const std::vector<std::vector<RecHitData> > &rechits, int windowHalfWidth, int windowHalfHeight)
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

            std::ofstream fout(fname.c_str());

            writeFloatTensor(fout, sizes, data);
        }

        /** use this e.g. for labels and weights */
        void writeFloatVector(const std::string &fname, const std::vector<float> &data)
        {
            std::vector<unsigned> sizes;
            sizes.push_back(data.size());
            std::ofstream fout(fname.c_str());

            writeFloatTensor(fout, sizes, data);
        }
    };

// ******************************************************************************************
// ******************************************************************************************

    TorchDumper::TorchDumper( const edm::ParameterSet &iConfig ):
        photonToken_( consumes<edm::View<flashgg::Photon> >( iConfig.getParameter<InputTag> ( "photons" ) ) )
    {
    }

    TorchDumper::~TorchDumper()
    {

    }

    void
    TorchDumper::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
    {
        Handle<edm::View<flashgg::Photon> > photons;
        iEvent.getByToken( photonToken_, photons );

        //for ( auto photon = photons.product()->begin(); photon != photons.product()->end(); ++photon)
        for ( auto photon : *photons.product()) {
            // TODO: should we take the square root of the event weight for photons ?
            // TODO: check if this corresponds to the final event weight ?!
            float weight = photon.centralWeight();

            addPhoton(photon, weight);
        }
    } // analyze

    void
    TorchDumper::beginJob()
    {
    }

    void
    TorchDumper::endJob()
    {
        // barrel
        writeRecHits("/tmp/barrel-photons-rechits.t7", barrelRecHits, barrelWindowHalfWidth, barrelWindowHalfHeight);
        writeFloatVector("/tmp/barrel-photons-labels.t7", barrelLabels);
        writeFloatVector("/tmp/barrel-photons-weights.t7", barrelWeights);

        // endcap
        writeRecHits("/tmp/endcap-photons-rechits.t7", endcapRecHits, endcapWindowHalfWidth, endcapWindowHalfHeight);
        writeFloatVector("/tmp/endcap-photons-labels.t7", endcapLabels);
        writeFloatVector("/tmp/endcap-photons-weights.t7", endcapWeights);
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

} // namespace flashgg

typedef flashgg::TorchDumper FlashggTorchDumper;
DEFINE_FWK_MODULE( FlashggTorchDumper );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

