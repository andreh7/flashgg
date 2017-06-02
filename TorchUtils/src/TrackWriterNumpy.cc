#include "flashgg/TorchUtils/interface/TrackWriterNumpy.h"

namespace flashgg 
{
    TrackWriterNumpy::TrackWriterNumpy(edm::ParameterSet const& iPS_, edm::ConsumesCollector && iC_) :
        TrackWriter(iPS_, std::move(iC_)) 
    {
    }

    void TrackWriterNumpy::writeOut(SimpleZipWriter &zip)
    {
        // fill tensors to be written out, similar to photon rechits sparse
        // format

        //----------
        // calculate first indices
        //----------
        vector<int32_t> firstIndex(trackpt.size()), numTracks(trackpt.size());

        // Torch's indexing is one based
        int32_t nextStartIndex = 1;

        for (unsigned i = 0; i < trackpt.size(); ++i)
            {
                firstIndex[i] = nextStartIndex;
                numTracks[i] = trackpt[i].size();
      
                nextStartIndex += trackpt[i].size();
            } // loop over photons

        //----------
    
        int32_t totNumTracks = nextStartIndex - 1;

        PhoIdWriterNumpy::writeTypeVector(zip, "tracks/firstIndex", firstIndex);
        PhoIdWriterNumpy::writeTypeVector(zip, "tracks/numTracks", numTracks);

        writeFlattenedVector(zip, "tracks/pt",          trackpt, totNumTracks);
        writeFlattenedVector(zip, "tracks/detaAtVertex",detaAtVertex, totNumTracks);
        writeFlattenedVector(zip, "tracks/dphiAtVertex",dphiAtVertex, totNumTracks);
        writeFlattenedVector(zip, "tracks/charge"      ,charge, totNumTracks);
        writeFlattenedVector(zip, "tracks/pdgId"       ,pdgId, totNumTracks);
        writeFlattenedVector(zip, "tracks/vtxX"        ,vtxX,  totNumTracks);
        writeFlattenedVector(zip, "tracks/vtxY"        ,vtxY,  totNumTracks);
        writeFlattenedVector(zip, "tracks/vtxZ"        ,vtxZ,  totNumTracks);
    }

    //----------------------------------------------------------------------

    // instantiate templated function in order to avoid missing
    // symbol at linking time
    void writeFlattenedVector(SimpleZipWriter &zip,
                              const std::string &key,
                              const std::vector<std::vector<float> > &values,
                              unsigned totNumItems
                              );

}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
