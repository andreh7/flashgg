#include "flashgg/TorchUtils/interface/TrackWriterTorch.h"

namespace flashgg 
{
    TrackWriterTorch::TrackWriterTorch(edm::ParameterSet const& iPS_, edm::ConsumesCollector && iC_) :
        TrackWriter(iPS_, std::move(iC_)) 
    {
    }

    void TrackWriterTorch::writeOut(TorchWriter &tw)
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

        // number of tensors to write out
        const unsigned tableSize = 10;

        tw.writeInt(tw.MAGIC_TABLE);
        tw.writeInt(tw.getNextObjectIndex());
        tw.writeInt(tableSize); // number of variables

        tw.writeInt(tw.MAGIC_STRING); tw.writeString("firstIndex");  tw.writeTypeVector(firstIndex);
        tw.writeInt(tw.MAGIC_STRING); tw.writeString("numTracks");  tw.writeTypeVector(numTracks);


        tw.writeInt(tw.MAGIC_STRING); tw.writeString("pt"          ); writeFlattenedVector(tw,trackpt, totNumTracks);
        tw.writeInt(tw.MAGIC_STRING); tw.writeString("detaAtVertex"); writeFlattenedVector(tw,detaAtVertex, totNumTracks);
        tw.writeInt(tw.MAGIC_STRING); tw.writeString("dphiAtVertex"); writeFlattenedVector(tw,dphiAtVertex, totNumTracks);
        tw.writeInt(tw.MAGIC_STRING); tw.writeString("charge"      ); writeFlattenedVector(tw,charge, totNumTracks);
        tw.writeInt(tw.MAGIC_STRING); tw.writeString("pdgId"      ); writeFlattenedVector(tw,pdgId, totNumTracks);
        tw.writeInt(tw.MAGIC_STRING); tw.writeString("vtxZ"       );  writeFlattenedVector(tw,vtxX,  totNumTracks);
        tw.writeInt(tw.MAGIC_STRING); tw.writeString("vtxY"       ); writeFlattenedVector(tw,vtxY,  totNumTracks);
        tw.writeInt(tw.MAGIC_STRING); tw.writeString("vtxZ"       ); writeFlattenedVector(tw,vtxZ,  totNumTracks);
    }

    //----------------------------------------------------------------------

    // instantiate templated function in order to avoid missing
    // symbol at linking time
    void writeFlattenedVector(TorchWriter &tw, 
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
