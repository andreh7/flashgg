#ifndef TorchUtils_TrackWriterTorch_h
#define TorchUtils_TrackWriterTorch_h

#include "flashgg/TorchUtils/interface/TrackWriter.h"

namespace flashgg
{
    class TrackWriterTorch : public TrackWriter 
    {
    public:
        TrackWriterTorch(edm::ParameterSet const& iPS, edm::ConsumesCollector && iC);

        /** called at the end to write out the collected data */
        void writeOut(TorchWriter &tw);

        //----------------------------------------

        template<typename DataType>
        void writeFlattenedVector(TorchWriter &tw, 
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


        //----------------------------------------



    };


} // namespace flashgg

#endif



// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
