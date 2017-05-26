#ifndef TorchUtils_TrackWriterNumpy_h
#define TorchUtils_TrackWriterNumpy_h

#include "flashgg/TorchUtils/interface/TrackWriter.h"
#include "flashgg/TorchUtils/interface/PhoIdWriterNumpy.h"

#include "flashgg/TorchUtils/interface/SimpleZipWriter.h"

namespace flashgg
{
    class TrackWriterNumpy : public TrackWriter 
    {
    public:
        TrackWriterNumpy(edm::ParameterSet const& iPS, edm::ConsumesCollector && iC);

        /** called at the end to write out the collected data */
        void writeOut(SimpleZipWriter &zip);

        //----------------------------------------

        template<typename DataType>
        void writeFlattenedVector(SimpleZipWriter &zip,
                                  const std::string &key,
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
            
            PhoIdWriterNumpy::writeTypeVector(zip, key, flatValues);
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
