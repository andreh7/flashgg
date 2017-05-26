#ifndef TorchUtils_PhoIdWriterNumpy_h
#define TorchUtils_PhoIdWriterNumpy_h

#include "flashgg/TorchUtils/interface/PhoIdWriter.h"
#include "flashgg/TorchUtils/interface/SimpleZipWriter.h"


#include <string>
#include <cassert>
#include <iostream>


namespace flashgg
{
    class PhoIdDumper;

    /** writes out photon id variables to a file */
    class PhoIdWriterNumpy : public PhoIdWriter 
    {
    public:
        virtual void writeTo(PhoIdDumper &dumper, const std::string &fname) override;

    protected:
        void static writeNdArrayHeader(std::ostream &os, const std::string &dataType, 
                                       const std::vector<unsigned> &dimensions);

        void writeRecHits(SimpleZipWriter &zip, const std::vector<std::vector<RecHitData> > &rechits, 
                          int windowHalfWidth, int windowHalfHeight);


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
        void writeRecHitsSparse(SimpleZipWriter &zip, const std::vector<std::vector<RecHitData> > &rechits);

        //----------------------------------------

    public:
        /** to write out nd arrays into the .npz file */
        template<typename DataType>
        static void writeTypeVector(SimpleZipWriter &zip, const std::string &path, 
                             std::vector<DataType> &values);

    protected:
        /** function to write out one element of each rechit as a table of tables */
        template<typename DataType, typename Func>
        void writeRecHitsValues(SimpleZipWriter &zip, 
                                const std::string &path, 
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
	
            assert(nextPos == totNumRecHits);
	
            writeTypeVector<DataType>(zip, path, values);
        }
        //----------------------------------------
    };
}

#endif


// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
