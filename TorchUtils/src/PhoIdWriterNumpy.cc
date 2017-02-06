#include "flashgg/TorchUtils/interface/PhoIdWriterNumpy.h"
#include "flashgg/TorchUtils/plugins/PhoIdDumper.h"

#warning IMPLEMENT THIS
// #include "flashgg/TorchUtils/interface/TrackWriterTorch.h"

#include <fstream>

#include <boost/lexical_cast.hpp>

namespace flashgg 
{
    //----------------------------------------------------------------------

    void PhoIdWriterNumpy::writeNdArrayHeader(std::ostream &os, const std::string &dataType, 
                                              const std::vector<unsigned> &dimensions) {

        const string magic = "\x93NUMPY";
        const unsigned char verMajor = 0x01;
        const unsigned char verMinor = 0x00;
  
        os.write(magic.c_str(), magic.size());
        os.write((char *)&verMajor, 1);
        os.write((char *)&verMinor, 1);
    
        // description of the data
        string desc = "{'descr': '<" + dataType + "', 'fortran_order': False, 'shape': (";

        for (unsigned i = 0; i < dimensions.size(); ++i) {
            if (i > 0)
                desc += ", ";

            desc += boost::lexical_cast<string>(dimensions[i]);
        }

        desc += "), }";

        // padding
        while ((desc.size() + magic.size() + 4) % 16 != 15)
            desc += " ";
    
        desc += "\n";

        // write header length little endian
        unsigned short int headerLen = desc.size();
        os.write((char *)&headerLen, 2);
    
        // write header out
        os.write(desc.c_str(), desc.size());
    }

    //----------------------------------------------------------------------

    template<>
    void PhoIdWriterNumpy::writeTypeVector<float>(SimpleZipWriter &zip, const std::string &path, 
                                                  std::vector<float> &values) {
  
        ostringstream os;
        vector<unsigned> dimensions;
        dimensions.push_back(values.size());
        writeNdArrayHeader(os, "f4", dimensions);

        // first index row, second index is column

        for (unsigned row = 0; row < values.size(); ++row) {
            os.write((char*)&values[row], sizeof(float));
        }

        // add this to the zip archive
        zip.addFile(path, os.str());
    }

    //----------------------------------------------------------------------

    template<>
    void PhoIdWriterNumpy::writeTypeVector<int32_t>(SimpleZipWriter &zip, const std::string &path, 
                                                    std::vector<int32_t> &values) {
  
        ostringstream os;
        vector<unsigned> dimensions;
        dimensions.push_back(values.size());
        writeNdArrayHeader(os, "i4", dimensions);

        for (unsigned row = 0; row < values.size(); ++row) {
            os.write((char*)&values[row], sizeof(int32_t));
        }

        // add this to the zip archive
        zip.addFile(path, os.str());
    }
                                
    //----------------------------------------------------------------------

    void PhoIdWriterNumpy::writeRecHits(SimpleZipWriter &zip, const std::vector<std::vector<RecHitData> > &rechits, 
                                        int windowHalfWidth, int windowHalfHeight)
    {
        throw std::logic_error("not implemented");
    }
  
    //----------------------------------------------------------------------

    void PhoIdWriterNumpy::writeRecHitsSparse(SimpleZipWriter &zip, const std::vector<std::vector<RecHitData> > &rechits)
    {
        std::vector<int32_t> firstIndex(rechits.size()), numRecHits(rechits.size());
        unsigned nextStartIndex = 1;

        for (unsigned i = 0; i < rechits.size(); ++i)
            {
                firstIndex[i] = nextStartIndex;
                numRecHits[i] = rechits[i].size();
                
                nextStartIndex += rechits[i].size();

            } // loop over photons
            
        unsigned totNumRecHits = nextStartIndex - 1;
    

        writeTypeVector(zip, "X/firstIndex", firstIndex);
        writeTypeVector(zip, "X/numRecHits", numRecHits);

        // we add one to the x and y indices to stick to torch's one based indexing
        // (note that in the dense writing routine this is not needed because
        // we calculate the address directly)
        writeRecHitsValues<float>(zip, "X/energy", rechits, totNumRecHits, [](const RecHitData &rh)  { return rh.energy; });
        writeRecHitsValues<int32_t>(zip, "X/x", rechits, totNumRecHits, [](const RecHitData &rh)  { return rh.dx + 1; });
        writeRecHitsValues<int32_t>(zip, "X/y", rechits, totNumRecHits, [](const RecHitData &rh)  { return rh.dy + 1; });
    }

    //----------------------------------------------------------------------

    void PhoIdWriterNumpy::writeTo(PhoIdDumper &dumper, const std::string &fname)
    {
        // write a dict/table with 
        //  X = rechits
        //  y = labels
        //  weight = weights
        //  mvaid (for comparison)
        //  charged isolation w.r.t. chosen vertex
        //  charged isolation w.r.t. worst vertex

        std::ofstream os(fname.c_str());
        SimpleZipWriter zip(os);

        if (dumper.writeRecHitsSparseFlag)
            writeRecHitsSparse(zip, dumper.recHits);
        else
            writeRecHits(zip, dumper.recHits, dumper.windowHalfWidth, dumper.windowHalfHeight);

    
        writeTypeVector(zip, "y",      dumper.labels);
        writeTypeVector(zip, "weight", dumper.weights);
        writeTypeVector(zip, "mvaid",  dumper.mvaID);
        writeTypeVector(zip, "genDR",  dumper.genDeltaR);

        writeTypeVector(zip, "chgIsoWrtChosenVtx", dumper.chgIsoWrtChosenVtx);
        writeTypeVector(zip, "chgIsoWrtWorstVtx",  dumper.chgIsoWrtWorstVtx);

        if (dumper.writePhotonIdInputVarsFlag)
            {
                writeTypeVector(zip, "phoIdInput/scRawE"         , dumper.scRawE          );
                writeTypeVector(zip, "phoIdInput/r9"             , dumper.r9              );
                writeTypeVector(zip, "phoIdInput/covIEtaIEta"    , dumper.covIEtaIEta     );
                writeTypeVector(zip, "phoIdInput/phiWidth"       , dumper.phiWidth        );
                writeTypeVector(zip, "phoIdInput/etaWidth"       , dumper.etaWidth        );
                writeTypeVector(zip, "phoIdInput/covIEtaIPhi"    , dumper.covIEtaIPhi     );
                writeTypeVector(zip, "phoIdInput/s4"             , dumper.s4              );
                writeTypeVector(zip, "phoIdInput/pfPhoIso03"     , dumper.pfPhoIso03      );
                writeTypeVector(zip, "phoIdInput/pfChgIso03"     , dumper.pfChgIso03      );
                writeTypeVector(zip, "phoIdInput/pfChgIso03worst", dumper.pfChgIso03worst );
                writeTypeVector(zip, "phoIdInput/scEta"          , dumper.scEta           );
                writeTypeVector(zip, "phoIdInput/rho"            , dumper.rho             );
                writeTypeVector(zip, "phoIdInput/esEffSigmaRR"   , dumper.esEffSigmaRR    );
            }

        //----------
        // other photon variables (e.g. Et for Et flattening)
        //----------
        {
            writeTypeVector(zip, "phoVars/phoEt", dumper.photonEt);
        }

        //----------


#warning IMPLEMENT HTIS
        // tw.writeInt(tw.MAGIC_STRING); tw.writeString("tracks");      
        //    dynamic_cast<TrackWriterNumpy*>(dumper.trackWriter)->writeOut(tw);
    }

    //----------------------------------------------------------------------

}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

