#include "flashgg/TorchUtils/interface/PhoIdWriterTorch.h"
#include "flashgg/TorchUtils/plugins/PhoIdDumper.h"

#include "flashgg/TorchUtils/interface/TrackWriterTorch.h"

#include <fstream>

namespace flashgg 
{
  void PhoIdWriterTorch::writeRecHits(TorchWriter &tw, const std::vector<std::vector<RecHitData> > &rechits, 
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
  
  //----------------------------------------------------------------------

                                
  //----------------------------------------------------------------------

  void PhoIdWriterTorch::writeRecHitsSparse(TorchWriter &tw, const std::vector<std::vector<RecHitData> > &rechits)
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

  //----------------------------------------------------------------------

  /** converts vector of unsigned to a vector of signed integers
      to write it out to Torch files which do not support unsigned types */
  std::vector<int> PhoIdWriterTorch::toSigned(const std::vector<unsigned> &values) {
    std::vector<int> retval(values.size());
    for (unsigned i = 0; i < values.size(); ++i) {
      retval[i] = (int)values[i];
    }

    return retval;
  }

  //----------------------------------------------------------------------

  /** converts vector of unsigned long long to a vector of signed long long
      to write it out to Torch files which do not support unsigned types */
  std::vector<long long> PhoIdWriterTorch::toSigned(const std::vector<unsigned long long> &values) {
    std::vector<long long> retval(values.size());
    for (unsigned i = 0; i < values.size(); ++i) {
      retval[i] = (long long)values[i];
    }

    return retval;
  }

  //----------------------------------------------------------------------

  void PhoIdWriterTorch::writeTo(PhoIdDumper &dumper, const std::string &fname)
  {
    // write a dict/table with 
    //  X = rechits
    //  y = labels
    //  weight = weights
    //  mvaid (for comparison)
    //  charged isolation w.r.t. chosen vertex
    //  charged isolation w.r.t. worst vertex
    unsigned tableSize = 7;

    if (dumper.writePhotonIdInputVarsFlag)
      // add sub-recrod with photon id input variables
      tableSize += 1;
            
    // for writing out tracks
    tableSize += 1;

    // for writing other photon variables
    tableSize += 1;

    // for writing run/ls/event
    tableSize += 3;

    //----------

    std::ofstream os(fname.c_str());
    TorchWriter tw(os);

    tw.writeInt(tw.MAGIC_TABLE);

    tw.writeInt(tw.getNextObjectIndex());
            
    tw.writeInt(tableSize);

    tw.writeInt(tw.MAGIC_STRING); tw.writeString("X");
    if (dumper.writeRecHitsSparseFlag)
      writeRecHitsSparse(tw, dumper.recHits);
    else
      writeRecHits(tw, dumper.recHits, dumper.windowHalfWidth, dumper.windowHalfHeight);

    // event identification
    tw.writeInt(tw.MAGIC_STRING); tw.writeString("run");    tw.writeTypeVector(toSigned(dumper.runNumber));
    tw.writeInt(tw.MAGIC_STRING); tw.writeString("ls");     tw.writeTypeVector(toSigned(dumper.lsNumber));
    tw.writeInt(tw.MAGIC_STRING); tw.writeString("event");  tw.writeTypeVector(toSigned(dumper.eventNumber));

    tw.writeInt(tw.MAGIC_STRING); tw.writeString("y");      tw.writeTypeVector(dumper.labels);
    tw.writeInt(tw.MAGIC_STRING); tw.writeString("weight"); tw.writeTypeVector(dumper.weights);
    tw.writeInt(tw.MAGIC_STRING); tw.writeString("mvaid");  tw.writeTypeVector(dumper.mvaID);
    tw.writeInt(tw.MAGIC_STRING); tw.writeString("genDR");  tw.writeTypeVector(dumper.genDeltaR);

    tw.writeInt(tw.MAGIC_STRING); tw.writeString("chgIsoWrtChosenVtx");  tw.writeTypeVector(dumper.chgIsoWrtChosenVtx);
    tw.writeInt(tw.MAGIC_STRING); tw.writeString("chgIsoWrtWorstVtx");   tw.writeTypeVector(dumper.chgIsoWrtWorstVtx);

    if (dumper.writePhotonIdInputVarsFlag)
      {
	tw.writeInt(tw.MAGIC_STRING); tw.writeString("phoIdInput");      

	tw.writeInt(tw.MAGIC_TABLE);
	tw.writeInt(tw.getNextObjectIndex());
	tw.writeInt(13); // number of input variables

	tw.writeInt(tw.MAGIC_STRING); tw.writeString("scRawE"          );  tw.writeTypeVector(dumper.scRawE          );
	tw.writeInt(tw.MAGIC_STRING); tw.writeString("r9"              );  tw.writeTypeVector(dumper.r9              );
	tw.writeInt(tw.MAGIC_STRING); tw.writeString("covIEtaIEta"     );  tw.writeTypeVector(dumper.covIEtaIEta     );
	tw.writeInt(tw.MAGIC_STRING); tw.writeString("phiWidth"        );  tw.writeTypeVector(dumper.phiWidth        );
	tw.writeInt(tw.MAGIC_STRING); tw.writeString("etaWidth"        );  tw.writeTypeVector(dumper.etaWidth        );
	tw.writeInt(tw.MAGIC_STRING); tw.writeString("covIEtaIPhi"     );  tw.writeTypeVector(dumper.covIEtaIPhi     );
	tw.writeInt(tw.MAGIC_STRING); tw.writeString("s4"              );  tw.writeTypeVector(dumper.s4              );
	tw.writeInt(tw.MAGIC_STRING); tw.writeString("pfPhoIso03"      );  tw.writeTypeVector(dumper.pfPhoIso03      );
	tw.writeInt(tw.MAGIC_STRING); tw.writeString("pfChgIso03"      );  tw.writeTypeVector(dumper.pfChgIso03      );
	tw.writeInt(tw.MAGIC_STRING); tw.writeString("pfChgIso03worst" );  tw.writeTypeVector(dumper.pfChgIso03worst );
	tw.writeInt(tw.MAGIC_STRING); tw.writeString("scEta"           );  tw.writeTypeVector(dumper.scEta           );
	tw.writeInt(tw.MAGIC_STRING); tw.writeString("rho"             );  tw.writeTypeVector(dumper.rho             );
	tw.writeInt(tw.MAGIC_STRING); tw.writeString("esEffSigmaRR"    );  tw.writeTypeVector(dumper.esEffSigmaRR    );
      }

    //----------
    // other photon variables (e.g. Et for Et flattening)
    //----------
    {
      tw.writeInt(tw.MAGIC_STRING); tw.writeString("phoVars");
      tw.writeInt(tw.MAGIC_TABLE);
      tw.writeInt(tw.getNextObjectIndex());
      tw.writeInt(1); // number of variables

      tw.writeInt(tw.MAGIC_STRING); tw.writeString("phoEt"          );  tw.writeTypeVector(dumper.photonEt);
    }

    //----------


    tw.writeInt(tw.MAGIC_STRING); tw.writeString("tracks");      
    dynamic_cast<TrackWriterTorch*>(dumper.trackWriter)->writeOut(tw);
  }

  //----------------------------------------------------------------------

}
