#ifndef TorchUtils_PhoIdWriterTorch_h
#define TorchUtils_PhoIdWriterTorch_h

#include "flashgg/TorchUtils/interface/PhoIdWriter.h"
#include "flashgg/TorchUtils/interface/TorchWriter.h"

#include <string>
#include <cassert>
#include <iostream>



namespace flashgg
{
  class PhoIdDumper;

  /** interface to write out photon id variables to a file */
  class PhoIdWriterTorch : public PhoIdWriter 
  {
  public:
    virtual void writeTo(PhoIdDumper &dumper, const std::string &fname) override;

  protected:
    void writeRecHits(TorchWriter &tw, const std::vector<std::vector<RecHitData> > &rechits, 
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
    void writeRecHitsSparse(TorchWriter &tw, const std::vector<std::vector<RecHitData> > &rechits);

    //----------------------------------------

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
	
	std::cout << "nextPos=" << nextPos << " totNumRecHits=" << totNumRecHits << std::endl;
	assert(nextPos == totNumRecHits);
	
	tw.writeTypeVector<DataType>(values);
      }
    //----------------------------------------
  };
}

#endif
