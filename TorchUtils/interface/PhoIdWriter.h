#ifndef TorchUtils_PhoIdWriter_h
#define TorchUtils_PhoIdWriter_h

#include "DataFormats/DetId/interface/DetId.h"

#include <string>

namespace flashgg
{
  class PhoIdDumper;

  /** interface to write out photon id variables to a file */
  class PhoIdWriter 
  {
  public:
    /** msut write the data to the given output file */
    virtual void writeTo(flashgg::PhoIdDumper &dumper, const std::string &fname) = 0;

    struct RecHitData
    {
      float energy;
      int dx; // eta in the barrel
      int dy; // phi in the barrel

      DetId detid;
    };

  };
}

#endif
