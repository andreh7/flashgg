#include "flashgg/DataFormats/interface/ZHLeptonicTag.h"
#include "flashgg/DataFormats/interface/WHLeptonicTag.h"
#include <algorithm>

using namespace flashgg;

ZHLeptonicTag::ZHLeptonicTag() : DiPhotonTagBase::DiPhotonTagBase()
{}

ZHLeptonicTag::~ZHLeptonicTag()
{}


ZHLeptonicTag::ZHLeptonicTag( edm::Ptr<DiPhotonCandidate> diPho, edm::Ptr<DiPhotonMVAResult> mvares ) : DiPhotonTagBase::DiPhotonTagBase( diPho, *mvares ) {}
ZHLeptonicTag::ZHLeptonicTag( edm::Ptr<DiPhotonCandidate> dipho, DiPhotonMVAResult mvares ) : DiPhotonTagBase::DiPhotonTagBase( dipho, mvares ) {}

double ZHLeptonicTag::transverseMass() const { return WHLeptonicTag::transverseMass(met(), electrons(), muons()); }


// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

