#include "flashgg/DataFormats/interface/TTHGenericTag.h"
#include <algorithm>

using namespace flashgg;

TTHGenericTag::TTHGenericTag() : DiPhotonTagBase::DiPhotonTagBase()
{}

TTHGenericTag::~TTHGenericTag()
{}

// N.B. Other attributes are set using methods in header file
TTHGenericTag::TTHGenericTag( edm::Ptr<DiPhotonCandidate> diPho, edm::Ptr<DiPhotonMVAResult> mvares ) : TTHGenericTag::TTHGenericTag( diPho, *mvares ) {}
TTHGenericTag::TTHGenericTag( edm::Ptr<DiPhotonCandidate> dipho, DiPhotonMVAResult mvares ) : DiPhotonTagBase::DiPhotonTagBase( dipho, mvares ) {}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

