#include "flashgg/DataFormats/interface/THQLeptonicTagTruth.h"
#include <iostream>

using namespace flashgg;

THQLeptonicTagTruth::THQLeptonicTagTruth() {}

THQLeptonicTagTruth::~THQLeptonicTagTruth() {}

/*
THQLeptonicTagTruth::THQLeptonicTagTruth(const THQLeptonicTagTruth &b) : TagTruthBase::TagTruthBase(b)
{
    std::cout << " Derived copy constructor!" << std::endl;
    setClosestGenJetToLeadingJet(b.closestGenJetToLeadingJet());
    setClosestGenJetToSubLeadingJet(b.closestGenJetToSubLeadingJet());
    setClosestParticleToLeadingJet(b.closestParticleToLeadingJet());
    setClosestParticleToSubLeadingJet(b.closestParticleToSubLeadingJet());
    setClosestParticleToLeadingPhoton(b.closestParticleToLeadingPhoton());
    setClosestParticleToSubLeadingPhoton(b.closestParticleToSubLeadingPhoton());
    setLeadingQuark(b.leadingParton());
    setSubLeadingQuark(b.subLeadingParton());
}
*/

THQLeptonicTagTruth *THQLeptonicTagTruth::clone() const
{
    //    return (new THQLeptonicTagTruth(*this));
    THQLeptonicTagTruth *result = new THQLeptonicTagTruth;
    result->setClosestGenJetToLeadingJet( closestGenJetToLeadingJet() );
    result->setClosestGenJetToSubLeadingJet( closestGenJetToSubLeadingJet() );
    result->setClosestGenJetToSubSubLeadingJet( closestGenJetToSubSubLeadingJet() );
    result->setClosestParticleToLeadingJet( closestParticleToLeadingJet() );
    result->setClosestParticleToSubLeadingJet( closestParticleToSubLeadingJet() );
    result->setClosestParticleToSubSubLeadingJet( closestParticleToSubSubLeadingJet() );
    result->setClosestPartonToLeadingJet( closestPartonToLeadingJet() );
    result->setClosestPartonToSubLeadingJet( closestPartonToSubLeadingJet() );
    result->setClosestPartonToSubSubLeadingJet( closestPartonToSubSubLeadingJet() );
    result->setClosestParticleToLeadingPhoton( closestParticleToLeadingPhoton() );
    result->setClosestParticleToSubLeadingPhoton( closestParticleToSubLeadingPhoton() );
    result->setClosestParticleToLeadingMuon( closestParticleToLeadingMuon() );
    result->setClosestParticleToSubLeadingMuon( closestParticleToSubLeadingMuon() );
    result->setClosestParticleToLeadingElectron( closestParticleToLeadingElectron() );
    result->setClosestParticleToSubLeadingElectron( closestParticleToSubLeadingElectron() );
    result->setClosestPromptParticleToLeadingMuon( closestPromptParticleToLeadingMuon() );
    result->setClosestPromptParticleToSubLeadingMuon( closestPromptParticleToSubLeadingMuon() );
    result->setClosestPromptParticleToLeadingElectron( closestPromptParticleToLeadingElectron() );
    result->setClosestPromptParticleToSubLeadingElectron( closestPromptParticleToSubLeadingElectron() );
    result->setLeadingJet( leadingJet() );
    result->setSubLeadingJet( subLeadingJet() );
    result->setSubSubLeadingJet( subSubLeadingJet() );
    result->setLeadingGenJet( leadingGenJet() );
    result->setSubLeadingGenJet( subLeadingGenJet() );
    result->setSubSubLeadingGenJet( subSubLeadingGenJet() );
    result->setLeadingParton( leadingParton() );
    result->setSubLeadingParton( subLeadingParton() );
    result->setSubSubLeadingParton( subSubLeadingParton() );
    result->setPtOrderedPartons( ptOrderedPartons() );
    result->setPtOrderedGenJets( ptOrderedGenJets() );
    result->setPtOrderedFggJets( ptOrderedFggJets() );
    result->setDiPhoton( diPhoton() );
    result->copyBaseInfo( *this );
    return result;

}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
