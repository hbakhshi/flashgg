#ifndef FLASHgg_THQLeptonicTag_h
#define FLASHgg_THQLeptonicTag_h

#include "flashgg/DataFormats/interface/DiPhotonTagBase.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Jet.h"
//#include "flashgg/DataFormats/interface/THQLeptonicMVAResult.h"

using namespace edm;

namespace flashgg {

    class THQLeptonicTag: public DiPhotonTagBase
    {
    public:
        THQLeptonicTag();
        THQLeptonicTag( edm::Ptr<DiPhotonCandidate>, edm::Ptr<DiPhotonMVAResult> );
        THQLeptonicTag( edm::Ptr<DiPhotonCandidate>, DiPhotonMVAResult );

        ~THQLeptonicTag();

        THQLeptonicTag *clone() const { return ( new THQLeptonicTag( *this ) ); }

        const std::vector<edm::Ptr<Muon> > muons() const { return Muons_;}
        const std::vector<edm::Ptr<flashgg::Electron> > electrons() const {return Electrons_;}
        const std::vector<edm::Ptr<Jet> > jets() const { return Jets_;}
        const std::vector<edm::Ptr<Jet> > bJets() const { return BJets_;}
        const std::vector<edm::Ptr<Jet> > nonbJets() const { return nonBJets_;}
        float thqleptonicMvaRes() const {return thqleptonicMvaRes_;}
        
        // const reco::LeafCandidate getLepton() const{
        //     if( electrons().size() == 1 && muons().size() == 0 ){
        //         return ( *(electrons()[0]) );
        //     }else if( electrons().size() == 0 && muons().size() == 1 ){
        //         return ( *(muons()[0]) );
        //     }
            
        //     reco::LeafCandidate ret;
        //     return ret;
        // }
        int getLeptonCharge() const{
            if( LeptonType_ == 1 )
                return electrons()[0]->charge();
            else if(LeptonType_ == 2 )
                return muons()[0]->charge();

            return 0;
        }
        reco::Candidate::LorentzVector getLeptonP4() const{
            if( LeptonType_ == 1 )
                return electrons()[0]->p4();
            else if(LeptonType_ == 2 )
                return muons()[0]->p4();

            reco::Candidate::LorentzVector ret;
            return ret;
        }
        

        int getLeptonType() const{
            return LeptonType_;
        }

        const edm::Ptr<Jet> get_bJet() const{
            if( bJets().size() > 0 )
                return bJets()[0];
            
            edm::Ptr<Jet> ret;
            return ret;
        }

        const Ptr<Jet> getFwdJet() const{
            if( fwdJet.size() > 0 )
                return fwdJet[0];

            edm::Ptr<Jet> ret;
            return ret;
        }
        
        const Ptr<Jet> getbJet() const{
            if( bJet.size() > 0 )
                return bJet[0];

            edm::Ptr<Jet> ret;
            return ret;
        }
        float getFoxWolframMoment_ONE()const{
            return FoxWolframMoment_ONE;
        }
        float getAplanarity() const{
            return Aplanarity;
        }
        float getTopMass() const{
            return TopMass;
        }
        float getMET() const{
            return MET;
        }
        float getMET_Phi() const{
            return MET_Phi;
        }

        void setJets( std::vector<edm::Ptr<Jet> > Jets ) { Jets_ = Jets; }
        void setBJets( std::vector<edm::Ptr<Jet> > BJets )  { BJets_ = BJets;}
        void setLightJets( std::vector<edm::Ptr<Jet> > Jets )  { nonBJets_ = Jets;}
        void setMuons( std::vector<edm::Ptr<Muon> > Muons ) {
            Muons_ = Muons;
            if( Muons.size() == 1 )
                includeWeights( *Muons_[0] );
        }
        void setElectrons( std::vector<edm::Ptr<Electron> > Electrons ) {
            Electrons_ = Electrons;
            if( Electrons_.size() == 1 )
                includeWeights( *Electrons_[0] );
        }
        void setMVAres(float val) {thqleptonicMvaRes_ = val;}
        void setFwdJet( Ptr<Jet> fwd ) { fwdJet.clear() ; fwdJet.push_back( fwd ) ;}
        void setbJet( Ptr<Jet> bj ) { bJet.clear() ; bJet.push_back( bj ) ;}


        void setValues( float fox , float aplan , float topMass , float met, float metPhi ){
            FoxWolframMoment_ONE = fox ;
            Aplanarity = aplan;

            TopMass = topMass;

            MET = met ;
            MET_Phi = metPhi ;
        }
        //void setTHQLeptonicMVA( THQLeptonicMVAResult THQLeptonicMVA ) {THQLeptonicMVA_ = THQLeptonicMVA;}
        void setLeptonType(int val){ LeptonType_ = val; }
    private:
        std::vector<edm::Ptr<Muon> > Muons_;
        std::vector<edm::Ptr<Electron> > Electrons_;
        std::vector<edm::Ptr<Jet> > Jets_;
        std::vector<edm::Ptr<Jet> > BJets_;
        std::vector<edm::Ptr<Jet> > nonBJets_;
        float thqleptonicMvaRes_;
        //THQLeptonicMVAResult THQLeptonicMVA_;
        int LeptonType_;
        std::vector< Ptr<Jet> > fwdJet ;
        std::vector< Ptr<Jet> > bJet;
        float FoxWolframMoment_ONE;
        float Aplanarity;
        float TopMass;
        float MET;
        float MET_Phi;
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

