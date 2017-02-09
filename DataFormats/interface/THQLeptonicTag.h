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

        const std::vector<edm::Ptr<reco::Vertex> > vertices() const { return vertices_;}
        const std::vector<edm::Ptr<Muon> > muons() const { return Muons_;}
        const std::vector<edm::Ptr<flashgg::Electron> > electrons() const {return Electrons_;}
        const std::vector<edm::Ptr<Jet> > jets() const { return Jets_;}
        const std::vector<edm::Ptr<Jet> > Jets_EtaSorted() const { return Jets_EtaSorted_;}
        const std::vector<edm::Ptr<Jet> > bJets() const { return BJets_;}

        
        float getElecAlpha(int eleIndex) const{
            float eleta =  electrons()[eleIndex]->eta();

            //for isolation recalculation        
            float Aeff = 0;
            //cmssw/RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt
            if( fabs(eleta) <= 1.0000 ){
                Aeff = 0.1703;
            } else if( fabs(eleta) <= 1.4790 ){
                Aeff = 0.1715;
            } else if( fabs(eleta) <= 2.0000 ){
                Aeff = 0.1213;
            } else if( fabs(eleta) <= 2.2000 ){
                Aeff = 0.1230;
            } else if( fabs(eleta) <= 2.3000 ){
                Aeff = 0.1635;
            } else if( fabs(eleta) <= 2.4000 ){
                Aeff = 0.1937;
            } else if( fabs(eleta) <= 5.0000 ){
                Aeff = 0.2393;
            }

            return Aeff;
            
        };
        
        float getrho() const{
            return rho_;
        }
        /*
        float getMuoDz(int muIndex) const{
            mouons()[muIndex]->muonBestTrack()->dz( diPhoton()->vtx()->position() ) ;
        }
        */
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

        int findIndex(string label) const{
            ptrdiff_t pos = find(bAssignmentLabels.begin(), bAssignmentLabels.end(), label) - bAssignmentLabels.begin();
            if( pos < int( bAssignmentLabels.size() ) )
                return pos;
            else{
                //cout << "label : " << label < " not found" << endl;
                return -1;
            }
        }
        float thqleptonicMvaRes(string label) const {
            int index = findIndex(label);
            if(index < 0 )
                return -100;
            return thqleptonicMvaRes_[index];
        }
        const Ptr<Jet> getFwdJet(string label) const{
            int index = findIndex(label);
            if(index < 0 ){
                edm::Ptr<Jet> ret;
                return ret;
            }
            return fwdJet[index];
        }
        const Ptr<Jet> getbJet(string label) const{
            int index = findIndex(label);
            if(index < 0 ){
                edm::Ptr<Jet> ret;
                return ret;
            }

            return bJet[index];
        }
        float getTopMass(string label) const{
            int index = findIndex(label);
            if(index < 0 )
                return -100;

            return TopMass[findIndex(label)];
        }

        float getFoxWolframMoment_ONE()const{
            return FoxWolframMoment_ONE;
        }
        float getAplanarity() const{
            return Aplanarity;
        }
        float getMET() const{
            return MET;
        }
        float getMET_Phi() const{
            return MET_Phi;
        }

        const Ptr<reco::Vertex> getVertex( int vtx_index ) {
            if(vtx_index < 0 ){
                edm::Ptr<reco::Vertex> ret;
                return ret;
            }
            return vertices_[vtx_index];
        }

        int getLeadingMuonVertexDxy( ) {
            if (vtx_dxy_.size()>0)
                return vtx_dxy_[0];
            else
                return -999;
        }
        int getLeadingMuonVertexDz( ) {
            if (vtx_dz_.size()>0)
                return vtx_dz_[0];
            else
                return -999;
        }
        int getSubleadingMuonVertexDxy( ) {
            if (vtx_dxy_.size()>1)
                return vtx_dxy_[1];
            else
                return -999;
        }
        int getSubLeadingMuonVertexDz( ) {
            if (vtx_dz_.size()>1)
                return vtx_dz_[1];
            else
                return -999;
        }
        int getLeadingMuonDiphoVertexDxy( ) {
            if (diphovtx_dxy_.size()>0)
                return diphovtx_dxy_[0];
            else
                return -999;
        }
        int getLeadingMuonDiphoVertexDz( ) {
            if (diphovtx_dz_.size()>0)
                return diphovtx_dz_[0];
            else
                return -999;
        }
        int getSubleadingMuonDiphoVertexDxy( ) {
            if (diphovtx_dxy_.size()>1)
                return diphovtx_dxy_[1];
            else
                return -999;
        }
        int getSubLeadingMuonDiphoVertexDz( ) {
            if (diphovtx_dz_.size()>1)
                return diphovtx_dz_[1];
            else
                return -999;
        }

        void setJets( std::vector<edm::Ptr<Jet> > Jets, std::vector<edm::Ptr<Jet> > Jets_Eta ) {
            Jets_ = Jets; 
            Jets_EtaSorted_ = Jets_Eta;
        }
        void setBJets( std::vector<edm::Ptr<Jet> > BJets )  { BJets_ = BJets;}
        void setVertices( std::vector<edm::Ptr<reco::Vertex> > vertices ) {
            vertices_ = vertices;
        }

        void setLeadingMuonVertices( float a, float b, float c, float d) {
            vtx_dxy_.push_back(a);
            vtx_dz_.push_back(b);
            diphovtx_dxy_.push_back(c);
            diphovtx_dz_.push_back(d);
        }
        void setSubleadingMuonVertices( float a, float b, float c, float d) {
            vtx_dxy_.push_back(a);
            vtx_dz_.push_back(b);
            diphovtx_dxy_.push_back(c);
            diphovtx_dz_.push_back(d);
        }


        void setMuons( std::vector<edm::Ptr<Muon> > Muons ) {
            Muons_ = Muons;
        }
        void setElectrons( std::vector<edm::Ptr<Electron> > Electrons ) {
            Electrons_ = Electrons;
        }
        void setMVAres(string label, float val , float topMass , Ptr<Jet> fwd , Ptr<Jet> bj ) {
            bAssignmentLabels.push_back( label);
            thqleptonicMvaRes_.push_back(val); 
            TopMass.push_back(topMass);
            fwdJet.push_back( fwd ) ;
            bJet.push_back( bj ) ;
        }


        void setValues( float fox , float aplan ,  float met, float metPhi ){
            FoxWolframMoment_ONE = fox ;
            Aplanarity = aplan;

            MET = met ;
            MET_Phi = metPhi ;
        }
        //void setTHQLeptonicMVA( THQLeptonicMVAResult THQLeptonicMVA ) {THQLeptonicMVA_ = THQLeptonicMVA;}
        void setLeptonType(int val){ LeptonType_ = val; }

        void setrho(float rho){
            rho_=rho;
        }
    private:
        std::vector<edm::Ptr<Muon> > Muons_;
        std::vector<edm::Ptr<Electron> > Electrons_;
        std::vector<edm::Ptr<Jet> > Jets_;
        std::vector<edm::Ptr<Jet> > Jets_EtaSorted_;
        std::vector<edm::Ptr<Jet> > BJets_;
        std::vector<float> thqleptonicMvaRes_;
        std::vector<edm::Ptr<reco::Vertex> > vertices_;
        std::vector<float> vtx_dxy_; std::vector<float> vtx_dz_;
        std::vector<float> diphovtx_dxy_; std::vector<float> diphovtx_dz_;
        //THQLeptonicMVAResult THQLeptonicMVA_;
        float rho_;
        int LeptonType_;
        std::vector< std::string > bAssignmentLabels;
        std::vector< Ptr<Jet> > fwdJet ;
        std::vector< Ptr<Jet> > bJet;
        float FoxWolframMoment_ONE;
        float Aplanarity;
        std::vector<float> TopMass;
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

