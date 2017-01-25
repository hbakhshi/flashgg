#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/THQLeptonicTag.h"
//#include "flashgg/DataFormats/interface/THQLeptonicMVAResult.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "DataFormats/PatCandidates/interface/MET.h"

#include "flashgg/DataFormats/interface/Met.h"

#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "flashgg/Taggers/interface/LeptonSelection.h"

#include "DataFormats/Math/interface/deltaR.h"

//#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "flashgg/DataFormats/interface/VBFTagTruth.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "flashgg/Taggers/interface/SemiLepTopQuark.h"
#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"
#include "flashgg/Taggers/interface/FoxWolfram.hpp"

#include <vector>
#include <algorithm>
#include <string>
#include <utility>
#include "TLorentzVector.h"
#include "TMath.h"
#include "TMVA/Reader.h"

// https://github.com/cms-analysis/flashgg/commit/f327ca16c29b4ced8eaf8c309cb9218fac265963 (fixing the tth taggers)
using namespace std;
using namespace edm;


namespace flashgg {
    class THQLeptonicTagProducer : public EDProducer
    {

    public:
        typedef math::XYZPoint Point;

        THQLeptonicTagProducer( const ParameterSet & );
    private:
        void produce( Event &, const EventSetup & ) override;

        std::vector<edm::EDGetTokenT<View<flashgg::Jet> > > tokenJets_;
        EDGetTokenT<View<DiPhotonCandidate> > diPhotonToken_;
        //EDGetTokenT<View<THQLeptonicMVAResult> > thqleptonicMvaResultToken_;
        //EDGetTokenT<View<Jet> > thejetToken_;
        std::vector<edm::InputTag> inputTagJets_;
        EDGetTokenT<View<Electron> > electronToken_;
        EDGetTokenT<View<flashgg::Muon> > muonToken_;
        EDGetTokenT<View<DiPhotonMVAResult> > mvaResultToken_;
        EDGetTokenT<View<Photon> > photonToken_;
        EDGetTokenT<View<reco::Vertex> > vertexToken_;
        //EDGetTokenT<View<pat::MET> > METToken_;
        EDGetTokenT<View<flashgg::Met> > METToken_;
        EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
        EDGetTokenT<View<reco::GenJet> > genJetToken_;
        EDGetTokenT<double> rhoTag_;
        string systLabel_;

        typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;

        //Thresholds
        double leptonPtThreshold_;
        double leptonEtaThreshold_;
        vector<double> electronEtaThresholds_;
        double leadPhoOverMassThreshold_;
        double subleadPhoOverMassThreshold_;
        double MVAThreshold_;
        double deltaRLepPhoThreshold_;
        double deltaRJetLepThreshold_;

        double deltaRJetLeadPhoThreshold_;
        double deltaRJetSubLeadPhoThreshold_;

        double jetsNumberThreshold_;
        double bjetsNumberThreshold_;
        double jetPtThreshold_;
        double jetEtaThreshold_;

        vector<double> bDiscriminator_;
        string bTag_;
        double muPFIsoSumRelThreshold_;
        double PhoMVAThreshold_;
        double DeltaRTrkElec_;
        double TransverseImpactParam_;
        double LongitudinalImpactParam_;

        double deltaRPhoElectronThreshold_;
        double Zmass_;
        double deltaMassElectronZThreshold_;

        vector<double> nonTrigMVAThresholds_;
        vector<double> nonTrigMVAEtaCuts_;
        double electronIsoThreshold_;
        double electronNumOfHitsThreshold_;
        vector<double>  electronEtaCuts_;

        bool hasGoodElec = false;  bool hasVetoElec = false;
        bool hasGoodMuons = false;

        unique_ptr<TMVA::Reader> thqLeptonicMva_;
        FileInPath thqLeptonicMVAweightfile_;
        string  MVAMethod_;
        float thqLeptonicMvaResult_value_;

        std::vector< TLorentzVector > particles_LorentzVector;
        std::vector< math::RhoEtaPhiVector > particles_RhoEtaPhiVector;
        
        struct particleinfo{
            float pt, eta, phi , other , w , another; //other : for photon id, for diphoton mass, for jets btagging vals
            unsigned short number;
            bool isSet;
            TLorentzVector lorentzVector_;
            std::map<std::string,float> info;
            particleinfo( double pt_=-999, double eta_=-999, double phi_=-999 , double other_= -999 , double W= 1.0 ){
                pt = pt_;
                eta = eta_;
                phi = phi_;
                other = other_;
                w = W;
                number = 255;
                isSet = false;
                lorentzVector_.SetPtEtaPhiM(pt,eta,phi,0.);
            };
            void set(double pt_=-999, double eta_=-999, double phi_=-999 , double other_= -999 , double W= 1.0 , double Another= -999 ){
                pt = pt_;
                eta = eta_;
                phi = phi_;
                other = other_;
                w = W;
                another = Another;
                isSet = true;
                lorentzVector_.SetPtEtaPhiM(pt,eta,phi,0.);
            };
            TLorentzVector LorentzVector(){
                return lorentzVector_;
            };
            void SetLorentzVector(TLorentzVector lorentzVector){
                lorentzVector_.SetPxPyPzE(lorentzVector.Px(),lorentzVector.Py(),lorentzVector.Pz(),lorentzVector.Energy());
            };
        };
        
        particleinfo G1 , G2 , DiG , lepton , bjet, jprime, eventshapes , eventshapesMet , met , THReco , Top ;
        particleinfo foxwolf1 , foxwolf2 , foxwolf1Met, foxwolf2Met ;

        TLorentzVector metL, bL,fwdJL;  //temp solution: make met, bjet & jprime global TLorentzVectors

        struct GreaterByPt
        {
        public:
            bool operator()( edm::Ptr<flashgg::Jet> lh, edm::Ptr<flashgg::Jet> rh ) const
            {
                return lh->pt() > rh->pt();
            };
        };
        
        struct GreaterByEta
        {
        public:
            bool operator()( edm::Ptr<flashgg::Jet> lh, edm::Ptr<flashgg::Jet> rh ) const
            {
                return fabs(lh->eta()) > fabs(rh->eta());
            };
        };
        
        struct GreaterByBTagging
        {
        public:
            GreaterByBTagging(std::string urName):
                urName(urName)
            {
            }

            bool operator()( edm::Ptr<flashgg::Jet> lh, edm::Ptr<flashgg::Jet> rh ) const
            {
                return lh->bDiscriminator( urName.data() ) > rh->bDiscriminator( urName.data() );
            };
        private:
            const std::string urName;
        };

        std::vector<float> jetsPt;
        std::vector<float> jetsEta;
        std::vector<float> jetsPhi;
        std::vector<float> jetsE;
        std::vector<int> jetsIndex;

        
         std::vector<edm::Ptr<flashgg::Jet> > SelJetVect; std::vector<edm::Ptr<flashgg::Jet> > SelJetVect_EtaSorted; std::vector<edm::Ptr<flashgg::Jet> > SelJetVect_PtSorted; std::vector<edm::Ptr<flashgg::Jet> > SelJetVect_BSorted;
        std::vector<edm::Ptr<flashgg::Jet> > LightJetVect; std::vector<edm::Ptr<flashgg::Jet> > LightJetVect_EtaSorted; std::vector<edm::Ptr<flashgg::Jet> > LightJetVect_PtSorted; std::vector<edm::Ptr<flashgg::Jet> > LightJetVect_BSorted;
        std::vector<edm::Ptr<flashgg::Jet> > LooseBJetVect; std::vector<edm::Ptr<flashgg::Jet> > LooseBJetVect_EtaSorted; std::vector<edm::Ptr<flashgg::Jet> > LooseBJetVect_PtSorted; std::vector<edm::Ptr<flashgg::Jet> > LooseBJetVect_BSorted;
        std::vector<edm::Ptr<flashgg::Jet> > MediumBJetVect; std::vector<edm::Ptr<flashgg::Jet> > MediumBJetVect_EtaSorted; std::vector<edm::Ptr<flashgg::Jet> > MediumBJetVect_PtSorted; std::vector<edm::Ptr<flashgg::Jet> > MediumBJetVect_BSorted;
        std::vector<edm::Ptr<flashgg::Jet> > TightBJetVect;  std::vector<edm::Ptr<flashgg::Jet> > TightBJetVect_EtaSorted; std::vector<edm::Ptr<flashgg::Jet> > TightBJetVect_PtSorted; std::vector<edm::Ptr<flashgg::Jet> > TightBJetVect_BSorted;
        float  n_jets = 0; float n_bjets=0; float  n_ljets = 0; int n_lbjets = 0; int n_mbjets = 0; int n_tbjets = 0;
        
        float jprime_eta,met_pt;

    };

    THQLeptonicTagProducer::THQLeptonicTagProducer( const ParameterSet &iConfig ) :
        diPhotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
        //thqleptonicMvaResultToken_( consumes<View<flashgg::THQLeptonicMVAResult> >( iConfig.getParameter<InputTag> ( "THQLeptonicMVAResultTag" ) ) ),
        //thejetToken_( consumes<View<flashgg::Jet> >( iConfig.getParameter<InputTag>( "JetTag" ) ) ),
        inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
        electronToken_( consumes<View<flashgg::Electron> >( iConfig.getParameter<InputTag>( "ElectronTag" ) ) ),
        muonToken_( consumes<View<flashgg::Muon> >( iConfig.getParameter<InputTag>( "MuonTag" ) ) ),
        mvaResultToken_( consumes<View<flashgg::DiPhotonMVAResult> >( iConfig.getParameter<InputTag> ( "MVAResultTag" ) ) ),
        vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
        //METToken_( consumes<View<pat::MET> >( iConfig.getParameter<InputTag> ( "METTag" ) ) ),
        METToken_( consumes<View<flashgg::Met> >( iConfig.getParameter<InputTag> ( "METTag" ) ) ),
        genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleTag" ) ) ),
        genJetToken_ ( consumes<View<reco::GenJet> >( iConfig.getParameter<InputTag> ( "GenJetTag" ) ) ),
        rhoTag_( consumes<double>( iConfig.getParameter<InputTag>( "rhoTag" ) ) ),
        systLabel_( iConfig.getParameter<string> ( "SystLabel" ) ),
        MVAMethod_    ( iConfig.getParameter<string> ( "MVAMethod"    ) )
    {

        double default_leptonPtThreshold_ = 20.;
        double default_leptonEtaThreshold_ = 2.4;
        double default_leadPhoOverMassThreshold_ = 0.5;
        double default_subleadPhoOverMassThreshold_ = 0.25;
        double default_MVAThreshold_ = -0.6;
        double default_deltaRLepPhoThreshold_ = 0.5;
        double default_deltaRJetLepThreshold_ = 0.5;
        double default_jetsNumberThreshold_ = 2;
        double default_bjetsNumberThreshold_ = 1;
        double default_jetPtThreshold_ = 30.;
        double default_jetEtaThreshold_ = 4.7;


        double default_deltaRJetLeadPhoThreshold_ = 0.5;
        double default_deltaRJetSubLeadPhoThreshold_ = 0.5;

        vector<double> default_bDiscriminator_;
        default_bDiscriminator_.push_back( 0.605 );
        default_bDiscriminator_.push_back( 0.890 );

        string default_bTag_ = "combinedInclusiveSecondaryVertexV2BJetTags";
        double default_muPFIsoSumRelThreshold_ = 0.25;
        double default_PhoMVAThreshold_ = -0.2;
        double default_DeltaRTrkElec_ = 1.;
        double default_TransverseImpactParam_ = 0.02;
        double default_LongitudinalImpactParam_ = 0.2;

        double default_deltaRPhoElectronThreshold_ = 1.;
        double default_Zmass_ = 91.9;
        double default_deltaMassElectronZThreshold_ = 10.;

        vector<double> default_nonTrigMVAThresholds_;
        default_nonTrigMVAThresholds_.push_back(0.913286);
        default_nonTrigMVAThresholds_.push_back(0.805013);
        default_nonTrigMVAThresholds_.push_back(0.358969);

        vector<double> default_nonTrigMVAEtaCuts_;
        default_nonTrigMVAEtaCuts_.push_back(0.8);
        default_nonTrigMVAEtaCuts_.push_back(1.479);
        default_nonTrigMVAEtaCuts_.push_back(2.5);

        double default_electronIsoThreshold_ = 0.15;
        double default_electronNumOfHitsThreshold_ = 1.;

        vector<double> default_electronEtaCuts_;
        default_electronEtaCuts_.push_back( 1.4442 );
        default_electronEtaCuts_.push_back( 1.566 );
        default_electronEtaCuts_.push_back( 2.5 );

        leptonPtThreshold_ = iConfig.getUntrackedParameter<double>( "leptonPtThreshold", default_leptonPtThreshold_ );
        leptonEtaThreshold_ = iConfig.getUntrackedParameter<double>( "leptonEtaThreshold", default_leptonEtaThreshold_ );
        electronEtaThresholds_ = iConfig.getParameter<vector<double > >( "electronEtaThresholds");
        leadPhoOverMassThreshold_ = iConfig.getUntrackedParameter<double>( "leadPhoOverMassThreshold", default_leadPhoOverMassThreshold_ );
        subleadPhoOverMassThreshold_ = iConfig.getUntrackedParameter<double>( "subleadPhoOverMassThreshold", default_subleadPhoOverMassThreshold_ );
        MVAThreshold_ = iConfig.getUntrackedParameter<double>( "MVAThreshold", default_MVAThreshold_ );
        deltaRLepPhoThreshold_ = iConfig.getUntrackedParameter<double>( "deltaRLepPhoThreshold", default_deltaRLepPhoThreshold_ );
        deltaRJetLepThreshold_ = iConfig.getUntrackedParameter<double>( "deltaRJetLepThreshold", default_deltaRJetLepThreshold_ );
        jetsNumberThreshold_ = iConfig.getUntrackedParameter<double>( "jetsNumberThreshold", default_jetsNumberThreshold_ );
        bjetsNumberThreshold_ = iConfig.getUntrackedParameter<double>( "bjetsNumberThreshold", default_bjetsNumberThreshold_ );
        jetPtThreshold_ = iConfig.getUntrackedParameter<double>( "jetPtThreshold", default_jetPtThreshold_ );
        jetEtaThreshold_ = iConfig.getUntrackedParameter<double>( "jetEtaThreshold", default_jetEtaThreshold_ );

        deltaRJetLeadPhoThreshold_ = iConfig.getUntrackedParameter<double>( "deltaRJetLeadPhoThreshold", default_deltaRJetLeadPhoThreshold_ );
        deltaRJetSubLeadPhoThreshold_ = iConfig.getUntrackedParameter<double>( "deltaRJetSubLeadPhoThreshold", default_deltaRJetSubLeadPhoThreshold_ );

        electronEtaCuts_ = iConfig.getUntrackedParameter<vector<double > >( "electronEtaCuts",default_electronEtaCuts_);
        bDiscriminator_ = iConfig.getUntrackedParameter<vector<double > >( "bDiscriminator", default_bDiscriminator_ );
        bTag_ = iConfig.getUntrackedParameter<string>( "bTag", default_bTag_ );

        muPFIsoSumRelThreshold_ = iConfig.getUntrackedParameter<double>( "muPFIsoSumRelThreshold", default_muPFIsoSumRelThreshold_ );
        PhoMVAThreshold_ = iConfig.getUntrackedParameter<double>( "PhoMVAThreshold", default_PhoMVAThreshold_ );
        DeltaRTrkElec_ = iConfig.getUntrackedParameter<double>( "DeltaRTrkElec", default_DeltaRTrkElec_ );
        TransverseImpactParam_ = iConfig.getUntrackedParameter<double>( "TransverseImpactParam", default_TransverseImpactParam_ );
        LongitudinalImpactParam_ = iConfig.getUntrackedParameter<double>( "LongitudinalImpactParam", default_LongitudinalImpactParam_ );

        deltaRPhoElectronThreshold_ = iConfig.getUntrackedParameter<double>( "deltaRPhoElectronThreshold", default_deltaRPhoElectronThreshold_ );
        Zmass_ = iConfig.getUntrackedParameter<double>( "Zmass_", default_Zmass_ );
        deltaMassElectronZThreshold_ = iConfig.getUntrackedParameter<double>( "deltaMassElectronZThreshold_", default_deltaMassElectronZThreshold_ );

        nonTrigMVAThresholds_ =  iConfig.getUntrackedParameter<vector<double > >( "nonTrigMVAThresholds", default_nonTrigMVAThresholds_ );
        nonTrigMVAEtaCuts_ =  iConfig.getUntrackedParameter<vector<double > >( "nonTrigMVAEtaCuts", default_nonTrigMVAEtaCuts_ );
        electronIsoThreshold_ = iConfig.getUntrackedParameter<double>( "electronIsoThreshold", default_electronIsoThreshold_ );
        electronNumOfHitsThreshold_ = iConfig.getUntrackedParameter<double>( "electronNumOfHitsThreshold", default_electronNumOfHitsThreshold_ );

     
        thqLeptonicMVAweightfile_ = iConfig.getParameter<edm::FileInPath>( "thqleptonicMVAweightfile" );

        if (MVAMethod_ != ""){
            thqLeptonicMva_.reset( new TMVA::Reader( "!Color:Silent" ) );

            thqLeptonicMva_->AddVariable( "nJets"              , &n_jets    );
            thqLeptonicMva_->AddVariable( "Max$(abs(jetsEta))" , &jprime_eta);  //jprime.eta 
            thqLeptonicMva_->AddVariable( "met.pt"             , &met_pt    );  //met.pt 
            thqLeptonicMva_->AddVariable( "lepton.charge"      , &lepton.another    );
            thqLeptonicMva_->AddVariable( "eventshapes.aplanarity", &eventshapes.pt);
            thqLeptonicMva_->AddVariable( "foxwolf1.ONE"       , &foxwolf1.another);
            
            thqLeptonicMva_->BookMVA( MVAMethod_.c_str() , thqLeptonicMVAweightfile_.fullPath() );
        }

        for (unsigned i = 0 ; i < inputTagJets_.size() ; i++) {
            auto token = consumes<View<flashgg::Jet> >(inputTagJets_[i]);
            tokenJets_.push_back(token);
        }
        produces<vector<THQLeptonicTag> >();
        produces<vector<VBFTagTruth> >();
    }

    void THQLeptonicTagProducer::produce( Event &evt, const EventSetup & )

    {

        //Handle<View<flashgg::Jet> > theJets;
        //evt.getByToken( thejetToken_, theJets );
        //const PtrVector<flashgg::Jet>& jetPointers = theJets->ptrVector();
        JetCollectionVector Jets( inputTagJets_.size() );
        for( size_t j = 0; j < inputTagJets_.size(); ++j ) {
            evt.getByToken( tokenJets_[j], Jets[j] );
        }

        edm::Handle<double>  rho;
        evt.getByToken(rhoTag_,rho);
        double rho_    = *rho;

        Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
        evt.getByToken( diPhotonToken_, diPhotons );

        //Handle<View<flashgg::THQLeptonicMVAResult> > thqleptonicMvaResults;
        //evt.getByToken( thqleptonicMvaResultToken_, thqleptonicMvaResults );

        Handle<View<flashgg::Muon> > theMuons;
        evt.getByToken( muonToken_, theMuons );

        Handle<View<flashgg::Electron> > theElectrons;
        evt.getByToken( electronToken_, theElectrons );

        Handle<View<flashgg::DiPhotonMVAResult> > mvaResults;
        evt.getByToken( mvaResultToken_, mvaResults );

        Handle<View<reco::GenParticle> > genParticles;
        Handle<View<reco::GenJet> > genJets;

        std::auto_ptr<vector<THQLeptonicTag> > thqltags( new vector<THQLeptonicTag> );

        Handle<View<reco::Vertex> > vertices;
        evt.getByToken( vertexToken_, vertices );

        //Handle<View<pat::MET> > METs;
        //evt.getByToken( METToken_, METs );

        Handle<View<flashgg::Met> > METs;
        evt.getByToken( METToken_, METs );

        /*
        std::auto_ptr<vector<VBFTagTruth> > truths( new vector<VBFTagTruth> );
        
        unsigned int index_leadq       = std::numeric_limits<unsigned int>::max();
        unsigned int index_subleadq    = std::numeric_limits<unsigned int>::max();
        unsigned int index_subsubleadq    = std::numeric_limits<unsigned int>::max();
        unsigned int index_leadjet       = std::numeric_limits<unsigned int>::max();
        unsigned int index_subleadjet    = std::numeric_limits<unsigned int>::max();
        float pt_leadq = 0., pt_subleadq = 0., pt_subsubleadq = 0.;
        float pt_leadjet = 0., pt_subleadjet = 0.; 
        Point higgsVtx;

        if( ! evt.isRealData() ) {
            evt.getByToken( genJetToken_, genJets );
            evt.getByToken( genParticleToken_, genParticles );
            for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
                int pdgid = genParticles->ptrAt( genLoop )->pdgId();
                if( pdgid == 25 || pdgid == 22 ) {
                    higgsVtx = genParticles->ptrAt( genLoop )->vertex();
                    break;
                }
            }
            for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
                edm::Ptr<reco::GenParticle> part = genParticles->ptrAt( genLoop );
                if( part->isHardProcess() ) {
                    if( abs( part->pdgId() ) <= 5 ) {
                        if( part->pt() > pt_leadq ) {
                            index_subleadq = index_leadq;
                            pt_subleadq = pt_leadq;
                            index_leadq = genLoop;
                            pt_leadq = part->pt();
                        } else if( part->pt() > pt_subleadq ) {
                            index_subsubleadq  = index_subleadq;
                            pt_subsubleadq     = pt_subleadq;
                            index_subleadq = genLoop;
                            pt_subleadq    = part->pt();
                        }else if( part->pt() > pt_subsubleadq ){
                            index_subsubleadq = genLoop;
                            pt_subleadq       = part->pt();
                        }
                    }
                }
            }
        }
        
        edm::RefProd<vector<VBFTagTruth> > rTagTruth = evt.getRefBeforePut<vector<VBFTagTruth> >();
        unsigned int idx = 0;
        */

        assert( diPhotons->size() == mvaResults->size() );

        bool photonSelection = false;
        double idmva1 = 0.;
        double idmva2 = 0.;


        
        for( unsigned int diphoIndex = 0; diphoIndex < diPhotons->size(); diphoIndex++ ) {

            hasGoodElec = false; hasVetoElec = false;
            hasGoodMuons = false;
            
            unsigned int jetCollectionIndex = diPhotons->ptrAt( diphoIndex )->jetCollectionIndex();
            
            std::vector<edm::Ptr<Muon> > tagMuons;
            std::vector<edm::Ptr<Electron> > tagElectrons;
            std::vector<edm::Ptr<Jet> > tagJets;
            std::vector<edm::Ptr<Jet> > tagBJets;
            std::vector<edm::Ptr<Jet> > nontagBJets;

            edm::Ptr<flashgg::DiPhotonCandidate> dipho = diPhotons->ptrAt( diphoIndex );
            edm::Ptr<flashgg::DiPhotonMVAResult> mvares = mvaResults->ptrAt( diphoIndex );

            //edm::Ptr<flashgg::THQLeptonicMVAResult> thqleptonic_mvares = thqleptonicMvaResults->ptrAt( diphoIndex );
            //THQLeptonicMVA_ = thqleptonicMvaResults->ptrAt( diphoIndex );
            //float thqleptonic_mva = thqleptonic_mvares->thqLeptonicMvaResult_value;

            //std::cout << "mva:" << thqleptonic_mva << std::endl;

            flashgg::THQLeptonicTag thqltags_obj( dipho, mvares );

            if( dipho->leadingPhoton()->pt() < ( dipho->mass() )*leadPhoOverMassThreshold_ ) { continue; }

            if( dipho->subLeadingPhoton()->pt() < ( dipho->mass() )*subleadPhoOverMassThreshold_ ) { continue; }


            idmva1 = dipho->leadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
            idmva2 = dipho->subLeadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );

            if( idmva1 <= PhoMVAThreshold_ || idmva2 <= PhoMVAThreshold_ ) { continue; }

            if( mvares->result < MVAThreshold_ ) { continue; }

            photonSelection = true;
            
        
            G1.set( diPhotons->ptrAt( diphoIndex )->leadingPhoton()->pt(),
                    diPhotons->ptrAt( diphoIndex )->leadingPhoton()->eta(),
                    diPhotons->ptrAt( diphoIndex )->leadingPhoton()->phi(),
                    diPhotons->ptrAt( diphoIndex )->leadingPhoton()->phoIdMvaDWrtVtx( diPhotons->ptrAt( diphoIndex )->vtx() ),
                    diPhotons->ptrAt( diphoIndex )->leadingPhoton()->centralWeight() );
            G1.number = diPhotons->ptrAt( diphoIndex )->leadingPhoton()->genMatchType();
            particles_LorentzVector.push_back(G1.LorentzVector());
            
            G2.set( diPhotons->ptrAt( diphoIndex )->subLeadingPhoton()->pt(),
                    diPhotons->ptrAt( diphoIndex )->subLeadingPhoton()->eta(),
                    diPhotons->ptrAt( diphoIndex )->subLeadingPhoton()->phi(),
                    diPhotons->ptrAt( diphoIndex )->subLeadingPhoton()->phoIdMvaDWrtVtx( diPhotons->ptrAt( diphoIndex )->vtx() ),
                    diPhotons->ptrAt( diphoIndex )->subLeadingPhoton()->centralWeight() );
            G2.number = diPhotons->ptrAt( diphoIndex )->subLeadingPhoton()->genMatchType() ;
            particles_LorentzVector.push_back(G2.LorentzVector());

            particles_RhoEtaPhiVector.push_back( math::RhoEtaPhiVector(G1.pt, G1.eta , G1.phi) );
            particles_RhoEtaPhiVector.push_back( math::RhoEtaPhiVector(G2.pt, G2.eta , G2.phi) );

            if( METs->size() != 1 ) { std::cout << "WARNING - #MET is not 1" << std::endl;}
            Ptr<flashgg::Met> theMET = METs->ptrAt( 0 );

            //const pat::MET &met_ = METs->front();
            //std::cout << met_.pt() <<std::endl;
            metL.SetPtEtaPhiE( theMET->getCorPt(),
                               theMET->eta(),
                               theMET->getCorPhi(),
                               theMET->energy()
                               ) ; 


            std::vector<edm::Ptr<flashgg::Muon> > goodMuons = selectMuons( theMuons->ptrs(), dipho, vertices->ptrs(), leptonEtaThreshold_ , leptonPtThreshold_,
                                                                           muPFIsoSumRelThreshold_, deltaRLepPhoThreshold_, deltaRLepPhoThreshold_ );


            //std::vector<edm::Ptr<Electron> > goodElectrons = selectElectrons( theElectrons->ptrs(), dipho, vertices->ptrs(), leptonPtThreshold_, 
            //                                                                  TransverseImpactParam_, LongitudinalImpactParam_, nonTrigMVAThresholds_, nonTrigMVAEtaCuts_,
            //                                                                  electronIsoThreshold_, electronNumOfHitsThreshold_, electronEtaCuts_ ,
            //                                                                 deltaRPhoElectronThreshold_,DeltaRTrkElec_,deltaMassElectronZThreshold_);
            
            std::vector<edm::Ptr<Electron> > vetoElectrons = selectStdElectrons(theElectrons->ptrs(), dipho, vertices->ptrs(), leptonPtThreshold_,  electronEtaThresholds_ ,
                                                                                0,0,
                                                                                deltaRPhoElectronThreshold_,DeltaRTrkElec_,deltaMassElectronZThreshold_ , rho_, evt.isRealData());
            

            std::vector<edm::Ptr<Electron> > goodElectrons = selectStdElectrons(theElectrons->ptrs(), dipho, vertices->ptrs(), leptonPtThreshold_,  electronEtaThresholds_ ,
                                                                                0,1,
                                                                                deltaRPhoElectronThreshold_,DeltaRTrkElec_,deltaMassElectronZThreshold_ , rho_, evt.isRealData());

            hasGoodElec = ( goodElectrons.size() == 1 ); hasVetoElec = ( vetoElectrons.size() > 0 );
            hasGoodMuons = ( goodMuons.size() == 1 );


            int LeptonType = 0; //1 : electron, 2:muon

            if( hasGoodMuons && !hasVetoElec) {

                LeptonType = 2;
                for( unsigned int muonIndex = 0; muonIndex < goodMuons.size(); muonIndex++ ) {

                    Ptr<flashgg::Muon> muon = goodMuons[muonIndex];

                    tagMuons.push_back( muon );

                    lepton.set( muon->pt(),
                                muon->eta() ,
                                muon->phi() ,
                                0. ,
                                1. ,
                                muon->charge() );
                    particles_LorentzVector.push_back(lepton.LorentzVector());
                    particles_RhoEtaPhiVector.push_back( math::RhoEtaPhiVector( lepton.pt, lepton.eta, lepton.phi ) );

                }//end of muons loop

            }

            if( hasGoodElec && !hasGoodMuons) {
                LeptonType = 1;

                for( unsigned int ElectronIndex = 0; ElectronIndex < goodElectrons.size(); ElectronIndex++ ) {

                    Ptr<Electron> Electron = goodElectrons[ElectronIndex];
                    tagElectrons.push_back( Electron );

                    lepton.set( Electron->pt(),
                                Electron->eta() ,
                                Electron->phi() ,
                                0. ,
                                1. ,
                                Electron->charge() );
                    particles_LorentzVector.push_back(lepton.LorentzVector());
                    particles_RhoEtaPhiVector.push_back( math::RhoEtaPhiVector( lepton.pt, lepton.eta, lepton.phi ) );


                }//end of electron loop

            }


            int deltaRLeptonJetcount = 0;

            for( unsigned int candIndex_outer = 0; candIndex_outer < Jets[jetCollectionIndex]->size() ; candIndex_outer++ ) {
                edm::Ptr<flashgg::Jet> thejet = Jets[jetCollectionIndex]->ptrAt( candIndex_outer );

                //std::cout << "prin: "<< Jets[jetCollectionIndex]->size() << " "<<thejet->pt() << " "<< thejet->eta()<<" "<< thejet->phi()<< " "<< thejet->energy() <<std::endl;

                if( !thejet->passesPuJetId( dipho ) ) { continue; }

                if( fabs( thejet->eta() ) > jetEtaThreshold_ ) { continue; }

                if( thejet->pt() < jetPtThreshold_ ) { continue; }

                float dRPhoLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->leadingPhoton()->superCluster()->eta(), dipho->leadingPhoton()->superCluster()->phi() ) ;
                float dRPhoSubLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->subLeadingPhoton()->superCluster()->eta(),
                                                dipho->subLeadingPhoton()->superCluster()->phi() );

                if( dRPhoLeadJet < deltaRJetLeadPhoThreshold_ || dRPhoSubLeadJet < deltaRJetSubLeadPhoThreshold_ ) { continue; }

                if( LeptonType != 0 ){
                    float dRJetLepton = deltaR( thejet->eta(), thejet->phi(), lepton.eta , lepton.phi ) ;

                    if( dRJetLepton < deltaRJetLepThreshold_ ) { continue; }
                }
                deltaRLeptonJetcount++;

                        

                tagJets.push_back( thejet );

                TLorentzVector jet_lorentzVector;
                jet_lorentzVector.SetPtEtaPhiE(  thejet->pt() , thejet->eta() , thejet->phi() , thejet->energy() );
                //std::cout <<  thejet->pt() << " "<< thejet->eta()<<" "<< thejet->phi()<< " "<< thejet->energy() <<std::endl;
                particles_LorentzVector.push_back( jet_lorentzVector );
                particles_RhoEtaPhiVector.push_back( math::RhoEtaPhiVector( thejet->pt(), thejet->eta(), thejet->phi() ) );
             

                double bDiscriminatorValue = thejet->bDiscriminator( bTag_.c_str() );

                if( bDiscriminatorValue > bDiscriminator_[1] ) {
                    tagBJets.push_back( thejet );
                    n_mbjets++;
                    MediumBJetVect.push_back( thejet ); 
                    MediumBJetVect_EtaSorted.push_back( thejet );
                    MediumBJetVect_PtSorted.push_back( thejet );
                    MediumBJetVect_BSorted.push_back( thejet );
                }
                else{
                    nontagBJets.push_back( thejet );
                    n_ljets++;
                    LightJetVect.push_back( thejet ); 
                    LightJetVect_EtaSorted.push_back( thejet );
                    LightJetVect_PtSorted.push_back( thejet );
                    LightJetVect_BSorted.push_back( thejet );
                }
                        
                SelJetVect.push_back( thejet ); 
                SelJetVect_EtaSorted.push_back( thejet );
                SelJetVect_PtSorted.push_back( thejet );
                SelJetVect_BSorted.push_back( thejet );
            }//end of jets loop
            std::sort(MediumBJetVect_EtaSorted.begin(),MediumBJetVect_EtaSorted.end(),GreaterByEta()); 
            std::sort(MediumBJetVect_PtSorted.begin(),MediumBJetVect_PtSorted.end(),GreaterByPt()); 
            std::sort(MediumBJetVect_BSorted.begin(),MediumBJetVect_BSorted.end(),GreaterByBTagging(bTag_.c_str())); 

            std::sort(LightJetVect_EtaSorted.begin(),LightJetVect_EtaSorted.end(),GreaterByEta()); 
            std::sort(LightJetVect_PtSorted.begin(),LightJetVect_PtSorted.end(),GreaterByPt()); 
            std::sort(LightJetVect_BSorted.begin(),LightJetVect_BSorted.end(),GreaterByBTagging(bTag_.c_str())); 

            std::sort(SelJetVect_EtaSorted.begin(),SelJetVect_EtaSorted.end(),GreaterByEta()); 
            std::sort(SelJetVect_PtSorted.begin(),SelJetVect_PtSorted.end(),GreaterByPt()); 
            std::sort(SelJetVect_BSorted.begin(),SelJetVect_BSorted.end(),GreaterByBTagging(bTag_.c_str())); 


            // std::cout << " THQLeptonicTagProducer tagBJets.size()=" << tagBJets.size()
            //           << " tagJets.size()=" << tagJets.size()
            //           << " photonSelection=" << photonSelection
            //           << " tagMuons.size()=" << tagMuons.size() << " muonJets=" << muonJets
            //           << " tagElectrons.size()="<< tagElectrons.size() << " ElectronJets=" << ElectronJets
            //           << std::endl;
            




            if( tagBJets.size() == bjetsNumberThreshold_ && tagJets.size() >= jetsNumberThreshold_ && photonSelection ){
                //&& ( ( (tagMuons.size() == 1 && muonJets) and  (tagElectrons.size() == 0 && !ElectronJets) )  || ( (tagMuons.size() == 0 && !muonJets)  and  (tagElectrons.size() == 1 && ElectronJets) ) ) ) 
                

                EventShapeVariables shapeVars(particles_RhoEtaPhiVector);
                //std::cout  << "aplanarity: "<<shapeVars.aplanarity()<<std::endl;
                eventshapes.set( shapeVars.aplanarity() ,
                                 shapeVars.C() ,
                                 shapeVars.circularity(),
                                 shapeVars.D() ,
                                 shapeVars.isotropy(),
                                 shapeVars.sphericity() );

                FoxWolfram fwam( particles_LorentzVector );
                std::vector< particleinfo*> allfoxwolfs = {&foxwolf1 , &foxwolf2 };
                for(uint ifw = 1 ; ifw < allfoxwolfs.size()+1 ; ifw++)
                    allfoxwolfs[ifw-1]->set( fwam.getMoment( FoxWolfram::SHAT , ifw ),
                                             fwam.getMoment( FoxWolfram::PT , ifw ),
                                             fwam.getMoment( FoxWolfram::ETA , ifw ),
                                             fwam.getMoment( FoxWolfram::PSUM , ifw ),
                                             fwam.getMoment( FoxWolfram::PZ , ifw ),
                                             fwam.getMoment( FoxWolfram::ONE , ifw ) );
                //std::cout<< "fox:" << foxwolf1.another<<std::endl;


                edm::Ptr<Jet> fwdJet = SelJetVect_EtaSorted[0];
                edm::Ptr<Jet> bJet = MediumBJetVect_BSorted[0];
                float topMass = -100.;
                if( deltaR( fwdJet->p4() , bJet->p4() ) < std::numeric_limits<double>::epsilon() )
                    fwdJet = SelJetVect_EtaSorted[1] ;
                if( LeptonType != 0){
                    
                    bL.SetPtEtaPhiE( bJet->pt(), bJet->eta(), bJet->phi(), bJet->energy());
                    fwdJL.SetPtEtaPhiE( fwdJet->pt(),fwdJet->eta(), fwdJet->phi(), fwdJet->energy());

                    
                    flashgg::SemiLepTopQuark singletop(bL, metL, lepton.LorentzVector(), fwdJL,fwdJL);
                    //met.SetLorentzVector(singletop.getMET());
                    n_jets = tagJets.size() ; 
                    metL = singletop.getMET() ;
                    jprime_eta  = fabs( fwdJL.Eta() );
                    met_pt = metL.Pt();

                    topMass = singletop.top().M() ;

                    if (MVAMethod_ != "") 
                        thqLeptonicMvaResult_value_ = thqLeptonicMva_->EvaluateMVA( MVAMethod_.c_str() );

                }            

                if( LeptonType == 1 )
                    thqltags_obj.includeWeights( *tagElectrons[0] );
                else if( LeptonType == 2 )
                    thqltags_obj.includeWeights( *tagMuons[0] );
                    

                thqltags_obj.setLeptonType(LeptonType);
                thqltags_obj.includeWeights( *dipho );
                thqltags_obj.setMVAres( thqLeptonicMvaResult_value_ );
                thqltags_obj.setFwdJet( fwdJet );
                thqltags_obj.setbJet( bJet );
                thqltags_obj.setJets( tagJets );
                thqltags_obj.setBJets( tagBJets );
                thqltags_obj.setLightJets( nontagBJets );
                thqltags_obj.setMuons( tagMuons );
                thqltags_obj.setElectrons( tagElectrons );
                thqltags_obj.setDiPhotonIndex( diphoIndex );
                thqltags_obj.setSystLabel( systLabel_ );
                thqltags_obj.setValues( foxwolf1.another , eventshapes.pt , topMass );

                thqltags->push_back( thqltags_obj );
            }//thq tag
            else {
                std::cout << " THQLeptonicTagProducer NO TAG " << std::endl;
            }

            n_jets = 0; n_bjets = 0; n_ljets = 0; n_lbjets = 0; n_mbjets = 0; n_tbjets = 0;

            jetsPhi.clear();
            jetsPt.clear();
            jetsEta.clear();
            jetsE.clear();
            jetsIndex.clear();
            particles_LorentzVector.clear();
            particles_RhoEtaPhiVector.clear();
            SelJetVect.clear(); SelJetVect_EtaSorted.clear(); SelJetVect_PtSorted.clear(); SelJetVect_BSorted.clear();
            LightJetVect.clear(); LightJetVect_EtaSorted.clear(); LightJetVect_PtSorted.clear(); LightJetVect_BSorted.clear();
            LooseBJetVect.clear(); LooseBJetVect_EtaSorted.clear(); LooseBJetVect_PtSorted.clear(); LooseBJetVect_BSorted.clear();
            MediumBJetVect.clear(); MediumBJetVect_EtaSorted.clear(); MediumBJetVect_PtSorted.clear(); MediumBJetVect_BSorted.clear();
            TightBJetVect.clear(); TightBJetVect_EtaSorted.clear(); TightBJetVect_PtSorted.clear(); TightBJetVect_BSorted.clear();

            
        }//diPho loop end !

        // if( thqltags->size() == 0 ){
        //     THQLeptonicTag empty;
        //     thqltags->push_back( empty );
            
        //     std::cout << " THQLeptonicTagProducer NO TAG is inserted in this event " << std::endl;
        // }
        
        evt.put( thqltags );
        //evt.put( truths );
            
        
    }
    
}
typedef flashgg::THQLeptonicTagProducer FlashggTHQLeptonicTagProducer;
DEFINE_FWK_MODULE( FlashggTHQLeptonicTagProducer );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

