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

#include "flashgg/DataFormats/interface/SemiLepTopQuark.h"
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
                return abs(lh->eta()) > abs(rh->eta());
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

        vector<int> numMuonJetsdR;
        vector<int> numElectronJetsdR;
        bool muonJets = false;
        bool ElectronJets = false;

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

            hasGoodElec = ( goodElectrons.size() > 0 ); hasVetoElec = ( vetoElectrons.size() > 0 );
            hasGoodMuons = ( goodMuons.size() > 0 );
            //if( !hasGoodElec && !hasGoodMuons ) { continue; }
            //if( hasGoodElec && hasGoodMuons ) { continue; }

            numMuonJetsdR.clear();
            numElectronJetsdR.clear();
            muonJets = false;
            ElectronJets = false;

            /*
            float pt1 = diPhotons->ptrAt( diphoIndex )->leadingPhoton()->phi();
            float phi1 = diPhotons->ptrAt( diphoIndex )->leadingPhoton()->phi();
            float eta1 = diPhotons->ptrAt( diphoIndex )->leadingPhoton()->eta();
            float pt2 = diPhotons->ptrAt( diphoIndex )->subLeadingPhoton()->phi();
            float phi2 = diPhotons->ptrAt( diphoIndex )->subLeadingPhoton()->phi();
            float eta2 = diPhotons->ptrAt( diphoIndex )->subLeadingPhoton()->eta();


            std::cout << "jet index!!!!!!!!!!!! "<< jetCollectionIndex<< std::endl;

            for( UInt_t jetLoop = 0; jetLoop < Jets[jetCollectionIndex]->size() ; jetLoop++ ) {
                Ptr<flashgg::Jet> jet  = Jets[jetCollectionIndex]->ptrAt( jetLoop );
                
                // within eta 4.7?
                if( fabs( jet->eta() ) > 4.7 ) { continue; }
                //pt less than 
                if ( jet->pt() <  30 ) { continue ; }
                if (!jet->passesJetID  ( flashgg::Loose ) ) continue;
                // close to lead photon?
                float dPhi = deltaPhi( jet->phi(), phi1 );
                float dEta = jet->eta() - eta1;
                if( sqrt( dPhi * dPhi + dEta * dEta ) < 0.4 ) { continue; }
                
                // close to sublead photon?
                dPhi = deltaPhi( jet->phi(), phi2 );
                dEta = jet->eta() - eta2;
                if( sqrt( dPhi * dPhi + dEta * dEta ) < 0.4 ) { continue; }

                tagJets.push_back( jet );

            }
            */


            if( hasGoodMuons && !hasVetoElec) {

                for( unsigned int muonIndex = 0; muonIndex < goodMuons.size(); muonIndex++ ) {

                    Ptr<flashgg::Muon> muon = goodMuons[muonIndex];

                    int deltaRMuonJetcount = 0;
                    double bDiscriminatorValue = -999.;


                    for( unsigned int candIndex_outer = 0; candIndex_outer < Jets[jetCollectionIndex]->size() ; candIndex_outer++ ) {
                        edm::Ptr<flashgg::Jet> thejet = Jets[jetCollectionIndex]->ptrAt( candIndex_outer );

                        std::cout << "prin: "<< Jets[jetCollectionIndex]->size() << " "<<thejet->pt() << " "<< thejet->eta()<<" "<< thejet->phi()<< " "<< thejet->energy() <<std::endl;

                        if( !thejet->passesPuJetId( dipho ) ) { continue; }

                        if( fabs( thejet->eta() ) > jetEtaThreshold_ ) { continue; }

                        if( thejet->pt() < jetPtThreshold_ ) { continue; }

                        float dRPhoLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->leadingPhoton()->superCluster()->eta(), dipho->leadingPhoton()->superCluster()->phi() ) ;
                        float dRPhoSubLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->subLeadingPhoton()->superCluster()->eta(),
                                                        dipho->subLeadingPhoton()->superCluster()->phi() );

                        if( dRPhoLeadJet < deltaRJetLeadPhoThreshold_ || dRPhoSubLeadJet < deltaRJetSubLeadPhoThreshold_ ) { continue; }

                        float dRJetMuon = deltaR( thejet->eta(), thejet->phi(), muon->eta(), muon->phi() ) ;

                        if( dRJetMuon < deltaRJetLepThreshold_ ) { continue; }
                        deltaRMuonJetcount++;

                        

                        tagJets.push_back( thejet );

                        TLorentzVector jet_lorentzVector;
                        jet_lorentzVector.SetPtEtaPhiE(  thejet->pt() , thejet->eta() , thejet->phi() , thejet->energy() );
                        std::cout <<  thejet->pt() << " "<< thejet->eta()<<" "<< thejet->phi()<< " "<< thejet->energy() <<std::endl;
                        particles_LorentzVector.push_back( jet_lorentzVector );
                        particles_RhoEtaPhiVector.push_back( math::RhoEtaPhiVector( thejet->pt(), thejet->eta(), thejet->phi() ) );
             

                        bDiscriminatorValue = thejet->bDiscriminator( bTag_.c_str() );

                        if( bDiscriminatorValue > bDiscriminator_[1] ) {
                            tagBJets.push_back( thejet );
                            n_mbjets++;
                            MediumBJetVect.push_back( thejet ); 
                            MediumBJetVect_EtaSorted.push_back( thejet );
                            MediumBJetVect_PtSorted.push_back( thejet );
                            MediumBJetVect_BSorted.push_back( thejet );
                            std::sort(MediumBJetVect_EtaSorted.begin(),MediumBJetVect_EtaSorted.end(),GreaterByEta()); 
                            std::sort(MediumBJetVect_PtSorted.begin(),MediumBJetVect_PtSorted.end(),GreaterByPt()); 
                            std::sort(MediumBJetVect_BSorted.begin(),MediumBJetVect_BSorted.end(),GreaterByBTagging(bTag_.c_str())); 
                        }
                        else{
                            nontagBJets.push_back( thejet );
                            n_ljets++;
                            LightJetVect.push_back( thejet ); 
                            LightJetVect_EtaSorted.push_back( thejet );
                            LightJetVect_PtSorted.push_back( thejet );
                            LightJetVect_BSorted.push_back( thejet );
                            std::sort(LightJetVect_EtaSorted.begin(),LightJetVect_EtaSorted.end(),GreaterByEta()); 
                            std::sort(LightJetVect_PtSorted.begin(),LightJetVect_PtSorted.end(),GreaterByPt()); 
                            std::sort(LightJetVect_BSorted.begin(),LightJetVect_BSorted.end(),GreaterByBTagging(bTag_.c_str())); 
                        }
                        
                        SelJetVect.push_back( thejet ); 
                        SelJetVect_EtaSorted.push_back( thejet );
                        SelJetVect_PtSorted.push_back( thejet );
                        SelJetVect_BSorted.push_back( thejet );
                        std::sort(SelJetVect_EtaSorted.begin(),SelJetVect_EtaSorted.end(),GreaterByEta()); 
                        std::sort(SelJetVect_PtSorted.begin(),SelJetVect_PtSorted.end(),GreaterByPt()); 
                        std::sort(SelJetVect_BSorted.begin(),SelJetVect_BSorted.end(),GreaterByBTagging(bTag_.c_str())); 

                    }//end of jets loop

                    numMuonJetsdR.push_back( deltaRMuonJetcount );
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

                std::vector<const flashgg::Photon *> photons;

                photons.push_back( dipho->leadingPhoton() );
                photons.push_back( dipho->subLeadingPhoton() );


                for( unsigned int ElectronIndex = 0; ElectronIndex < goodElectrons.size(); ElectronIndex++ ) {

                    Ptr<Electron> Electron = goodElectrons[ElectronIndex];

                    int deltaRElectronJetcount = 0;
                    double bDiscriminatorValue = -999.;

                    for( unsigned int candIndex_outer = 0; candIndex_outer < Jets[jetCollectionIndex]->size() ; candIndex_outer++ ) {
                        edm::Ptr<flashgg::Jet> thejet = Jets[jetCollectionIndex]->ptrAt( candIndex_outer );

                        if( !thejet->passesPuJetId( dipho ) ) { continue; }

                        //https://github.com/h2gglobe/h2gglobe/blob/master/PhotonAnalysis/src/PhotonAnalysis.cc#L5367
                        if( fabs( thejet->eta() ) > jetEtaThreshold_ ) { continue; }

                        //https://github.com/h2gglobe/h2gglobe/blob/master/PhotonAnalysis/src/PhotonAnalysis.cc#L5371
                        if( thejet->pt() < jetPtThreshold_ ) { continue; }

                        float dRJetElectron = deltaR( thejet->eta(), thejet->phi(), Electron->eta(), Electron->phi() ) ;

                        //https://github.com/njets_btagmediumh2gglobe/h2gglobe/blob/master/PhotonAnalysis/src/PhotonAnalysis.cc#L5370
                        if( dRJetElectron < deltaRJetLepThreshold_ ) { continue; }
                        deltaRElectronJetcount++;

                        tagJets.push_back( thejet );

                        TLorentzVector jet_lorentzVector;
                        jet_lorentzVector.SetPtEtaPhiE(  thejet->pt() , thejet->eta() , thejet->phi() , thejet->energy() );
                        particles_LorentzVector.push_back( jet_lorentzVector );
                        particles_RhoEtaPhiVector.push_back( math::RhoEtaPhiVector( thejet->pt(), thejet->eta(), thejet->phi() ) );

                        bDiscriminatorValue = thejet->bDiscriminator( bTag_.c_str() );

                        if( bDiscriminatorValue > bDiscriminator_[1] ) {
                            tagBJets.push_back( thejet );
                            n_mbjets++;
                            MediumBJetVect.push_back( thejet ); 
                            MediumBJetVect_EtaSorted.push_back( thejet );
                            MediumBJetVect_PtSorted.push_back( thejet );
                            MediumBJetVect_BSorted.push_back( thejet );
                            std::sort(MediumBJetVect_EtaSorted.begin(),MediumBJetVect_EtaSorted.end(),GreaterByEta()); 
                            std::sort(MediumBJetVect_PtSorted.begin(),MediumBJetVect_PtSorted.end(),GreaterByPt()); 
                            std::sort(MediumBJetVect_BSorted.begin(),MediumBJetVect_BSorted.end(),GreaterByBTagging(bTag_.c_str())); 
                        }
                        else{
                            nontagBJets.push_back( thejet );
                            n_ljets++;
                            LightJetVect.push_back( thejet ); 
                            LightJetVect_EtaSorted.push_back( thejet );
                            LightJetVect_PtSorted.push_back( thejet );
                            LightJetVect_BSorted.push_back( thejet );
                            std::sort(LightJetVect_EtaSorted.begin(),LightJetVect_EtaSorted.end(),GreaterByEta()); 
                            std::sort(LightJetVect_PtSorted.begin(),LightJetVect_PtSorted.end(),GreaterByPt()); 
                            std::sort(LightJetVect_BSorted.begin(),LightJetVect_BSorted.end(),GreaterByBTagging(bTag_.c_str())); 
                        }

                        SelJetVect.push_back( thejet ); 
                        SelJetVect_EtaSorted.push_back( thejet );
                        SelJetVect_PtSorted.push_back( thejet );
                        SelJetVect_BSorted.push_back( thejet );
                        std::sort(SelJetVect_EtaSorted.begin(),SelJetVect_EtaSorted.end(),GreaterByEta()); 
                        std::sort(SelJetVect_PtSorted.begin(),SelJetVect_PtSorted.end(),GreaterByPt()); 
                        std::sort(SelJetVect_BSorted.begin(),SelJetVect_BSorted.end(),GreaterByBTagging(bTag_.c_str())); 

                    }//end of jets loop
                    numElectronJetsdR.push_back( deltaRElectronJetcount );
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


            for( unsigned num = 0; num < numMuonJetsdR.size(); num++ ) {
                int check = numMuonJetsdR.at( num );
                if( check >= jetsNumberThreshold_ ) {muonJets = true;}
            }

            for( unsigned num = 0; num < numElectronJetsdR.size(); num++ ) {
                int check = numElectronJetsdR.at( num );
                if( check >= jetsNumberThreshold_ ) {ElectronJets = true;}
            }

            
            std::cout << " THQLeptonicTagProducer tagBJets.size()=" << tagBJets.size()
                      << " tagJets.size()=" << tagJets.size()
                      << " photonSelection=" << photonSelection
                      << " tagMuons.size()=" << tagMuons.size() << " muonJets=" << muonJets
                      << " tagElectrons.size()="<< tagElectrons.size() << " ElectronJets=" << ElectronJets
                      << std::endl;
            
            EventShapeVariables shapeVars(particles_RhoEtaPhiVector);
            std::cout  << "aplanarity: "<<shapeVars.aplanarity()<<std::endl;
            eventshapes.set( shapeVars.aplanarity() ,
                             shapeVars.C() ,
                             shapeVars.circularity(),
                             shapeVars.D() ,
                             shapeVars.isotropy(),
                             shapeVars.sphericity() );

            n_jets = tagJets.size() ; 

            if( tagBJets.size() >= bjetsNumberThreshold_ && tagJets.size() >= jetsNumberThreshold_ && photonSelection
                && ( ( (tagMuons.size() == 1 && muonJets) and  (tagElectrons.size() == 0 && !ElectronJets) )  || ( (tagMuons.size() == 0 && !muonJets)  and  (tagElectrons.size() == 1 && ElectronJets) ) ) ) {
                if( tagElectrons.size() > 0 && ElectronJets ) {
                    std::cout << "including electron weights" << std::endl;
                    thqltags_obj.includeWeights( *tagElectrons[0] );
                } else if( tagMuons.size() > 0 && muonJets ) {
                    std::cout << "including muon weights" << std::endl;
                    thqltags_obj.includeWeights( *tagMuons[0] );
                }

                //TLorentzVector bL,fwdJL;
                bL.SetPtEtaPhiE( MediumBJetVect_BSorted[0]->pt(), MediumBJetVect_BSorted[0]->eta(), MediumBJetVect_BSorted[0]->phi(), MediumBJetVect_BSorted[0]->energy());
                fwdJL.SetPtEtaPhiE( SelJetVect_EtaSorted[0]->pt(),SelJetVect_EtaSorted[0]->eta(), SelJetVect_EtaSorted[0]->phi(), SelJetVect_EtaSorted[0]->energy());
                
                if (std::abs(fwdJL.Mag() - bL.Mag()) < std::numeric_limits<float>::epsilon() )
                    fwdJL.SetPtEtaPhiE( SelJetVect_EtaSorted[1]->pt(),SelJetVect_EtaSorted[1]->eta(), SelJetVect_EtaSorted[1]->phi(), SelJetVect_EtaSorted[1]->energy());   

                    
                flashgg::SemiLepTopQuark singletop(bL, metL, lepton.LorentzVector(), fwdJL,fwdJL);
                //met.SetLorentzVector(singletop.getMET());
                metL = singletop.getMET() ;
                jprime_eta  = fwdJL.Eta();
                met_pt = metL.Pt();
                //met.SetLorentzVector(flashgg::SemiLepTopQuark(bL, met.LorentzVector(), lepton.LorentzVector(), fwdJL,fwdJL).getMET());
                //if (LightJetVect_BSorted.size()>1)
                //std::cout<< " met meta: " << met.LorentzVector().X() << " "<< met.LorentzVector().E()<< " "<< LightJetVect_EtaSorted[0]->eta()<< " "<<LightJetVect_EtaSorted[1]->eta()<< std::endl; 
                
                
                //TVector3 jprimev3;
                //jprime.set( LightJetVect_EtaSorted[0]->pt() , 
                //            LightJetVect_EtaSorted[0]->eta(),
                //            LightJetVect_EtaSorted[0]->phi());
                
                //jprime.set( 0.,
                //            0.,
                //            0.);

                //particles_LorentzVector.push_back( met.LorentzVector() ) ;

                
                FoxWolfram fwam( particles_LorentzVector );

                std::vector< particleinfo*> allfoxwolfs = {&foxwolf1 , &foxwolf2 };
                for(uint ifw = 1 ; ifw < allfoxwolfs.size()+1 ; ifw++)
                    allfoxwolfs[ifw-1]->set( fwam.getMoment( FoxWolfram::SHAT , ifw ),
                                             fwam.getMoment( FoxWolfram::PT , ifw ),
                                             fwam.getMoment( FoxWolfram::ETA , ifw ),
                                             fwam.getMoment( FoxWolfram::PSUM , ifw ),
                                             fwam.getMoment( FoxWolfram::PZ , ifw ),
                                             fwam.getMoment( FoxWolfram::ONE , ifw ) );
                std::cout<< "fox:" << foxwolf1.another<<std::endl;
                
                
                if (MVAMethod_ != "") 
                    thqLeptonicMvaResult_value_ = thqLeptonicMva_->EvaluateMVA( MVAMethod_.c_str() );
                
                //std::cout<< "fox:" << foxwolf1.another<< " "<< thqLeptonicMvaResult_value_ << std::endl;
                
                thqltags_obj.includeWeights( *dipho );
                //thqltags_obj.setTHQLeptonicMVA( thqleptonic_mvares );
                thqltags_obj.setJets( tagJets );
                thqltags_obj.setBJets( tagBJets );
                thqltags_obj.setLightJets( nontagBJets );
                thqltags_obj.setMuons( tagMuons );
                thqltags_obj.setElectrons( tagElectrons );
                thqltags_obj.setDiPhotonIndex( diphoIndex );
                thqltags_obj.setSystLabel( systLabel_ );
                thqltags->push_back( thqltags_obj );
                
                //std::vector<edm::Ptr<flashgg::Jet> > nontagBJets;   
                //std::vector<edm::Ptr<flashgg::Jet> >::iterator it;
                //std::sort (tagBJets.begin(), tagBJets.end());
                //std::sort (tagJets.begin(),tagJets.end());   // 10 20 30 40 50
  

                //it=std::set_difference (tagJets.begin(), tagJets.end(), tagBJets.begin(), tagBJets.end(), nontagBJets.begin());
                /*
                if( ! evt.isRealData()  )
                {
                    
                    for( unsigned int jetCollectionIndex = 0; jetCollectionIndex<nontagBJets.size(); ++jetCollectionIndex ) 
                    {
                        //edm::Ptr<flashgg::Jet> thejet = Jets[jetCollectionIndex]->ptrAt( candIndex_outer );
                        edm::Ptr<flashgg::Jet> thejet = tagJets[jetCollectionIndex];
                        //std::cout<< (*it)->pt()<<std::endl;()
                        if( thejet->pt() > pt_leadjet ) 
                        {
                            index_subleadjet = index_leadjet;
                            pt_subleadjet = pt_leadjet;
                            index_leadjet = jetCollectionIndex;
                            pt_leadjet = thejet->pt();
                        } 
                        else if( thejet->pt() > pt_subleadjet )
                        {
                            index_subleadjet = jetCollectionIndex;
                            pt_subleadjet    = thejet->pt();
                        }
                        }
                    std::cout<< index_subleadjet <<" "<< index_leadjet <<" "<< pt_subleadjet << pt_leadjet<<" "<<index_leadq<<  " "<< index_subleadq << " "<< index_subsubleadq << " "<<pt_leadq << " "<<pt_subleadq << " "<<pt_subsubleadq<<std::endl;
                    //edm::Ptr<flashgg::Jet> theleadjet = Jets[jetCollectionIndex]->ptrAt( index_leadjet);
                    
                    //edm::Ptr<flashgg::Jet> thesubleadjet = Jets[jetCollectionIndex]->ptrAt( index_subleadjet );
                    edm::Ptr<flashgg::Jet>  theleadjet; 
                    if (nontagBJets.size()>0)  theleadjet = nontagBJets[ index_leadjet];
                    else   theleadjet = tagJets[ 0 ];
                    //edm::Ptr<flashgg::Jet> thesubleadjet = nontagBJets[ index_subleadjet ];       
                    VBFTagTruth truth_obj;
                    unsigned int index_gp_leadjet = std::numeric_limits<unsigned int>::max();
                    unsigned int index_gp_subleadjet = std::numeric_limits<unsigned int>::max();
                    unsigned int index_gp_leadphoton = std::numeric_limits<unsigned int>::max();
                    unsigned int index_gp_subleadphoton = std::numeric_limits<unsigned int>::max();
                    unsigned int index_gj_leadjet = std::numeric_limits<unsigned int>::max();
                    unsigned int index_gj_subleadjet = std::numeric_limits<unsigned int>::max();
                    
                    float dr_gp_leadjet = 999.;
                    float dr_gp_subleadjet = 999.;
                    float dr_gp_leadphoton = 999.;
                    float dr_gp_subleadphoton = 999.;
                    float dr_gj_leadjet = 999.;
                    float dr_gj_subleadjet = 999.;
                    for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) 
                    {
                        edm::Ptr<reco::GenParticle> part = genParticles->ptrAt( genLoop );
                        if( part->isHardProcess() ) 
                        {
                            
                            float dr = deltaR( theleadjet->eta(), theleadjet ->phi(), part->eta(), part->phi() );
                            if( dr < dr_gp_leadjet ) 
                            {
                                dr_gp_leadjet = dr;
                                index_gp_leadjet = genLoop;
                            }
                            //dr = deltaR( thesubleadjet->eta(), thesubleadjet->phi(), part->eta(), part->phi() );
                            //if( dr < dr_gp_subleadjet ) 
                            // {
                            //   dr_gp_subleadjet = dr;
                            //   index_gp_subleadjet = genLoop;
                            //}
                            dr = deltaR( thqltags_obj.diPhoton()->leadingPhoton()->eta(), thqltags_obj.diPhoton()->leadingPhoton()->phi(), part->eta(), part->phi() );
                            if( dr < dr_gp_leadphoton ) 
                            {
                                dr_gp_leadphoton = dr;
                                index_gp_leadphoton = genLoop;
                            }
                            dr = deltaR( thqltags_obj.diPhoton()->subLeadingPhoton()->eta(), thqltags_obj.diPhoton()->subLeadingPhoton()->phi(), part->eta(), part->phi() );
                            if( dr < dr_gp_subleadphoton ) 
                            {
                                dr_gp_subleadphoton = dr;
                                index_gp_subleadphoton = genLoop;
                            }
                        }
                    }
                    
                    std::cout<< index_gp_leadjet <<" "<< index_gp_subleadjet <<" "<< index_gp_leadphoton << index_gp_subleadphoton <<" "<<index_gj_leadjet <<  " "<< index_gj_subleadjet << " "<< dr_gp_leadjet << " "<<dr_gp_subleadjet << " "<<dr_gp_leadphoton << " "<<dr_gp_subleadphoton<<" "<< dr_gj_leadjet <<" "<< dr_gj_subleadjet<<std::endl;
                    
                    if( index_gp_leadjet < std::numeric_limits<unsigned int>::max() ) { truth_obj.setClosestParticleToLeadingJet( genParticles->ptrAt( index_gp_leadjet ) ); }
                    if( index_gp_subleadjet < std::numeric_limits<unsigned int>::max() ) { truth_obj.setClosestParticleToSubLeadingJet( genParticles->ptrAt( index_gp_subleadjet ) ); }
                    if( index_gp_leadphoton < std::numeric_limits<unsigned int>::max() ) { truth_obj.setClosestParticleToLeadingPhoton( genParticles->ptrAt( index_gp_leadphoton ) ); }
                    if( index_gp_subleadphoton < std::numeric_limits<unsigned int>::max() ) { truth_obj.setClosestParticleToSubLeadingPhoton( genParticles->ptrAt( index_gp_subleadphoton ) ); }
                    
                    for( unsigned int gjLoop = 0 ; gjLoop < genJets->size() ; gjLoop++ ) 
                    {
                        edm::Ptr <reco::GenJet> gj = genJets->ptrAt( gjLoop );
                        float dr = deltaR( theleadjet->eta(), theleadjet->phi(), gj->eta(), gj->phi() );
                        if( dr < dr_gj_leadjet ) 
                        {
                            dr_gj_leadjet = dr;
                            index_gj_leadjet = gjLoop;
                        }
                        //dr = deltaR( thesubleadjet->eta(), thesubleadjet->phi(), gj->eta(), gj->phi() );
                        //if( dr < dr_gj_subleadjet ) 
                        //{
                        //   dr_gj_subleadjet = dr;
                        //index_gj_subleadjet = gjLoop;
                        //}
                    }
                    if( index_gj_leadjet < std::numeric_limits<unsigned int>::max() ) { truth_obj.setClosestGenJetToLeadingJet( genJets->ptrAt( index_gj_leadjet ) ); }
                    if( index_gj_subleadjet < std::numeric_limits<unsigned int>::max() ) { truth_obj.setClosestGenJetToSubLeadingJet( genJets->ptrAt( index_gj_subleadjet ) ); }
                    
                    if( index_leadq < std::numeric_limits<unsigned int>::max() ) { truth_obj.setLeadingParton( genParticles->ptrAt( index_leadq ) ); }
                    if( index_subleadq < std::numeric_limits<unsigned int>::max() ) { truth_obj.setSubLeadingParton( genParticles->ptrAt( index_subleadq ) ); }
                    if( index_subsubleadq < std::numeric_limits<unsigned int>::max()) { truth_obj.setSubSubLeadingParton( genParticles->ptrAt( index_subsubleadq ));}
                    // --------
                    //Partons
                    //Lead
                    std::vector<edm::Ptr<reco::GenParticle>> ptOrderedPartons;
                    for (unsigned int genLoop(0);genLoop < genParticles->size();genLoop++) 
                    {
                        edm::Ptr<reco::GenParticle> gp = genParticles->ptrAt(genLoop);
                        bool isGluon = abs( gp->pdgId() ) < 7 && gp->numberOfMothers() == 0;
                        bool isQuark = gp->pdgId() == 21 && gp->numberOfMothers() == 0;
                        if (isGluon || isQuark) {
                            unsigned int insertionIndex(0);
                            for (unsigned int parLoop(0);parLoop<ptOrderedPartons.size();parLoop++) {
                                if (gp->pt() < ptOrderedPartons[parLoop]->pt()) { insertionIndex = parLoop + 1; }
                            }
                            ptOrderedPartons.insert( ptOrderedPartons.begin() + insertionIndex, gp);
                        }
                    }
                    if ( ptOrderedPartons.size() > 0 ) 
                    {
                        float dr(999.0);
                        unsigned pIndex(0);
                        for (unsigned partLoop(0);partLoop<ptOrderedPartons.size();partLoop++) {
                            float deltaR_temp = deltaR(theleadjet->eta(),theleadjet->phi(),
                                                       ptOrderedPartons[partLoop]->eta(),ptOrderedPartons[partLoop]->phi());
                            if (deltaR_temp < dr) {dr = deltaR_temp; pIndex = partLoop;}
                        }
                        truth_obj.setClosestPartonToLeadingJet( ptOrderedPartons[pIndex] );
                    }
                    
                    //Sublead
                                        truth_obj.setGenPV( higgsVtx );
                    truths->push_back( truth_obj );
                    thqltags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<VBFTagTruth> >( rTagTruth, idx++ ) ) );
                }
                */
            }//thq tag
            else {
                std::cout << " THQLeptonicTagProducer NO TAG " << std::endl;
            }
            
        }//diPho loop end !
        evt.put( thqltags );
        //evt.put( truths );

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

