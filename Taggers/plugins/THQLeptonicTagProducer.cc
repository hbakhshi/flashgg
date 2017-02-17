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
    std::vector<edm::InputTag> inputTagJets_;
    EDGetTokenT<View<Electron> > electronToken_;
    EDGetTokenT<View<flashgg::Muon> > muonToken_;
    EDGetTokenT<View<DiPhotonMVAResult> > mvaResultToken_;
    EDGetTokenT<View<Photon> > photonToken_;
    EDGetTokenT<View<reco::Vertex> > vertexToken_;
    EDGetTokenT<View<flashgg::Met> > METToken_;
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

    double deltaRPhoElectronThreshold_;
    double Zmass_;
    double deltaMassElectronZThreshold_;

    bool hasGoodElec = false;  bool hasVetoElec = false;
    bool hasGoodMuons = false;

    unique_ptr<TMVA::Reader> thqLeptonicMva_;
    FileInPath thqLeptonicMVAweightfile_;
    string  MVAMethod_;
    float thqLeptonicMvaResult_value_, topMass;

    std::vector< TLorentzVector > particles_LorentzVector;
    std::vector< math::RhoEtaPhiVector > particles_RhoEtaPhiVector;
        
    TLorentzVector metL, bL,fwdJL, G1, G2;  //temp solution: make met, bjet & jprime global TLorentzVectors

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

        
    int LeptonType;
    std::vector<edm::Ptr<flashgg::Jet> > SelJetVect; std::vector<edm::Ptr<flashgg::Jet> > SelJetVect_EtaSorted; std::vector<edm::Ptr<flashgg::Jet> > SelJetVect_PtSorted; std::vector<edm::Ptr<flashgg::Jet> > SelJetVect_BSorted;
    std::vector<edm::Ptr<flashgg::Jet> > MediumBJetVect, MediumBJetVect_PtSorted;
    std::vector<edm::Ptr<flashgg::Jet> > LooseBJetVect, LooseBJetVect_PtSorted ;
    std::vector<edm::Ptr<flashgg::Jet> > TightBJetVect, TightBJetVect_PtSorted;


    edm::Ptr<flashgg::Jet> fwdJet;
    edm::Ptr<flashgg::Jet> bJet  ;
    void topReco( std::vector<edm::Ptr<flashgg::Jet> >* bjets ){
      topMass = -100.;
      thqLeptonicMvaResult_value_ = -100.;

      if ( bjets->size() < 1 || SelJetVect.size() < 2 || LeptonType == 0){
	return ;
      }
      fwdJet = SelJetVect_EtaSorted[0];
      bJet = bjets->at(0);
      if( fwdJet == bJet )
	fwdJet = SelJetVect_EtaSorted[1] ;

      
      bL.SetPtEtaPhiE( bJet->pt(), bJet->eta(), bJet->phi(), bJet->energy());
      fwdJL.SetPtEtaPhiE( fwdJet->pt(),fwdJet->eta(), fwdJet->phi(), fwdJet->energy());


      flashgg::SemiLepTopQuark singletop(bL, metL, lepton.LorentzVector(), fwdJL,fwdJL);
      n_jets = SelJetVect.size();
      metL = singletop.getMET() ;
      jprime_eta  = fabs( fwdJL.Eta() );
      met_pt = metL.Pt();

      topMass = singletop.top().M() ;

      if (MVAMethod_ != "") 
	thqLeptonicMvaResult_value_ = thqLeptonicMva_->EvaluateMVA( MVAMethod_.c_str() );

    };

    
    //MVA INPUTS
    float  n_jets = 0;
    float jprime_eta,met_pt;

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
	lorentzVector_.SetPtEtaPhiM(pt,eta,phi,other_);
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
        
    particleinfo lepton ,  eventshapes;
    particleinfo foxwolf1 ; // foxwolf2 , foxwolf1Met, foxwolf2Met ;
  };

  THQLeptonicTagProducer::THQLeptonicTagProducer( const ParameterSet &iConfig ) :
    diPhotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
    inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
    electronToken_( consumes<View<flashgg::Electron> >( iConfig.getParameter<InputTag>( "ElectronTag" ) ) ),
    muonToken_( consumes<View<flashgg::Muon> >( iConfig.getParameter<InputTag>( "MuonTag" ) ) ),
    mvaResultToken_( consumes<View<flashgg::DiPhotonMVAResult> >( iConfig.getParameter<InputTag> ( "MVAResultTag" ) ) ),
    vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
    METToken_( consumes<View<flashgg::Met> >( iConfig.getParameter<InputTag> ( "METTag" ) ) ),
    rhoTag_( consumes<double>( iConfig.getParameter<InputTag>( "rhoTag" ) ) ),
    systLabel_( iConfig.getParameter<string> ( "SystLabel" ) ),
    MVAMethod_    ( iConfig.getParameter<string> ( "MVAMethod"    ) )
  {
    double default_Zmass_ = 91.9;
    double default_deltaMassElectronZThreshold_ = 10.;


    vector<double> default_electronEtaCuts_;
    default_electronEtaCuts_.push_back( 1.4442 );
    default_electronEtaCuts_.push_back( 1.566 );
    default_electronEtaCuts_.push_back( 2.5 );

    leptonEtaThreshold_ = iConfig.getParameter<double>( "leptonEtaThreshold" );
    leptonPtThreshold_ = iConfig.getParameter<double>( "leptonPtThreshold" );
    electronEtaThresholds_ = iConfig.getParameter<vector<double > >( "electronEtaThresholds");
    leadPhoOverMassThreshold_ = iConfig.getParameter<double>( "leadPhoOverMassThreshold" );
    subleadPhoOverMassThreshold_ = iConfig.getParameter<double>( "subleadPhoOverMassThreshold" );
    MVAThreshold_ = iConfig.getParameter<double>( "MVAThreshold" );
    deltaRLepPhoThreshold_ = iConfig.getParameter<double>( "deltaRLepPhoThreshold" );
    deltaRJetLepThreshold_ = iConfig.getParameter<double>( "deltaRJetLepThreshold" );
    jetsNumberThreshold_ = iConfig.getParameter<double>( "jetsNumberThreshold" );
    bjetsNumberThreshold_ = iConfig.getParameter<double>( "bjetsNumberThreshold" );
    jetPtThreshold_ = iConfig.getParameter<double>( "jetPtThreshold" );
    jetEtaThreshold_ = iConfig.getParameter<double>( "jetEtaThreshold" );

    deltaRJetLeadPhoThreshold_ = iConfig.getParameter<double>( "deltaRJetLeadPhoThreshold" );
    deltaRJetSubLeadPhoThreshold_ = iConfig.getParameter<double>( "deltaRJetSubLeadPhoThreshold" );

    electronEtaThresholds_ = iConfig.getUntrackedParameter<vector<double > >( "electronEtaCuts",default_electronEtaCuts_);
    bDiscriminator_ = iConfig.getParameter<vector<double > >( "bDiscriminator" );
    bTag_ = iConfig.getParameter<string>( "bTag" );

    muPFIsoSumRelThreshold_ = iConfig.getParameter<double>( "muPFIsoSumRelThreshold" );
    PhoMVAThreshold_ = iConfig.getParameter<double>( "PhoMVAThreshold" );
    DeltaRTrkElec_ = iConfig.getParameter<double>( "DeltaRTrkElec" );

    deltaRPhoElectronThreshold_ = iConfig.getParameter<double>( "deltaRPhoElectronThreshold" );
    Zmass_ = iConfig.getUntrackedParameter<double>( "Zmass_", default_Zmass_ );
    deltaMassElectronZThreshold_ = iConfig.getUntrackedParameter<double>( "deltaMassElectronZThreshold_", default_deltaMassElectronZThreshold_ );

     
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
  }

  void THQLeptonicTagProducer::produce( Event &evt, const EventSetup & )

  {
    JetCollectionVector Jets( inputTagJets_.size() );
    for( size_t j = 0; j < inputTagJets_.size(); ++j ) {
      evt.getByToken( tokenJets_[j], Jets[j] );
    }

    edm::Handle<double>  rho;
    evt.getByToken(rhoTag_,rho);
    float rho_    = *rho;

    Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
    evt.getByToken( diPhotonToken_, diPhotons );

    Handle<View<flashgg::Muon> > theMuons;
    evt.getByToken( muonToken_, theMuons );

    Handle<View<flashgg::Electron> > theElectrons;
    evt.getByToken( electronToken_, theElectrons );

    Handle<View<flashgg::DiPhotonMVAResult> > mvaResults;
    evt.getByToken( mvaResultToken_, mvaResults );

    std::auto_ptr<vector<THQLeptonicTag> > thqltags( new vector<THQLeptonicTag> );

    Handle<View<reco::Vertex> > vertices;
    evt.getByToken( vertexToken_, vertices );

    Handle<View<flashgg::Met> > METs;
    evt.getByToken( METToken_, METs );


    assert( diPhotons->size() == mvaResults->size() );

    bool photonSelection = false;
    double idmva1 = 0.;
    double idmva2 = 0.;

    for( unsigned int diphoIndex = 0; diphoIndex < diPhotons->size(); diphoIndex++ ) {

      hasGoodElec = false; hasVetoElec = false;
      hasGoodMuons = false;
            
      unsigned int jetCollectionIndex = diPhotons->ptrAt( diphoIndex )->jetCollectionIndex();
            
      edm::Ptr<flashgg::DiPhotonCandidate> dipho = diPhotons->ptrAt( diphoIndex );
      edm::Ptr<flashgg::DiPhotonMVAResult> mvares = mvaResults->ptrAt( diphoIndex );


      flashgg::THQLeptonicTag thqltags_obj( dipho, mvares );

      if( dipho->leadingPhoton()->pt() < ( dipho->mass() )*leadPhoOverMassThreshold_ ) { continue; }

      if( dipho->subLeadingPhoton()->pt() < ( dipho->mass() )*subleadPhoOverMassThreshold_ ) { continue; }


      idmva1 = dipho->leadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
      idmva2 = dipho->subLeadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );

      if( idmva1 <= PhoMVAThreshold_ || idmva2 <= PhoMVAThreshold_ ) { continue; }

      if( mvares->result < MVAThreshold_ ) { continue; }

      photonSelection = true;
            
        
      G1.SetPtEtaPhiM( diPhotons->ptrAt( diphoIndex )->leadingPhoton()->pt(),
		       diPhotons->ptrAt( diphoIndex )->leadingPhoton()->eta(),
		       diPhotons->ptrAt( diphoIndex )->leadingPhoton()->phi() , 
		       0 );
      particles_LorentzVector.push_back( G1 );
            
      G2.SetPtEtaPhiM( diPhotons->ptrAt( diphoIndex )->subLeadingPhoton()->pt(),
		       diPhotons->ptrAt( diphoIndex )->subLeadingPhoton()->eta(),
		       diPhotons->ptrAt( diphoIndex )->subLeadingPhoton()->phi(),
		       0 );
      particles_LorentzVector.push_back(G2);

      particles_RhoEtaPhiVector.push_back( math::RhoEtaPhiVector(G1.Pt(), G1.Eta() , G1.Phi() ) );
      particles_RhoEtaPhiVector.push_back( math::RhoEtaPhiVector(G2.Pt(), G2.Eta() , G2.Phi() ) );

      if( METs->size() != 1 ) { std::cout << "WARNING - #MET is not 1" << std::endl;}
      Ptr<flashgg::Met> theMET = METs->ptrAt( 0 );

      //const pat::MET &met_ = METs->front();
      //std::cout << met_.pt() <<std::endl;
      metL.SetPtEtaPhiE( theMET->getCorPt(),
			 theMET->eta(),
			 theMET->getCorPhi(),
			 theMET->energy()
			 ) ; 


      std::vector<edm::Ptr<flashgg::Muon> > goodLooseMuons = selectMuons( theMuons->ptrs(), dipho, vertices->ptrs(), leptonEtaThreshold_ , leptonPtThreshold_,
									  2.0, deltaRLepPhoThreshold_, deltaRLepPhoThreshold_ , true);

      std::vector<edm::Ptr<flashgg::Muon> > goodMuons = selectMuons( theMuons->ptrs(), dipho, vertices->ptrs(), leptonEtaThreshold_ , leptonPtThreshold_,
								     muPFIsoSumRelThreshold_, deltaRLepPhoThreshold_, deltaRLepPhoThreshold_ , false);

      std::vector<edm::Ptr<Electron> > vetoElectronsNonIso = selectStdElectrons(theElectrons->ptrs(), dipho, vertices->ptrs(), leptonPtThreshold_,  electronEtaThresholds_ ,
									  0,0,
									  deltaRPhoElectronThreshold_,DeltaRTrkElec_,deltaMassElectronZThreshold_ , rho_, true ); //evt.isRealData()

      
      std::vector<edm::Ptr<Electron> > vetoElectrons = selectStdElectrons(theElectrons->ptrs(), dipho, vertices->ptrs(), leptonPtThreshold_,  electronEtaThresholds_ ,
									  0,1,
									  deltaRPhoElectronThreshold_,DeltaRTrkElec_,deltaMassElectronZThreshold_ , rho_, true ); //evt.isRealData()
            

      std::vector<edm::Ptr<Electron> > goodElectrons = selectStdElectrons(theElectrons->ptrs(), dipho, vertices->ptrs(), leptonPtThreshold_,  electronEtaThresholds_ ,
									  0,2,
									  deltaRPhoElectronThreshold_,DeltaRTrkElec_,deltaMassElectronZThreshold_ , rho_, true);


      hasGoodElec = ( goodElectrons.size() == 1 ); hasVetoElec = ( vetoElectrons.size() > 0 );
      hasGoodMuons = ( goodMuons.size() == 1 );


      LeptonType = 0; //1 : electron, 2:muon

      
      if( hasGoodMuons && !hasVetoElec){
	LeptonType = 2;
      }
      for( unsigned int muonIndex = 0; muonIndex < goodMuons.size(); muonIndex++ ) {
                
	Ptr<flashgg::Muon> muon = goodMuons[muonIndex];

	thqltags_obj.includeWeights( *goodMuons[muonIndex] );
	
	lepton.set( muon->pt(),
		    muon->eta() ,
		    muon->phi() ,
		    muon->energy(),
		    1. ,
		    muon->charge() );
	particles_LorentzVector.push_back(lepton.LorentzVector());
	particles_RhoEtaPhiVector.push_back( math::RhoEtaPhiVector( lepton.pt, lepton.eta, lepton.phi ) );
	
      }//end of muons loop


      if( hasGoodElec && !hasGoodMuons){
	LeptonType = 1;
      }

      for( unsigned int ElectronIndex = 0; ElectronIndex < goodElectrons.size(); ElectronIndex++ ) {

	thqltags_obj.includeWeights( *goodElectrons[ElectronIndex] );
                
	Ptr<Electron> Electron = goodElectrons[ElectronIndex];
	lepton.set( Electron->pt(),
		    Electron->eta() ,
		    Electron->phi() ,
		    Electron->energy(),
		    1. ,
		    Electron->charge() );
	particles_LorentzVector.push_back(lepton.LorentzVector());
	particles_RhoEtaPhiVector.push_back( math::RhoEtaPhiVector( lepton.pt, lepton.eta, lepton.phi ) );
                
                
      }//end of electron loop





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


	TLorentzVector jet_lorentzVector;
	jet_lorentzVector.SetPtEtaPhiE(  thejet->pt() , thejet->eta() , thejet->phi() , thejet->energy() );
	//std::cout <<  thejet->pt() << " "<< thejet->eta()<<" "<< thejet->phi()<< " "<< thejet->energy() <<std::endl;
	particles_LorentzVector.push_back( jet_lorentzVector );
	particles_RhoEtaPhiVector.push_back( math::RhoEtaPhiVector( thejet->pt(), thejet->eta(), thejet->phi() ) );
                

	double minDrLepton = 999.;
	for(auto mu : goodMuons){
	  float dRJetLepton = deltaR( thejet->eta(), thejet->phi(), mu->eta() , mu->phi() );
	  if( dRJetLepton < minDrLepton ) { minDrLepton = dRJetLepton; }
	}
	for(auto ele : goodElectrons){
	  float dRJetLepton = deltaR( thejet->eta(), thejet->phi(), ele->eta() , ele->phi() );
	  if( dRJetLepton < minDrLepton ) { minDrLepton = dRJetLepton; }
	}

	if( minDrLepton < deltaRJetLepThreshold_) continue;

	double bDiscriminatorValue = thejet->bDiscriminator( bTag_.c_str() );

	if( bDiscriminatorValue > bDiscriminator_[0] ) {
	  LooseBJetVect_PtSorted.push_back( thejet ); 
	  LooseBJetVect.push_back( thejet );
	}

	if( bDiscriminatorValue > bDiscriminator_[1] ) {
	  MediumBJetVect.push_back( thejet ); 
	  MediumBJetVect_PtSorted.push_back( thejet );
	}

	if( bDiscriminatorValue > bDiscriminator_[2] ) {
	  TightBJetVect_PtSorted.push_back( thejet ); 
	  TightBJetVect.push_back( thejet );
	}

                        
	SelJetVect.push_back( thejet ); 
	SelJetVect_EtaSorted.push_back( thejet );
	SelJetVect_PtSorted.push_back( thejet );
	SelJetVect_BSorted.push_back( thejet );
      }//end of jets loop
      std::sort(LooseBJetVect_PtSorted.begin(),LooseBJetVect_PtSorted.end(),GreaterByPt()); 
      std::sort(LooseBJetVect.begin(),LooseBJetVect.end(),GreaterByBTagging(bTag_.c_str())); 

      std::sort(MediumBJetVect_PtSorted.begin(),MediumBJetVect_PtSorted.end(),GreaterByPt()); 
      std::sort(MediumBJetVect.begin(),MediumBJetVect.end(),GreaterByBTagging(bTag_.c_str())); 

      std::sort(TightBJetVect_PtSorted.begin(),TightBJetVect_PtSorted.end(),GreaterByPt()); 
      std::sort(TightBJetVect.begin(),TightBJetVect.end(),GreaterByBTagging(bTag_.c_str())); 

      std::sort(SelJetVect_EtaSorted.begin(),SelJetVect_EtaSorted.end(),GreaterByEta()); 
      std::sort(SelJetVect_PtSorted.begin(),SelJetVect_PtSorted.end(),GreaterByPt()); 
      std::sort(SelJetVect.begin(),SelJetVect.end(),GreaterByBTagging(bTag_.c_str())); 


            



      if( photonSelection ){
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
	std::vector< particleinfo*> allfoxwolfs = {&foxwolf1 };
	for(uint ifw = 1 ; ifw < allfoxwolfs.size()+1 ; ifw++)
	  allfoxwolfs[ifw-1]->set( fwam.getMoment( FoxWolfram::SHAT , ifw ),
				   fwam.getMoment( FoxWolfram::PT , ifw ),
				   fwam.getMoment( FoxWolfram::ETA , ifw ),
				   fwam.getMoment( FoxWolfram::PSUM , ifw ),
				   fwam.getMoment( FoxWolfram::PZ , ifw ),
				   fwam.getMoment( FoxWolfram::ONE , ifw ) );
	//std::cout<< "fox:" << foxwolf1.another<<std::endl;

	thqltags_obj.setrho(rho_);

	thqltags_obj.setLeptonType(LeptonType);
	thqltags_obj.includeWeights( *dipho );

	thqltags_obj.photonWeights = dipho->leadingPhoton()->centralWeight()*dipho->subLeadingPhoton()->centralWeight() ;

	thqltags_obj.setJets( SelJetVect_PtSorted , SelJetVect_EtaSorted);
	thqltags_obj.setBJets( SelJetVect_BSorted );

	thqltags_obj.bTagWeight = 1.0;



	for( auto j : SelJetVect_PtSorted )
	  //for(auto itr = j->weightListBegin() ; itr != j->weightListEnd() ; itr++)
	  //  cout << *itr << endl;
	  if( j->hasWeight("JetBTagCutWeightCentral") )
	      thqltags_obj.bTagWeight *= j->weight( "JetBTagCutWeightCentral" );
	  else
	    cout << "BTag weight is not set in jet" << endl;


	thqltags_obj.setVertices( vertices->ptrs() );

	std::vector <float> a; std::vector <float> b; std::vector <float> c; std::vector <float> d;
	for( unsigned int muonIndex = 0; muonIndex < goodLooseMuons.size(); muonIndex++ ) {
                
	  Ptr<flashgg::Muon> muon = goodLooseMuons[muonIndex];
	  
	  int vtxInd = -1;
	  double dzmin = 9999;
	  for( size_t ivtx = 0 ; ivtx < vertices->ptrs().size(); ivtx++ ) {
	    Ptr<reco::Vertex> vtx = vertices->ptrs()[ivtx];
	    if( !muon->innerTrack() ) continue; 
	    if( fabs( muon->innerTrack()->vz() - vtx->position().z() ) < dzmin ) {                    
	      dzmin = fabs( muon->innerTrack()->vz() - vtx->position().z() );
	      vtxInd = ivtx;
	    }
	  }
	  Ptr<reco::Vertex> best_vtx = vertices->ptrs()[vtxInd]; 
	  a.push_back(muon->muonBestTrack()->dxy(best_vtx->position()));
	  b.push_back(muon->muonBestTrack()->dz(best_vtx->position()));
	  c.push_back(muon->muonBestTrack()->dxy(dipho->vtx()->position()));
	  d.push_back(muon->muonBestTrack()->dz(dipho->vtx()->position()));
	}//end of muons loop

	thqltags_obj.setLeptonVertices( "muon", a, b, c, d) ;
	
	//std::cout << "new vertex !! "<< thqltags_obj.getSubLeadingLeptonVertexDxy( "muon") << std::endl;

	thqltags_obj.setMuons( goodLooseMuons );
	//cout << "nLooseMuons : " << goodLooseMuons.size() << " and nTightMuons : " << goodMuons.size() << " out of : " << theMuons->ptrs().size() << endl;

	a.clear();b.clear();c.clear();d.clear();
	for( unsigned int ElectronIndex = 0; ElectronIndex < vetoElectronsNonIso.size(); ElectronIndex++ ) {
                
	  Ptr<flashgg::Electron> electron = vetoElectronsNonIso[ElectronIndex];
	  
	  int vtxInd = -1;
	  double dzmin = 9999;
	  for( size_t ivtx = 0 ; ivtx < vertices->ptrs().size(); ivtx++ ) {
	    Ptr<reco::Vertex> vtx = vertices->ptrs()[ivtx];
	    if( fabs( electron->gsfTrack()->dz(vtx->position()) ) < dzmin ) {                    
	      dzmin = fabs(electron->gsfTrack()->dz( vtx->position() )); 
	      vtxInd = ivtx;
	    }
	  }
	  Ptr<reco::Vertex> best_vtx = vertices->ptrs()[vtxInd]; 
	  a.push_back(electron->gsfTrack()->dxy(best_vtx->position()));
	  b.push_back(electron->gsfTrack()->dz(best_vtx->position()));
	  c.push_back(electron->gsfTrack()->dxy(dipho->vtx()->position()));
	  d.push_back(electron->gsfTrack()->dz(dipho->vtx()->position()));
	  int elMissedHits = electron->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS);
	  thqltags_obj.setElectronMisHits(elMissedHits);
	}//end of electrons loop

	thqltags_obj.setLeptonVertices( "electron", a, b, c, d) ;

	thqltags_obj.setElectrons( vetoElectronsNonIso );

	thqltags_obj.setDiPhotonIndex( diphoIndex );
	thqltags_obj.setSystLabel( systLabel_ );

	thqltags_obj.setValues( foxwolf1.another , eventshapes.pt , metL.Pt(), metL.Phi()  );

	topReco( &SelJetVect_BSorted );
	thqltags_obj.setMVAres("HighestBTagVal" ,  thqLeptonicMvaResult_value_ , topMass , fwdJet , bJet);

	topReco( &MediumBJetVect_PtSorted );
	thqltags_obj.setMVAres("Medium" ,  thqLeptonicMvaResult_value_ , topMass , fwdJet , bJet);
	thqltags_obj.nMedium_bJets = MediumBJetVect_PtSorted.size();

	topReco( &LooseBJetVect_PtSorted );
	thqltags_obj.setMVAres("Loose" ,  thqLeptonicMvaResult_value_ , topMass , fwdJet , bJet);
	thqltags_obj.nLoose_bJets = LooseBJetVect_PtSorted.size();

	topReco( &TightBJetVect_PtSorted );
	thqltags_obj.setMVAres("Tight" ,  thqLeptonicMvaResult_value_ , topMass , fwdJet , bJet);
	thqltags_obj.nTight_bJets = TightBJetVect_PtSorted.size();

	thqltags->push_back( thqltags_obj );
      }//thq tag
      else {
	if(false)
	  std::cout << " THQLeptonicTagProducer NO TAG " << std::endl;
      }

      n_jets = 0;

      particles_LorentzVector.clear();
      particles_RhoEtaPhiVector.clear();
      SelJetVect.clear(); SelJetVect_EtaSorted.clear(); SelJetVect_PtSorted.clear(); SelJetVect_BSorted.clear();
      LooseBJetVect.clear(); LooseBJetVect_PtSorted.clear(); 
      MediumBJetVect.clear(); MediumBJetVect_PtSorted.clear();
      TightBJetVect.clear(); TightBJetVect_PtSorted.clear();

            
    }//diPho loop end !
        
    evt.put( thqltags );
  }
    
}
typedef flashgg::THQLeptonicTagProducer FlashggTHQLeptonicTagProducer;
DEFINE_FWK_MODULE( FlashggTHQLeptonicTagProducer );

