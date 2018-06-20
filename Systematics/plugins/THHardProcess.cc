#include "flashgg/Systematics/interface/THHardProcess.h"

THHardProcess::THHardProcess(const edm::ParameterSet& iConfig) :
  genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag>( "genParticleTag" ) ) ),
  genInfoToken_( consumes<GenEventInfoProduct>( iConfig.getParameter<InputTag>( "genInfoTag" ) ) ),
  debug_( iConfig.getParameter<bool>( "debug" ) ),
  token_lhe( consumes<LHEEventProduct>( InputTag( "externalLHEProducer" ) ) )
{
  npass = 0;
  nfail = 0;
  ntot = 0;
  wpass = 0.;
  wfail = 0.;
  wtot = 0.;

  
}

void THHardProcess::beginJob(){
  CTCVWeightedVariables["jprimeeta"] = new CTCVWeightedVariable("jPrimeEta" , "#eta_{j'}" , 40 , -5,  5 );
  CTCVWeightedVariables["jprimept"] = new CTCVWeightedVariable("jPrimePt" , "p^{j'}_{T}" , 50 , 0 , 200 );
  
  CTCVWeightedVariables["lightbeta"] = new CTCVWeightedVariable("lightbEta" , "#eta_{light-b}" , 40 , -5,  5 );
  CTCVWeightedVariables["lightbpt"] = new CTCVWeightedVariable("lightbPt" , "light-b p_{T}" , 50 , 0 , 200 );
  
  CTCVWeightedVariables["higgseta"] = new CTCVWeightedVariable("higgsEta" , "#eta_{h}" , 40 , -5,  5 );
  CTCVWeightedVariables["higgspt"] = new CTCVWeightedVariable("higgsPt" , "higgs p_{T}" , 50 , 0 , 200 );
  
  CTCVWeightedVariables["topeta"] = new CTCVWeightedVariable("topeta" , "#eta_{top}" , 40 , -5,  5 );
  CTCVWeightedVariables["toppt"] = new CTCVWeightedVariable("toppt" , "top p_{T}" , 50 , 0 , 200 );
}

THHardProcess::~THHardProcess()
{
}

void THHardProcess::endJob()
{
  for( auto& a : CTCVWeightedVariables)
    a.second->Write();
  std::cout << "THHardProcess endJob npassHiggsTop/nfailFWD/ntot: " << npass << "/" << nfail << "/" << ntot << " wpassHiggsTop/wfailFWD/wtot: " << wpass << "/" << wfail << "/" << wtot << std::endl;
}


// member functions
//

// ------------ method called on each new Event  ------------
bool
THHardProcess::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<LHEEventProduct> product_lhe;
  vector< double > CtCvWeights ;
  iEvent.getByToken(token_lhe, product_lhe);
  for (uint i = 446 ; i < product_lhe->weights().size() ; i++)
    CtCvWeights.push_back(product_lhe->weights()[i].wgt/product_lhe->originalXWGTUP () );
  

  float weight = 1.;
  edm::Handle<GenEventInfoProduct> genInfo;
  iEvent.getByToken(genInfoToken_, genInfo);
  if( genInfo.isValid() ) {
    const auto &weights = genInfo->weights();
    if( ! weights.empty() ) {
      weight *= weights[0];
    }
  }


  int numGoodPhotons = 0;

  Handle<View<reco::GenParticle> > genParticles;
  iEvent.getByToken( genParticleToken_, genParticles );
  
  Ptr<reco::GenParticle> top ;
  Ptr<reco::GenParticle> higgs ;
  Ptr<reco::GenParticle> fwdJet ;
  Ptr<reco::GenParticle> light_b;

  bool higgsSet = false;
  bool topSet = false;
  bool fwdSet = false;
  bool haslightb = false;
  for( unsigned int i = 0 ; i < genParticles->size(); i++ ) {
    Ptr<reco::GenParticle> gen = genParticles->ptrAt(i);

    if( gen->isHardProcess() ){
      const reco::GenParticle* mother = GetMother( gen.get() );
      const reco::GenParticle* grandmom = GetMother( mother );
      int mother_id = 0 ;
      if( mother != NULL )
	mother_id = mother->pdgId() ;
      int gmother_id = 0 ;
      if( grandmom != NULL )
	gmother_id = grandmom->pdgId() ;
      //cout << "\t" << gmother_id << " --> " << mother_id << " --> " << gen->pdgId() << endl;

      if( gmother_id == 0 && mother_id == 0 )
	{
	  if( abs(gen->pdgId()) == 5 ){
	    light_b = gen;
	    haslightb = true;
	  }else if( abs(gen->pdgId()) < 5 ){
	    fwdJet = gen ;
	    fwdSet = true;
	  }
	}
    }
    

    if( gen->isLastCopy() && abs( gen->pdgId() )==6 ){
      topSet = true;
      top = gen;
      continue;
    }
    if( gen->isLastCopy() && abs( gen->pdgId() )==25) {
      higgsSet = true;
      higgs = gen;
      continue;
    }


  }


  if( topSet ){
    CTCVWeightedVariables["topeta"]->Fill( top->eta() , CtCvWeights );
    CTCVWeightedVariables["toppt"]->Fill( top->pt() , CtCvWeights );
  }

  if( higgsSet ){
    CTCVWeightedVariables["higgseta"]->Fill( higgs->eta() , CtCvWeights );
    CTCVWeightedVariables["higgspt"]->Fill( higgs->pt() , CtCvWeights );
  }

  if( haslightb ){
    CTCVWeightedVariables["lightbeta"]->Fill( light_b->eta() , CtCvWeights );
    CTCVWeightedVariables["lightbpt"]->Fill( light_b->pt() , CtCvWeights );
  }
  if( fwdSet ){
    CTCVWeightedVariables["jprimeeta"]->Fill( fwdJet->eta() , CtCvWeights );
    CTCVWeightedVariables["jprimept"]->Fill( fwdJet->pt() , CtCvWeights );
  }
  

  if (debug_) {
    std::cout << " THHardProcess debug top pt = " << top->pt() << " status=" << top->status() 
	      << "\n higgs pt = " << higgs->pt() << " status=" << higgs->status() << std::endl;
    if( fwdSet )
      std::cout << " fwdJet pt = " << fwdJet->pt() << endl;
  }


  // counting
  ntot += 1;
  wtot += weight;

  if (debug_) {
    std::cout << "THHardProcess debug numGoodPhotons=" << numGoodPhotons << " pass=" << (numGoodPhotons==2) << std::endl;
  }

  return (fwdSet);
}

DEFINE_FWK_MODULE(THHardProcess);

