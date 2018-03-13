#ifndef THHARDPROCESS_h
#define THHARDPROCESS_h

// -*- C++ -*-
//
// Package:    THHardProcess
// Class:      THHardProcess
// 

// system include files
#include <memory>
#include <map>

#include "TCanvas.h"
#include "TH1.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//
// class decleration
//

using namespace std;
using namespace edm;

class CTCVWeightedVariable {
public:
  CTCVWeightedVariable( string name , string title , int nBins , double min , double max ){
    Name = name;
    edm::Service<TFileService> fs;
    Directory = fs->mkdir( name ) ; 
    for (uint i = 0 ; i < 70 ; i++){
      Histos.push_back( Directory.make< TH1D >( ("ctcv_"+to_string(i)).c_str() , (title + "," + to_string(i)).c_str() , nBins , min, max ) );
    }
  };

  void Fill( double value , std::vector<double> weights){
    Histos[0]->Fill( value );
    for( uint i = 0 ; i < weights.size() ; i++)
      Histos[i+1]->Fill( value , weights[i] );
  };

  void Write(){
    Directory.make< TCanvas >( ("Canvas_"+Name).c_str() );
    map< int , int > colors ;
    colors[0] = 2;
    colors[1] = 4;
    colors[12] = 8;
    /*for( auto h : Histos )
      if( h->Integral() > 0 )
      h->DrawNormalized();*/
    Histos[0]->SetTitle("ct/cv = -1");
    Histos[0]->SetLineColor( colors[0] );
    Histos[0]->SetMarkerColor( colors[0] );
    Histos[0]->DrawNormalized();

    Histos[1]->SetTitle("ct/cv = -3");
    Histos[1]->SetLineColor( colors[1] );
    Histos[1]->SetMarkerColor( colors[1] );
    Histos[1]->DrawNormalized("same");

    Histos[12]->SetTitle("ct/cv = 1");
    Histos[12]->SetLineColor( colors[12] );
    Histos[12]->SetMarkerColor( colors[12] );
    Histos[12]->DrawNormalized("same");
  }

  TFileDirectory Directory;
  vector< TH1* > Histos ;
  string Name;
};

class THHardProcess : public edm::EDFilter {
 public:
  map< string , CTCVWeightedVariable* > CTCVWeightedVariables;
  explicit THHardProcess(const edm::ParameterSet&);
  ~THHardProcess();

  virtual void endJob() override;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
 private:
  // ----------member data ---------------------------      
  EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
  EDGetTokenT<vector<pat::PackedGenParticle> > packedGenParticleToken_;
  EDGetTokenT<GenEventInfoProduct> genInfoToken_;

  bool usePacked_;
  bool debug_;
  edm::EDGetTokenT< LHEEventProduct > token_lhe;

  int npass,nfail,ntot;
  float wpass,wfail,wtot;
  virtual void beginJob() override;

  const reco::GenParticle* GetMother( const reco::GenParticle* gen ){
    if( gen == NULL )
      return NULL;
    if( gen->numberOfMothers() == 0 )
      return NULL;
    
    const reco::GenParticle* mom = dynamic_cast<const reco::GenParticle*>(gen->mother());

    while( true ){
      if( mom == NULL ){
	cout << "ERROR : mom is not convertible to GenParticle*" << endl;
	return NULL;
      }
      if( mom->pdgId() != gen->pdgId()  ){
	return mom;
      }
      if( mom->numberOfMothers() == 0 ){
	return NULL;
      }
      mom = dynamic_cast<const reco::GenParticle*>( mom->mother() );
    }
  }

  const reco::GenParticle* GetStableParticle( const reco::GenParticle* input ){
    if( input == NULL )
      return NULL;
    if( input->isLastCopy() )
      return input ;

    for( unsigned int i =0 ; i < input->numberOfDaughters() ; i++){
      const reco::GenParticle* gen = dynamic_cast< const reco::GenParticle*>( input->daughter(i) );
      if( gen->pdgId() == input->pdgId() )
	return GetStableParticle( gen );
    }
    return NULL;
  }
};
#endif
