#include <iostream>
#include <cmath>

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif


double jetFunction(double pt) {
    return 0.85 * tanh(0.0025 * pt) * (25.0 / (1 + 0.063 * pt));
}


void plotUp() {
  gSystem->Load("libDelphes");
  
    // Create chain of root trees for background file
    TChain backgroundChain("Delphes");
    backgroundChain.Add("finalb.root");

    // Create chain of root trees for signal file
    TChain signalChain("Delphes");
    signalChain.Add("/home/phodiso/Downloads/tag_1_delphes_events_1_2.root");

    // Create object of class ExRootTreeReader for background
    ExRootTreeReader *backgroundReader = new ExRootTreeReader(&backgroundChain);

    // Create object of class ExRootTreeReader for signal
    ExRootTreeReader *signalReader = new ExRootTreeReader(&signalChain);

    // Get pointers to branches used in this analysis for background
    TClonesArray *backgroundJet = backgroundReader->UseBranch("Jet");
    TClonesArray *backgroundPhoton = backgroundReader->UseBranch("Photon");
    // Get pointers to branches used in this analysis for signal
    TClonesArray *signalJet = signalReader->UseBranch("Jet");
    TClonesArray *signalPhoton = signalReader->UseBranch("Photon");

  //const double backgroundCrossSection = 4.081;
  const double backgroundCrossSection = 0.53;
  const double signalCrossSection = 0.004;

  // Integrated luminosity (in pb^-1)
  const double luminosity = 139000; // Example luminosity, adjust as per your data

  // Calculate event weights
  double backgroundWeight = 7.7 * backgroundCrossSection * luminosity / backgroundChain.GetEntries();
  double signalWeight = signalCrossSection * luminosity / signalChain.GetEntries();

  // Create histogram
  TH1F *histInvariantMassDiphotonBackground = new TH1F("InvariantMassDiphoton", "Invariant Mass of Diphoton; m_{#gamma#gamma} [GeV]; Events/2.5GeV", 22, 105, 160);
  TH1F *histInvariantMassDiphotonSignal = new TH1F("InvariantMassDiphoton", "Invariant Mass of Diphoton; m#gamma#gamma [GeV]", 22, 105, 160);
  TH1F *histJetFunctionBackground = new TH1F("JetFunctionBackground", "Jet Function (Background); pT [GeV]; Value", 100, 0, 1000);
  TH1F *histJetFunctionSignal = new TH1F("JetFunctionSignal", "Jet Function (Signal); pT [GeV]; Value", 100, 0, 1000);
  TGraph *graphJetFunctionBackground = new TGraph();
  TGraph *graphJetFunctionSignal = new TGraph();  

  // Initialize counters for cut flow
  double NumberWeightedEvents = 0.0;
  double NumberWeightedEventsignal = 0.0;
  double weightedEventsInRange = 0.0;
  double weightedEventsInRanges = 0.0;
  int passJetPtCut = 0;
  int Bee = 0;
  int Bees = 0;
  int bTaggedJetss = 0;
  int weights = 0;
  int totalEvents = 0;
  int totalEventsignal = 0;
  int passPTandEtaCuts = 0;
  int passPTandEtaCut = 0;
  int passETcut = 0;
  int passETcuts = 0;
  int passMassRangeCut = 0;
  int passMassRangeCuts = 0;
  int passPtRatioCut1 = 0;
  int passPtRatioCut1s = 0;
  int passPtRatioCut2 = 0;
  int passPtRatioCut2s = 0;  
  int passJetEtaCut = 0;
  int passJetEtaCuts = 0;
  int eventsInRange = 0;
  int eventsInRanges = 0;
  int passBTaggedJet = 0;
  int passBTaggedJets = 0;
  double significance = 0.0;
  double scale = 0.0;
  double sumJetPt = 0.0;
  int countJets = 0;

  // Loop over events for background
  Long64_t backgroundEntries = backgroundReader->GetEntries();
  for (Long64_t entry = 0; entry < backgroundEntries; ++entry) {
    backgroundReader->ReadEntry(entry);

    // Access objects from branches
    TClonesArray *jets = backgroundJet;
    TClonesArray *photons = backgroundPhoton;


    // Apply event weights
    double weight = backgroundWeight;
    NumberWeightedEvents += weight;   
    totalEvents++;

    // Select photons with pT > 25 GeV and |eta| < 2.37
    bool passPTandEtaCutEvent = true;
    for (int i = 0; i < photons->GetEntries(); ++i) {
      Photon *photon = (Photon *)photons->At(i);
      if (photon->PT <= 25 || fabs(photon->Eta) >= 2.37 || (fabs(photon->Eta) > 1.37 && fabs(photon->Eta) < 1.52)) {
        passPTandEtaCutEvent = false;
        break;
      }
    }
    if (!passPTandEtaCutEvent) {
      continue; // Skip events failing photon pT and eta cut
    }
    passPTandEtaCut++;

    // Apply other cuts as before...
    // Cut 1: Select events with exactly 2 photons
    if (photons->GetEntries() != 2) {
        continue; // Skip events with incorrect number of photons
    }
    passETcut++;

    // Cut 2: Apply ET cut
    Photon *photon1 = (Photon *)photons->At(0);
    Photon *photon2 = (Photon *)photons->At(1);
    if (photon1->PT <= 35) {
        continue; // Skip events failing ET cut
    }
    passMassRangeCut++;

    // Cut 3: Apply invariant mass range cut
    TLorentzVector photon1Vector, photon2Vector;
    photon1Vector.SetPtEtaPhiM(photon1->PT, photon1->Eta, photon1->Phi, 0);
    photon2Vector.SetPtEtaPhiM(photon2->PT, photon2->Eta, photon2->Phi, 0);
    TLorentzVector diphotonVector = photon1Vector + photon2Vector;
    if (diphotonVector.M() <= 105 || diphotonVector.M() >= 160) {
        continue; // Skip events outside mass range
    }
    passPtRatioCut1++;

    // Cut 4: Apply pT/mγγ ratio cut for photon 1
    if (photon1->PT / diphotonVector.M() <= 0.35) {
        continue; // Skip events failing pT/mγγ ratio cut for photon 1
    }
    passPtRatioCut2++;

    // Cut 5: Apply pT/mγγ ratio cut for photon 2
    if (photon2->PT / diphotonVector.M() <= 0.25) {
        continue; // Skip events failing pT/mγγ ratio cut for photon 2
    }
    //passJetPtCut++;
    
   
        
     for (int i = 0; i < jets->GetEntries(); ++i) {
    Jet *jet = (Jet *)jets->At(i);
    double pt = jet->PT;

    // Check if the jet passes your additional cuts (if any)
    // For example, you might want to check the jet's pT or eta here
        // Additional cuts: Jet Eta and b-tagging
    bool passJetEtaCutEvent = false;
    for (int i = 0; i < jets->GetEntries(); ++i) {
        Jet *jet = (Jet *)jets->At(i);
        if (fabs(jet->Eta) < 2.5) {
            passJetEtaCutEvent = true;
            passJetEtaCut++;
            break;
        }
    }

    if (!passJetEtaCutEvent) {
        continue; // Skip events failing jet PT cut
    }
    for (int i = 0; i < jets->GetEntries(); ++i) {
            Jet *jet = (Jet *)jets->At(i);
            double pt = jet->PT;
            double value = jetFunction(pt);
            histJetFunctionBackground->Fill(pt, value * backgroundWeight);
            //graphJetFunctionBackground->SetPoint(graphJetFunctionBackground->GetN(), pt, value);
            int pointIndex = graphJetFunctionBackground->GetN();
            graphJetFunctionBackground->SetPoint(pointIndex, pt, jetFunction(pt)*backgroundWeight);
             sumJetPt += pt;
             countJets++;
        }
    // If the jet passes the cuts, add its pT to the sum
    //sumJetPt += pt;
    // Increment the count of jets that passed the cuts
    //countJets++;
}

std::vector<Jet*> bTaggedJets;
        for (int i = 0; i < jets->GetEntries(); ++i) {
            Jet *jet = (Jet *)jets->At(i);
            if (jet->BTag) {
                bTaggedJets.push_back(jet);
            }
        }

        if (bTaggedJets.size() < 2) {
            continue; // Skip events with incorrect number of b-tagged jets
        }

        //passBTaggedJet++;

        // Calculate di-btag jet invariant mass
        TLorentzVector jet1Vector, jet2Vector;
        jet1Vector.SetPtEtaPhiM(bTaggedJets[0]->PT, bTaggedJets[0]->Eta, bTaggedJets[0]->Phi, bTaggedJets[0]->Mass);
        jet2Vector.SetPtEtaPhiM(bTaggedJets[1]->PT, bTaggedJets[1]->Eta, bTaggedJets[1]->Phi, bTaggedJets[1]->Mass);
        TLorentzVector diBtagJetVector = jet1Vector + jet2Vector;
// Calculate the average jet pT
//double avgJetPt = sumJetPt / countJets;   

    // Additional cuts: Jet Eta and b-tagging
    Bee++;
    //if (diBtagJetVector.M() <= 63.7918 || diBtagJetVector.M() >= 127.2082) {
      //continue;
          //eventsInRanges++;
          //histInvariantMassDibjetSignal->Fill(dibjetVector.M(), weight);
          //weightedEventsInRanges += signalWeight;
      //}

    //int bTaggedJets = 0;
    //for (int i = 0; i < jets->GetEntries(); ++i) {
        //Jet *jet = (Jet *)jets->At(i);
        //if (jet->BTag) {
        //    bTaggedJets++;
      //  }
    //}
    //bTaggedJetss++;
    //if (bTaggedJets != 2) {
      //  continue; // Skip events failing b-tagged jet cut
    //}

    passBTaggedJet++;
   // scale = (double)passBTaggedJet/2.5;
    histInvariantMassDiphotonBackground->Fill(diphotonVector.M(), weight);
    
    // Count events within invariant mass range (147.524 - 155.476)
    //TLorentzVector diphotonVector = photon1Vector + photon2Vector;
    if (diphotonVector.M() >= 147.524 && diphotonVector.M() <= 155.476) {
      eventsInRange++;
      //histInvariantMassDiphotonBackground->Fill(diphotonVector.M(), weight);
      weightedEventsInRange += weight;
    }

    // Fill histogram with weighted events passing all cuts

  }
  
  // Loop over events for signal
  Long64_t signalEntries = signalReader->GetEntries();
  for (Long64_t entry = 0; entry < signalEntries; ++entry) {
    signalReader->ReadEntry(entry);

    // Access objects from branches
    TClonesArray *jets = signalJet;
    TClonesArray *photons = signalPhoton;

    totalEventsignal++;

    // Apply event weights
    double weight = signalWeight;
    NumberWeightedEventsignal += weight;   

    // Select photons with pT > 25 GeV and |eta| < 2.37
    bool passPTandEtaCutEvents = true;
    for (int i = 0; i < photons->GetEntries(); ++i) {
      Photon *photon = (Photon *)photons->At(i);
      if (photon->PT <= 25 || fabs(photon->Eta) >= 2.37 || (fabs(photon->Eta) > 1.37 && fabs(photon->Eta) < 1.52)) {
        passPTandEtaCutEvents = false;
        break;
      }
    }
    if (!passPTandEtaCutEvents) {
      continue; // Skip events failing photon pT and eta cut
    }
    passPTandEtaCuts++;

    // Apply other cuts as before...
    // Cut 1: Select events with exactly 2 photons
    if (photons->GetEntries() != 2) {
        continue; // Skip events with incorrect number of photons
    }
    passETcuts++;

    // Cut 2: Apply ET cut
    Photon *photon1 = (Photon *)photons->At(0);
    Photon *photon2 = (Photon *)photons->At(1);
    if (photon1->PT <= 35) {
        continue; // Skip events failing ET cut
    }
    passMassRangeCuts++;

    // Cut 3: Apply invariant mass range cut
    TLorentzVector photon1Vector, photon2Vector;
    photon1Vector.SetPtEtaPhiM(photon1->PT, photon1->Eta, photon1->Phi, 0);
    photon2Vector.SetPtEtaPhiM(photon2->PT, photon2->Eta, photon2->Phi, 0);
    TLorentzVector diphotonVector = photon1Vector + photon2Vector;
    if (diphotonVector.M() <= 105 || diphotonVector.M() >= 160) {
        continue; // Skip events outside mass range
    }
    passPtRatioCut1s++;

    // Cut 4: Apply pT/mγγ ratio cut for photon 1
    if (photon1->PT / diphotonVector.M() <= 0.35) {
        continue; // Skip events failing pT/mγγ ratio cut for photon 1
    }
    passPtRatioCut2s++;

    // Cut 5: Apply pT/mγγ ratio cut for photon 2
    if (photon2->PT / diphotonVector.M() <= 0.25) {
        continue; // Skip events failing pT/mγγ ratio cut for photon 2
    }
    passJetPtCut++;

    // Additional cuts: Jet Eta and b-tagging
    bool passJetEtaCutEvents = false;
    for (int i = 0; i < jets->GetEntries(); ++i) {
        Jet *jet = (Jet *)jets->At(i);
        if (fabs(jet->Eta) < 2.5) {
            passJetEtaCutEvents = true;
            passJetEtaCuts++;
            break;
        }
    }

    if (!passJetEtaCutEvents) {
        continue; // Skip events failing jet PT cut
    }
    
    std::vector<Jet*> bTaggedJets;
        for (int i = 0; i < jets->GetEntries(); ++i) {
            Jet *jet = (Jet *)jets->At(i);
            if (jet->BTag) {
                bTaggedJets.push_back(jet);
            }
        }

        if (bTaggedJets.size() < 2) {
            continue; // Skip events with incorrect number of b-tagged jets
        }

        //passBTaggedJets++;

        // Calculate di-btag jet invariant mass
        TLorentzVector jet1Vector, jet2Vector;
        jet1Vector.SetPtEtaPhiM(bTaggedJets[0]->PT, bTaggedJets[0]->Eta, bTaggedJets[0]->Phi, bTaggedJets[0]->Mass);
        jet2Vector.SetPtEtaPhiM(bTaggedJets[1]->PT, bTaggedJets[1]->Eta, bTaggedJets[1]->Phi, bTaggedJets[1]->Mass);
        TLorentzVector diBtagJetVector = jet1Vector + jet2Vector;
        Bees++;

//if (diBtagJetVector.M() <= 63.7918 || diBtagJetVector.M() >= 127.2082) {
      //continue;
          //eventsInRanges++;
          //histInvariantMassDibjetSignal->Fill(dibjetVector.M(), weight);
          //weightedEventsInRanges += signalWeight;
     // }
    //int bTaggedJets = 0;
    //for (int i = 0; i < jets->GetEntries(); ++i) {
        //Jet *jet = (Jet *)jets->At(i);
        //if (jet->BTag) {
        //    bTaggedJets++;
      //  }
    //}
    //if (bTaggedJets != 2) {
      //  continue; // Skip events failing b-tagged jet cut
    //}

    passBTaggedJets++;

    histInvariantMassDiphotonSignal->Fill(diphotonVector.M(), weight);
    
    // Count events within invariant mass range (147.524 - 155.476)
    //TLorentzVector diphotonVector = photon1Vector + photon2Vector;
    if (diphotonVector.M() >= 147.524 && diphotonVector.M() <= 155.476) {
      eventsInRanges++;
      //histInvariantMassDiphotonSignal->Fill(diphotonVector.M(), weight);
      weightedEventsInRanges += weight;
    }

    // Fill histogram with weighted events passing all cuts

  }

  double efficiencyETcut = (double)passETcut / totalEvents;
  double efficiencyMassRangeCut = (double)passMassRangeCut / totalEvents;
  double efficiencyPtRatioCut1 = (double)passPtRatioCut1 / totalEvents;
  double efficiencyPtRatioCut2 = (double)passPtRatioCut2 / totalEvents;
  double efficiencyJetEtaCut = (double)passJetEtaCut / totalEvents;
  double efficiencyBTaggedJet = (double)passBTaggedJet / totalEvents;
  double beefs = (double)Bee /totalEvents;
  
  double efficiencyETcuts = (double)passETcuts / totalEventsignal;
  double efficiencyMassRangeCuts = (double)passMassRangeCuts / totalEventsignal;
  double efficiencyPtRatioCut1s = (double)passPtRatioCut1s / totalEventsignal;
  double efficiencyPtRatioCut2s = (double)passPtRatioCut2s / totalEventsignal;
  double efficiencyJetEtaCuts = (double)passJetEtaCuts / totalEventsignal;
  double efficiencyBTaggedJets = (double)eventsInRanges /totalEventsignal;
  double beefss = (double)Bees /totalEventsignal;

  // Print cut flow
  cout << "Background: " << endl;
  cout << "Total events: " << totalEvents << endl;
  cout << "Weights: " << backgroundWeight << endl;
  cout << "Number of weighted events: " << NumberWeightedEvents << endl;
  cout << "Photon pT > 25 GeV and |eta| < 2.37: " << passPTandEtaCut << " Efficiency: " << (double)passPTandEtaCut / totalEvents << endl;
  // Print other cuts...
  cout << "ET(γ1) > 35 GeV: " << passETcut << " Efficiency: " << efficiencyETcut << endl;
  cout << "105 GeV < mγγ < 160 GeV: " << passMassRangeCut << " Efficiency: " << efficiencyMassRangeCut << endl;
  cout << "pT/mγγ > 0.35 for γ1: " << passPtRatioCut1 << " Efficiency: " << efficiencyPtRatioCut1 << endl;
  cout << "pT/mγγ > 0.25 for γ2: " << passPtRatioCut2 << " Efficiency: " << efficiencyPtRatioCut2 << endl;
  //cout << "Jet |Eta| < 2.5: " << passJetEtaCut << " Efficiency: " << efficiencyJetEtaCut << endl;
  cout << "B-jet: " << Bee << " Efficiency: " << beefs << endl;
  cout << "At least 2 b-tagged jets and 2 photons: " << passBTaggedJet << " Efficiency: " << efficiencyBTaggedJet << endl;
  // Print events within invariant mass range and expected number of signal events
  cout << "Events within invariant mass range (147.524 - 155.476 GeV): " << eventsInRange << endl;
  cout << "Weighted events within invariant mass range: " << weightedEventsInRange << endl;
  double expectedBackgroundEvents =backgroundCrossSection * luminosity * efficiencyBTaggedJet;
  cout << "Expected number of background events: " << expectedBackgroundEvents << endl;
  double avgJetPt = sumJetPt / countJets;
  cout << "Average pT: " << avgJetPt << endl;

  cout << "Signal: " << endl;  
  cout << "Total events: " << totalEventsignal << endl;
  cout << "Weights: " << signalWeight << endl;
  cout << "Number of weighted events: " << NumberWeightedEventsignal << endl;
  cout << "Photon pT > 25 GeV and |eta| < 2.37: " << passPTandEtaCuts << " Efficiency: " << (double)passPTandEtaCuts / totalEventsignal << endl;
  // Print other cuts...
  cout << "ET(γ1) > 35 GeV: " << passETcuts << " Efficiency: " << efficiencyETcuts << endl;
  cout << "105 GeV < mγγ < 160 GeV: " << passMassRangeCuts << " Efficiency: " << efficiencyMassRangeCuts << endl;
  cout << "pT/mγγ > 0.35 for γ1: " << passPtRatioCut1s << " Efficiency: " << efficiencyPtRatioCut1s << endl;
  cout << "pT/mγγ > 0.25 for γ2: " << passPtRatioCut2s << " Efficiency: " << efficiencyPtRatioCut2s << endl;
  //cout << "Jet |Eta| < 2.5: " << passJetEtaCut << " Efficiency: " << efficiencyJetEtaCut << endl;
  cout << "B-jet selection: " << Bees << " Efficiency: " << beefss << endl;
  cout << "At least 2 b-tagged jets and 2 photons: " << passBTaggedJets << " Efficiency: " << efficiencyBTaggedJets << endl;
  // Print events within invariant mass range and expected number of signal events
  cout << "Events within invariant mass range (147.524 - 155.476 GeV): " << eventsInRanges << endl;
  cout << "Weighted events within invariant mass range: " << weightedEventsInRanges << endl;
  double expectedSignalEvents =signalCrossSection * luminosity * efficiencyBTaggedJets;
  cout << "Expected number of signal events: " << expectedSignalEvents << endl;
  significance = weightedEventsInRanges / sqrt(weightedEventsInRange);
  cout << "Significance: " << significance << endl;

  // Create canvas and draw histogram
  TCanvas *canvas = new TCanvas("canvas", "Histogram", 800, 600);
  histInvariantMassDiphotonBackground->SetFillColor(kCyan);
  histInvariantMassDiphotonSignal->SetFillColor(kPink);
  histInvariantMassDiphotonBackground->Draw("HIST");
  histInvariantMassDiphotonSignal->Draw("HISTSAME");
  
  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(histInvariantMassDiphotonSignal, "g g > H > #gamma #gamma b b~ (Signal)", "f");
  legend->AddEntry(histInvariantMassDiphotonBackground, "p p > #gamma #gamma b b~ (Background)", "f");
  //legend->AddEntry(histInvariantMassDiphotonBackground, Form("Background (%.2f pb)", backgroundCrossSection), "f");
  //legend->AddEntry(nullptr, Form("Luminosity: %.2f pb^-1", luminosity), "");
  //legend->AddEntry(nullptr, Form("Weight: %.5f", backgroundWeight), "");
  //legend->AddEntry(nullptr, Form("Events passed: %d", passBTaggedJet), "");
  //legend->AddEntry((const char*)nullptr, Form("Background (%.2f pb)", backgroundCrossSection), "");
  //legend->AddEntry((const char*)nullptr, Form("Luminosity: %.2f pb^-1", luminosity), "");
  //legend->AddEntry((const char*)nullptr, Form("Weight: %.5f", backgroundWeight), "");
  legend->SetBorderSize(0);
  legend->Draw();
  
    //TCanvas *canvas = new TCanvas("canvas", "Jet Function", 800, 600);
    //histJetFunctionBackground->SetFillColor(kCyan);
    //histJetFunctionSignal->SetLineColor(kRed);
    //histJetFunctionBackground->Draw();
    //histJetFunctionSignal->Draw("SAME");

    // Add legend
    //TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    //legend->AddEntry(histJetFunctionBackground, "Background", "f");
    //legend->AddEntry(histJetFunctionSignal, "Signal", "l");
    //legend->Draw();
    
    //TCanvas *canvas = new TCanvas("canvas", "Jet Function", 800, 600);
    //graphJetFunctionBackground->SetLineColor(kCyan);
    //graphJetFunctionSignal->SetLineColor(kRed);
    //graphJetFunctionBackground->Draw("AP");
    //graphJetFunctionSignal->Draw("L SAME");
}

