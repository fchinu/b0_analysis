#include <TFile.h>
#include <TNtuple.h>
#include <TRandom3.h>
#include <TDatabasePDG.h>
#include <TROOT.h>
#include <TMath.h>
#include <TF1.h>

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include "Pythia8/Pythia.h"

//__________________________________________________________________________________________________
void simulateBtoDKa(int pdgB = 511, int nDecays = 1000000, int seed = 42) {

  gROOT->SetBatch(true);

  // init pythia object
  Pythia8::Pythia pythia;
  pythia.readString("SoftQCD:inelastic = on");

  pythia.readString("511:onMode = off");
  pythia.readString("521:onMode = off");
  pythia.readString("531:onMode = off");
  pythia.readString("5122:onMode = off");
  pythia.readString("411:onMode = off");
  pythia.readString("421:onMode = off");
  pythia.readString("431:onMode = off");
  pythia.readString("4122:onMode = off");
  pythia.readString("511:onIfMatch = -411 321");
  pythia.readString("411:onIfMatch = 211 -321 211");
  pythia.readString("Tune:pp = 14");
  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed %d", seed));
  pythia.init();

  gRandom->SetSeed(seed);

  std::map<int, int> pdgC = {
    {511, 411},
    {521, 421},
    {531, 431},
    {5122, 4122}
  };

  TFile* fonllFile = TFile::Open("simulations/b_to_dka/fonll.root");
  TH1D* hFonll = (TH1D*)fonllFile->Get("h_fonll");
  hFonll->SetDirectory(0);
  fonllFile->Close();

  TFile* massFile = TFile::Open("fit/outputs/dummy_template_DKa/B0_mass23_24_full_dataset.root");
  TH1D* hWidth = (TH1D*)massFile->Get("h_sigmas");
  hWidth->SetDirectory(0);
  massFile->Close();

  auto tree = new TNtuple("treeBtoDKa", "treeBtoDKa", "mB:mBSmeared:mBcorr:mBcorrSmeared:pTB:pTKa:selected");

  //__________________________________________________________
  // perform the simulation
  float massB = TDatabasePDG::Instance()->GetParticle(pdgB)->Mass();
  float massC = TDatabasePDG::Instance()->GetParticle(pdgC[pdgB])->Mass();
  float massPi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  float massKa = TDatabasePDG::Instance()->GetParticle(321)->Mass();
  for (auto iDecay{0}; iDecay < nDecays; ++iDecay) {
    auto ptB = hFonll->GetRandom(gRandom);
    auto yB = gRandom->Uniform(-1., 1.);
    auto phiB = gRandom->Rndm() * 2 * TMath::Pi();
    auto pxB = ptB * TMath::Cos(phiB);
    auto pyB = ptB * TMath::Sin(phiB);
    auto mtB = TMath::Sqrt(massB * massB + ptB * ptB);
    auto pzB = mtB * TMath::SinH(yB);
    auto pB = TMath::Sqrt(ptB * ptB + pzB * pzB);
    auto eB = TMath::Sqrt(massB * massB + pB * pB);

    // Bhad
    Pythia8::Particle Bhad;
    Bhad.id(pdgB);
    Bhad.status(81);
    Bhad.m(massB);
    Bhad.xProd(0.);
    Bhad.yProd(0.);
    Bhad.zProd(0.);
    Bhad.tProd(0.);
    Bhad.e(eB);
    Bhad.px(pxB);
    Bhad.py(pyB);
    Bhad.pz(pzB);
    Bhad.tau(0.f);

    pythia.event.reset();
    pythia.event.append(Bhad);
    auto idPart = pythia.event[1].id();
    pythia.particleData.mayDecay(idPart, true);
    pythia.moreDecays();

    ROOT::Math::PxPyPzMVector pVecC, pVecKa, pVecKaCorr;
    double massBSmeared{0.f}, massBCorrSmeared{0.f};
    int pdgLept{0};
    bool selected{true};
    for (int iPart{1}; iPart < pythia.event.size(); ++iPart)
    {
      auto absPdg = std::abs(pythia.event[iPart].id());
      vector<int> mothers = pythia.event[iPart].motherList();
      bool fromD = false;
      for (auto iMother : mothers) {
        if (std::abs(pythia.event[iMother].id()) == pdgC[pdgB]) {
          fromD = true;
          break;
        }
      }
      if (absPdg == pdgC[pdgB]) {
        pVecC = ROOT::Math::PxPyPzMVector(pythia.event[iPart].px(), pythia.event[iPart].py(), pythia.event[iPart].pz(), massC);
      } else if (absPdg == 321 && !fromD) {
        pVecKa = ROOT::Math::PxPyPzMVector(pythia.event[iPart].px(), pythia.event[iPart].py(), pythia.event[iPart].pz(), massPi); // we assume it to be a pion as in the analysis
        pVecKaCorr = ROOT::Math::PxPyPzMVector(pythia.event[iPart].px(), pythia.event[iPart].py(), pythia.event[iPart].pz(), TDatabasePDG::Instance()->GetParticle(absPdg)->Mass());
      }
      if ((absPdg == 321 || absPdg == 211) && pythia.event[iPart].pT() < 1.0) {
        selected = false;
      }
    }
    ROOT::Math::PxPyPzMVector pVecB = pVecC + pVecKa;
    ROOT::Math::PxPyPzMVector pVecBcorr = pVecC + pVecKaCorr;
    if (pVecB.Pt() > 23.5 || pVecB.Pt() < 1.0) continue;
    massBSmeared = gRandom->Gaus(pVecB.M(), hWidth->GetBinContent(hWidth->FindBin(pVecB.Pt())));
    massBCorrSmeared = gRandom->Gaus(pVecBcorr.M(), hWidth->GetBinContent(hWidth->FindBin(pVecB.Pt())));
    tree->Fill(pVecB.M(), massBSmeared, pVecBcorr.M(), massBCorrSmeared, pVecBcorr.Pt(), pVecKaCorr.Pt(), selected);
  }

  TFile ouputFile(Form("simulations/b_to_dka/MDPi_%d_from_DKa_decay_with_pt_cuts.root", pdgB), "recreate");
  tree->Write();
  hFonll->Write();
  ouputFile.Close();
}