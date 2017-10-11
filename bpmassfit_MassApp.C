#include "TTree.h"
using namespace RooFit ;
using namespace std;

void bpmassfit_MassApp() 
{            
  
  float JpsiPhiKMassmin=5.15, JpsiPhiKMassmax=5.45;


  RooRealVar varJpsiPhiKMass("varJpsiPhiKMass", "m(J/#\psi#phiK^{+}) GeV", JpsiPhiKMassmin,  JpsiPhiKMassmax) ;
//  RooRealVar varJpsiPhiMass("varJpsiPhiMass", "m(J/#\psi#phi) GeV", 4.1047,  4.9047) ; 
 RooRealVar varJpsiPhiMass("varJpsiPhiMass", "m(J/#psi#phi) GeV", 0.,  99.) ;
  RooRealVar varPhiKMass("varPhiKMass", " ", 0.,  4.) ;
  RooRealVar varKaon1Pt("varKaon1Pt","Kaon1Pt", .7, 99.);  //1.5
  RooRealVar varKaon2Pt("varKaon2Pt","Kaon2Pt", .7, 99.);
  RooRealVar varKaon3Pt("varKaon3Pt","Kaon3Pt", .7, 99.);
  RooRealVar varBvtx("varBvtx", "Bvtx", 0.05, 999.); //0.2
  RooRealVar varPhiMassWindow("varPhiMassWindow", "PhiMassWindow", 0.0 , 99.); //0.006
  RooRealVar varBLTSig("varBLTSig", "BLTSig", 1., 999.); //8.
  RooRealVar varJpsiMass("varJpsiMass", "JpsiMassWindow", 0.0, 99.); //0.06
  RooRealVar varDeltaR1("varDeltaR1", "DeltaR1", 0., 1.5); //do not change.
  RooRealVar varDeltaR2("varDeltaR2", "DeltaR2", 0., 1.5);
  RooRealVar varDeltaR3("varDeltaR3", "DeltaR3", 0., 1.5);
  RooRealVar varJPsiFL("varJPsiFL","JPsiFLsig", 8., 99.); //6 
  RooRealVar varBpt("varBpt", "varBpt", 0., 99.);
  RooRealVar varcosEtaPV("varcosEtaPV", "cosEtaPV", 0., 1.01);
  RooRealVar varJPsiPt("JPsiPt", "JPsiPt", 8., 999. );
  RooRealVar varBeta("varBeta", "Beta", 0., 2.2);
  RooRealVar varJPsiEta("varJPsiEta", "JPsiEta", 0., 2.2);
  RooRealVar varPsiTwoSMass("varPsiTwoSMass", "varPsiTwoSMass", 0.0, 99.);
  RooRealVar varJPsicosHel("varJPsicosHel", "JPsiHelicity", -1., 1.);
  RooRealVar varPhicosHel("varPhicosHel", "PhiHelicity", -1., 1.);

  RooRealVar   Bp_peak  ("Bp_peak"    ,"Bp_peak"  ,5.2789, 5.26, 5.3 ); //org
  RooRealVar   Bp_width ("Bp_width"    ,"Bp_width"  ,0.0086181, 0.001, 0.1 ); //org
  
  RooAbsPdf* pdfBp = new RooGaussian( "pdfBp", "pdfBp",  varJpsiPhiKMass, Bp_peak ,Bp_width ) ;
  RooRealVar nSigBp ("nSigBp", "Number of signal 1 candidates ", 100,  -10.0, 1000000.0); 

  Bp_peak.setConstant();
  Bp_width.setConstant();

  RooRealVar c1("c1", "c1", -4.66666e-02,  -10.0, 10.0); //org
  RooRealVar c2("c2", "c2", -3.54591e-01,  -10., 1.0); //org
  RooAbsPdf *  BkgPolPdf = new RooChebychev("BkgPolPdf","BkgPolPdf",varJpsiPhiKMass,RooArgSet(c1,c2));
  
  RooRealVar nBckPol("nBckPol", "Number of signal 2 candidates ", 10000,  0.0, 9000000.0);

  RooExtendPdf *  extendpdfSigBp = new RooExtendPdf("extendpdfSigBp","Signal 1 PDF",*pdfBp, nSigBp);
  RooExtendPdf *  extendpdfBkgPol = new RooExtendPdf("extendpdfBkgPol","Signal 1 PDF",*BkgPolPdf, nBckPol);

  RooAddPdf mytotalPdf("mytotPdf", "mytotPdf", RooArgList(*extendpdfSigBp, *extendpdfBkgPol), RooArgList(nSigBp,nBckPol) ) ;

    RooArgSet s(varJpsiPhiKMass,  varKaon1Pt, varKaon2Pt, varKaon3Pt, varBLTSig, varcosEtaPV, varBvtx, varJpsiMass);
  s.add(RooArgSet( varPhiMassWindow, varBpt, varDeltaR1, varDeltaR2, varDeltaR3, varPsiTwoSMass));
 //   RooArgSet s(varJpsiPhiKMass, varKaon1Pt, varKaon2Pt, varKaon3Pt, varBLTSig, varBvtx, varJpsiMass);
//  s.add(RooArgSet( varPhiMassWindow, varPsiTwoSMass));
  RooDataSet *dataSet = (RooDataSet::read("B-kp1-kp2-kp3-BLT-PA-Bvtx-JMass-PhiM-Bpt-D1-D2-D3-Psi2SM_AllforTMVAB.txt", s, "Q")); 
// RooDataSet *dataSet = (RooDataSet::read("B-kp1-kp2-kp3-BLT-Bvtx-JMass-PhiM-Psi2SM_AllforTMVAB.txt", s, "Q"));

/* varPsiTwoSMass.setRange("selectionA", 0., 3.606);
  varPsiTwoSMass.setRange("selectionB", 3.766, 99.);
  RooDataSet * dataSet = (RooDataSet*)psiTwoSMassL->reduce(CutRange("selectionA"));
  RooDataSet * psiTwoSMassH = (RooDataSet*)psiTwoSMassL->reduce(CutRange("selectionB"));
  dataSet->append(*psiTwoSMassH);*/

  mytotalPdf.fitTo(*dataSet,"mer");  
 
  TCanvas * c=new TCanvas("c","c",800,600);
  RooPlot *frame = varJpsiPhiKMass.frame(60);
  dataSet->plotOn(frame);  
  mytotalPdf.plotOn(frame, Components(RooArgSet(*extendpdfBkgPol)),LineStyle(kDashed),LineColor(kBlue),Range(JpsiPhiKMassmin,JpsiPhiKMassmax) );
  mytotalPdf.paramOn(frame, Format("NE",FixedPrecision(5)));
  mytotalPdf.plotOn(frame);


  frame->SetTitle("");  
  frame->GetXaxis()->SetTitleSize(0.05);
  frame->GetXaxis()->SetLabelSize(0.045);
  frame->GetXaxis()->SetTitleOffset(0.85);
  frame->GetYaxis()->SetTitleSize(0.05);
  frame->GetYaxis()->SetTitle ("Candidates / 0.005 GeV");
  frame->GetYaxis()->SetTitleOffset(0.85);
  frame->SetMarkerColor(1);
  frame->SetMarkerSize(2.0);
  frame->Draw();

/*  double myBGausSigma = Bp_width.getVal();
  varJpsiPhiKMass.setRange("selection",JpsiPhiKMassmin,JpsiPhiKMassmax);
  RooAbsReal*  pdfBptotal=pdfBp->createIntegral(RooArgSet(varJpsiPhiKMass),"selection");
  RooAbsReal* BkgPolPdftotal=BkgPolPdf->createIntegral(RooArgSet(varJpsiPhiKMass),"selection");*/

   TFile f("forTMVApp.root","recreate");
    TTree* chain =  RooStats::GetAsTTree("TreeS","TreeS",*dataSet); //ds2
    chain->Fill();
    f.Write();
   f.Close();

}

