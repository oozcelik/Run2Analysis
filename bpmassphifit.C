#include "TColor.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TLatex.h"

using namespace std;
using namespace RooFit ;

void bpmassphifit(char FileName[128])
{

  float JpsiPhiKMassmin=5.15, JpsiPhiKMassmax=5.45;
  std::ofstream myoutfile("PhiM.txt");
   TH1F *varPhiMass = new TH1F ("PhiMass","m(#phi{+}) GeV", 6, 0.986,1.058); // 18

   varPhiMass-> GetXaxis()->SetTitle("PhiMass");
   varPhiMass->GetYaxis()->SetTitle("Candidates / 4 MeV");
   varPhiMass->Sumw2();

   TCanvas *c11 = new TCanvas("c1", "c1", 700, 700);
   c11->Divide(2, 2);
   TCanvas        *c33 = new TCanvas("c33", "c33", 700, 700);
   c33->Divide(2, 2);
   TCanvas        *c44 = new TCanvas("c44", "c44", 700, 700);
   c44->Divide(2, 2);
   TCanvas        *c55 = new TCanvas("c55", "c55", 700, 700);
   c55->Divide(2, 2);
   TCanvas        *c66 = new TCanvas("c66", "c66", 700, 700);
   c66->Divide(2, 1);

        TLatex         *t = new TLatex();
        t->SetNDC();
        t->SetTextSize(0.05);
        t->SetTextAlign(12);
        t->SetTextFont(42);
        t->SetTextColor(kRed);

       char DeltaM[28];
      for (int ii = 0; ii < 6; ii++) {

        float dmbinlow = 0.986 + (float) ii * 0.012; // 0.004
        float dmbinhigh = 0.986 + (float) (ii + 1) * 0.012; //0.004

        RooRealVar varJpsiPhiKMass("varJpsiPhiKMass", "m(J/#\psi#phiK^{+}) GeV", JpsiPhiKMassmin,  5.45);
        RooRealVar varphimass("varphimass", "MassOfPhi " , dmbinlow, dmbinhigh );
        
        sprintf(DeltaM, " %.2f GeV < #Deltam  < %.2f GeV", dmbinlow, dmbinhigh);

        RooRealVar   Bp_peak  ("Bp_peak"    ,"Bp_peak"  ,5.2789, 5.26, 5.3 ); //org
        RooRealVar   Bp_width ("Bp_width"    ,"Bp_width"  ,0.0086181, 0.001, 0.1 ); //org

        RooAbsPdf* pdfBp = new RooGaussian( "pdfBp", "pdfBp",  varJpsiPhiKMass, Bp_peak ,Bp_width ) ;
        RooRealVar nSigBp ("nSigBp", "Number of signal 1 candidates ", 100,  -10.0, 1000000.0);

        RooRealVar c1("c1", "c1", -4.66666e-02,  -10.0, 10.0); //org
        RooRealVar c2("c2", "c2", -3.54591e-01,  -10., 1.0); //org
        RooAbsPdf *  BkgPolPdf = new RooChebychev("BkgPolPdf","BkgPolPdf",varJpsiPhiKMass,RooArgSet(c1,c2));
        RooRealVar nBckPol("nBckPol", "Number of signal 2 candidates ", 10000,  0.0, 9000000.0);

        RooExtendPdf *  extendpdfSigBp = new RooExtendPdf("extendpdfSigBp","Signal 1 PDF",*pdfBp, nSigBp);
        RooExtendPdf *  extendpdfBkgPol = new RooExtendPdf("extendpdfBkgPol","Signal 1 PDF",*BkgPolPdf, nBckPol);

        RooAddPdf mytotalPdf("mytotPdf", "mytotPdf", RooArgList(*extendpdfSigBp, *extendpdfBkgPol), RooArgList(nSigBp,nBckPol) ) ;
        RooDataSet * dataSet= (RooDataSet::read(FileName,RooArgSet(varJpsiPhiKMass, varphimass),"Q"));
        gROOT->SetStyle("Plain");
        mytotalPdf.fitTo(*dataSet,"mer"); 
        
        if (ii < 4) c11->cd(ii + 1);
        else if (ii < 8)  c33->cd(ii + 1 - 4);
        else if (ii < 12) c44->cd(ii + 1 - 8);
        else if (ii < 16) c55->cd(ii + 1 - 12);
        else if (ii < 18) c66->cd(ii + 1 - 18);

        RooPlot *frame = varJpsiPhiKMass.frame(75);
        dataSet->plotOn(frame);
        mytotalPdf.plotOn(frame, Components(RooArgSet(*extendpdfBkgPol)),LineStyle(kDashed),LineColor(kBlue),Range(JpsiPhiKMassmin,JpsiPhiKMassmax) );
        mytotalPdf.plotOn(frame, Components(RooArgSet(*extendpdfSigBp)), LineStyle(kSolid), LineColor(kRed), Range(JpsiPhiKMassmin, JpsiPhiKMassmax));
        // mytotalPdf.plotOn(frame, Components(RooArgSet(*extendpdfBkgPol)),LineStyle(kDashed),LineColor(kBlue),Range(JpsiPhiKMassmin,JpsiPhiKMassmax) );
        // mytotalPdf.paramOn(frame, Format("NE",FixedPrecision(5)));
        mytotalPdf.plotOn(frame);

        frame->SetTitle("");
        frame->GetXaxis()->SetTitleSize(0.05);
        frame->GetXaxis()->SetLabelSize(0.045);
        frame->GetXaxis()->SetTitleOffset(0.85);
        frame->GetYaxis()->SetTitleSize(0.05);
        frame->GetYaxis()->SetTitle ("Candidates / 0.004 GeV");
        frame->GetYaxis()->SetTitleOffset(0.85);
        frame->SetMarkerColor(1);
        frame->SetMarkerSize(2.0);
        frame->Draw();

       myoutfile  << ii << " " << nSigBp.getVal() << " " << nSigBp.getError() << endl;
              //Const[ii] = frame->chiSquare(4);
       sprintf(DeltaM, " %2.0f GeV < #Deltam  < %2.0f GeV", dmbinlow, dmbinhigh);
       t->DrawLatex(0.2, 0.8, DeltaM);
       varPhiMass->SetBinContent(ii+1, nSigBp.getVal());
       varPhiMass->SetBinContent(ii+1, nSigBp.getError());

 }
     TCanvas * c=new TCanvas("c","c",800,600);
     varPhiMass->Draw();
 //  varPhiMass->GetYaxis()->SetRangeUser(0., 4000.);

}
