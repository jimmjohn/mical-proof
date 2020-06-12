#define Selector_cxx
// The class definition in Selector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("Selector.C")
// root> T->Process("Selector.C","some options")
// root> T->Process("Selector.C+")
//


#include "Selector.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>
#include "TPostScript.h"
#include <TParameter.h>
#include "StraightLineFit.h"
#include <algorithm>


void Selector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   Reset();

   //print the option specified in the Process function.
   TString option = GetOption();
   Info("Begin", "starting : %s", option.Data());



}

void Selector::SlaveBegin(TTree  *tree)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   //initialize the Tree branch addresses
   Init(tree);

   // We may be creating a dataset or a merge file: check it
   TNamed *nm = dynamic_cast<TNamed *>(fInput->FindObject("SimpleNtuple.root"));
   if (nm) {
      // Just create the object
      UInt_t opt = TProofOutputFile::kRegister | TProofOutputFile::kOverwrite | TProofOutputFile::kVerify;
      fProofFile = new TProofOutputFile("SimpleNtuple.root",
                                        TProofOutputFile::kDataset, opt, nm->GetTitle());
   } else {
      // For the ntuple, we use the automatic file merging facility
      // Check if an output URL has been given
      TNamed *out = (TNamed *) fInput->FindObject("PROOF_OUTPUTFILE_LOCATION");
      Info("SlaveBegin", "PROOF_OUTPUTFILE_LOCATION: %s", (out ? out->GetTitle() : "undef"));
      fProofFile = new TProofOutputFile("SimpleNtuple.root", (out ? out->GetTitle() : "M"));
      out = (TNamed *) fInput->FindObject("PROOF_OUTPUTFILE");
      if (out) fProofFile->SetOutputFileName(out->GetTitle());
   }

   // Open the file
   fFile = fProofFile->OpenFile("RECREATE");
   if (fFile && fFile->IsZombie()) SafeDelete(fFile);

   // Cannot continue
   if (!fFile) {
      Info("SlaveBegin", "could not create '%s': instance is invalid!", fProofFile->GetName());
      return;
   }


   //Histograms
   char title[100];
   for (int ij=0; ij<nlayer; ij++) {
     sprintf(title, "xlayer_occu_l%i", ij);
     xlayer_occu[ij] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);
     xlayer_occu[ij]->SetDirectory(fFile);
     sprintf(title, "ylayer_occu_l%i", ij);
     ylayer_occu[ij] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);
     sprintf(title, "xlayer_mult_l%i", ij);
     xlayer_mult[ij] = new TH1F(title, title, nstrip+1, -0.5, nstrip+0.5);
     sprintf(title, "ylayer_mult_l%i", ij);
     ylayer_mult[ij] = new TH1F(title, title, nstrip+1, -0.5, nstrip+0.5);
     sprintf(title, "rawhits_corr_xymul_l%i", ij);
     rawhits_corr_xymul[ij] = new TH2F(title, title, nstrip+1, -0.5, nstrip+0.5,  nstrip+1, -0.5, nstrip+0.5);
     sprintf(title, "raw_occu_l%i", ij);
     raw_occu[ij] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
     for (int jk=ij+1; jk<nlayer; jk++) {
       sprintf(title, "rawhits_xlay_corr_mul_l%i_l%i", ij, jk);
       rawhits_xlay_corr_mul[ij][jk] = new TH2F(title, title, nstrip+1, -0.5, nstrip+0.5,  nstrip+1, -0.5, nstrip+0.5);
       sprintf(title, "rawhits_ylay_corr_mul_l%i_l%i", ij, jk);
       rawhits_ylay_corr_mul[ij][jk] = new TH2F(title, title, nstrip+1, -0.5, nstrip+0.5,  nstrip+1, -0.5, nstrip+0.5);
     }
     double xposmx = 6.0;
     sprintf(title, "xlayer_reso_l%i", ij);
     xlayer_reso[ij] = new TH1F(title, title, 150, -xposmx, xposmx);
     sprintf(title, "ylayer_reso_l%i", ij);
     ylayer_reso[ij] = new TH1F(title, title, 150, -xposmx, xposmx);
     sprintf(title, "totalentry_l%i", ij);
     totalentry[ij]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

   }

   //create histograms
   h1 = new TH1F("h1","h1",1000,10000,240000);
   fOutput->Add(h1);
   //print the option specified in the Process function.
   TString option = GetOption();
}

Bool_t Selector::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fChain->GetTree()->GetEntry(entry);
   h1->Fill(EveTS[0]->GetTime());
   //Reading from Main root file.


   //Position Analysis
   //-----------------

   //X-side
   int    xhits[nlayer],xfitfailed,xndof;
   double xpos[nlayer],xxerr[nlayer],xintersect,xslope,xerrinter,xerrslope,xerrcovar,xchisquare;
   double xext[nlayer],xdev[nlayer],xexter[nlayer],xextloc[nlayer],xposinstr[nlayer];
   bool   xusedpos[nlayer];

   //Y-side
   int    yhits[nlayer],yfitfailed,yndof;
   double ypos[nlayer],yyerr[nlayer],yintersect,yslope,yerrinter,yerrslope,yerrcovar,ychisquare;
   double yext[nlayer],ydev[nlayer],yexter[nlayer],yextloc[nlayer],yposinstr[nlayer];
   bool   yusedpos[nlayer];

   //Additional variables
   vector<int> xpts[nlayer];
   vector<int> ypts[nlayer];
   vector<int> xptsall[nlayer]; //For trigger criteria
   vector<int> yptsall[nlayer];
   vector<int> xptsalltdc[nlayer][nTDCpLayer]; //signal for each TDC
   vector<int> yptsalltdc[nlayer][nTDCpLayer];


   for(int jk=0;jk<nlayer;jk++) {
     xpts[jk].clear(); ypts[jk].clear(); xptsall[jk].clear(); yptsall[jk].clear();
     for (int kl=0; kl<nTDCpLayer; kl++) {
       xptsalltdc[jk][kl].clear();  yptsalltdc[jk][kl].clear();
     }
     for(int kl=0; kl<nstrip; kl++) {
       if(xLayer[jk]->TestBitNumber(kl)) {
         xpts[jk].push_back(kl);
         xptsalltdc[jk][kl%8].push_back(kl);
       }
       if(yLayer[jk]->TestBitNumber(kl)) {
         ypts[jk].push_back(kl);
         yptsalltdc[jk][kl%8].push_back(kl);
       }
     }
     for (unsigned int ix=0; ix<xpts[jk].size(); ix++) {
       xlayer_occu[jk]->Fill(xpts[jk][ix]);
     }
     for (unsigned int iy=0; iy<ypts[jk].size(); iy++) {
       ylayer_occu[jk]->Fill(ypts[jk][iy]);
     }
     xlayer_mult[jk]->Fill(xpts[jk].size());
     ylayer_mult[jk]->Fill(ypts[jk].size());
     rawhits_corr_xymul[jk]->Fill(xpts[jk].size(), ypts[jk].size());
     for (unsigned int ix=0; ix<xpts[jk].size(); ix++) {
       for (unsigned int iy=0; iy<ypts[jk].size(); iy++) {
         raw_occu[jk]->Fill(xpts[jk][ix], ypts[jk][iy]);
       }
     }
     for (int ij=jk+1; ij<nlayer; ij++) {
       rawhits_xlay_corr_mul[jk][ij]->Fill(xpts[jk].size(), xpts[ij].size());
       rawhits_ylay_corr_mul[jk][ij]->Fill(ypts[jk].size(), ypts[ij].size());
     }
     xhits[jk] = xpts[jk].size();
     yhits[jk] = ypts[jk].size();
     xdev[jk] = 100; xpos[jk]=0.0;
     ydev[jk] = 100; ypos[jk]=0.0;
     //X-Side
     if (xhits[jk]<=0 || xhits[jk]>nmxhits) {
       xpos[jk]= -100;
     } else {
       for (int ix=0; ix<xhits[jk]; ix++) {
         // Only first cluster will be used.... Need to modify to big cluster
         if (ix<xhits[jk]-1 && abs(xpts[jk][ix]-xpts[jk][ix+1])>1) { xpos[jk]=-100; break;}
         xpos[jk] +=xpts[jk][ix];
       }
     }
     xxerr[jk] = errxco[jk]*errxco[jk];
     if (xpos[jk]>=0.0) {
       xpos[jk]  = xpos[jk]/xhits[jk] + 0.5 - xoff[jk];
       //Aligmnent Correction
       xpos[jk] -= cal_slope2(xpos[jk], &align_xstr_xdev[jk][0]);
       xxerr[jk] = xposerrsq[xhits[jk]-1][jk];
     }
     //Y-Side
     if (yhits[jk]<=0 || yhits[jk]>nmxhits) {
       ypos[jk]= -100;
     } else {
       for (int ix=0; ix<yhits[jk]; ix++) {
         //Only first cluster will be used.... Need to modify to big cluster
         if (ix<yhits[jk]-1 && abs(ypts[jk][ix]-ypts[jk][ix+1])>1) { ypos[jk]=-100; break;}
         ypos[jk] +=ypts[jk][ix];
       }
     }
     yyerr[jk] = errxco[jk]*errxco[jk];
     if (ypos[jk]>=0.0) {
       ypos[jk]  = ypos[jk]/yhits[jk] + 0.5 - yoff[jk];
       //Aligmnent Correction
       ypos[jk] -= cal_slope2(ypos[jk], &align_ystr_ydev[jk][0]);
       yyerr[jk] = yposerrsq[yhits[jk]-1][jk];
     }
     //Sort out hits, which can be used for fit
     xusedpos[jk] = (xpos[jk]>-99 && xhits[jk]<=nmxusedhits) ? true : false;
     yusedpos[jk] = (ypos[jk]>-99 && yhits[jk]<=nmxusedhits) ? true : false;
   }//Filled everything for fitting

   //Position StraightLineFit

   //Temporary storing of position is required. Otherwise after each fitting alignment correction is done
   double tempxpos[nlayer], tempypos[nlayer];
   copy(begin(xpos), end(xpos), begin(tempxpos));
   copy(begin(ypos), end(ypos), begin(tempypos));

   for(int occulyr=0;occulyr<=nlayer;occulyr++)
   {
    double zval[nlayer];
    for (unsigned int ix=0; ix<nlayer; ix++) {
      zval[ix]=layerzpos[ix];
      xext[ix]= xextloc[ix] = xexter[ix] =xposinstr[ix] =  100;
      yext[ix]= yextloc[ix] = yexter[ix] =yposinstr[ix] =  100;
    }
    copy(begin(tempxpos), end(tempxpos), begin(xpos));
    copy(begin(tempypos), end(tempypos), begin(ypos));

    xndof=0;
    xfitfailed=0;
    xchisquare=0;
    //X-Side fit
    StraightLineFit xposfit(1, zval, xpos,  xxerr, xusedpos, occulyr, occulyr, layfirst, laylast, xyPosDev);
    xposfit.GetParameters(xfitfailed, xintersect, xslope);
    xposfit.GetFitValues(xext, xdev, xexter);

    yndof=0;
    yfitfailed=0;
    ychisquare=0;
   //Y-Side fit
   StraightLineFit yposfit(1, zval, ypos,  yyerr, yusedpos, occulyr, occulyr, layfirst, laylast, xyPosDev);
   yposfit.GetParameters(yfitfailed, yintersect, yslope);
   yposfit.GetFitValues(yext, ydev, yexter);
   for (int ix=0; ix<nlayer; ix++) {
       xpos[ix] -=cal_slope2(yext[ix], &align_ystr_xdev[ix][0]);
       ypos[ix] -=cal_slope2(xext[ix], &align_xstr_ydev[ix][0]);
   }

   //X-Side fit
   xposfit = StraightLineFit(1, zval, xpos,  xxerr, xusedpos, occulyr, occulyr, layfirst, laylast, xyPosDev);
   xposfit.GetParameters(xfitfailed, xintersect, xslope);
   xposfit.GetFitValues(xext, xdev, xexter);
   xposfit.GetChisqure(xndof,xchisquare);
   GetPosInStrip(0, xext, yext, xoff, xposinstr, xextloc);
   //Y-Side fit
   yposfit = StraightLineFit(1, zval, ypos,  yyerr, yusedpos, occulyr, occulyr, layfirst, laylast, xyPosDev);
   yposfit.GetParameters(yfitfailed, yintersect, yslope);
   yposfit.GetFitValues(yext, ydev, yexter);
   yposfit.GetChisqure(yndof,ychisquare);
   GetPosInStrip(1, yext, xext, yoff, yposinstr, yextloc);
   //if(occulyr==0){cout << "layer 0 \t xext = " << xext[occulyr] << "\txpos = " << xpos[occulyr] << endl;}
   //if(occulyr==nlayer){cout << "combined \t xext = " << xext[0] << "\txpos = " << xpos[0] << endl;}

   if(occulyr<nlayer)
   {
   // if (abs(xdev[occulyr])<6.0 && xusedpos[occulyr] && xndof>=nmnhits && xchisquare/(xndof-2)<mxchisq && xfitfailed==0) {
   //   xlayer_reso[occulyr]->Fill(xdev[occulyr]);
   // }
   // if (abs(ydev[occulyr])<6.0 && yusedpos[occulyr] && yndof>=nmnhits && ychisquare/(yndof-2)<mxchisq && yfitfailed==0) {
   //   ylayer_reso[occulyr]->Fill(ydev[occulyr]);
   // }
   // totalentry[occulyr]->Fill(xextloc[occulyr],yextloc[occulyr]);

   }// occulyr < nlayer
   else{
     for(int jk=0;jk<nlayer;jk++)
     {
       if (abs(xdev[jk])<6.0 && xusedpos[jk] && xndof>=nmnhits && xchisquare/(xndof-2)<mxchisq && xfitfailed==0) {
         xlayer_reso[jk]->Fill(xdev[jk]);
       }
       if (abs(ydev[jk])<6.0 && yusedpos[jk] && yndof>=nmnhits && ychisquare/(yndof-2)<mxchisq && yfitfailed==0) {
         ylayer_reso[jk]->Fill(ydev[jk]);
       }
       totalentry[jk]->Fill(xextloc[jk],yextloc[jk]);
     }
   } //occulyr = nlayer

   }//occulyr loop end


   //Timing analysis
   //---------------

   //X-time
   int    xtfitfailed,xtndof;
   double xtintersect,xtslope,xterrinter,xterrslope,xterrcovar,xtchisquare;
   double xtext[nlayer],xtdev[nlayer],xtexter[nlayer];
   //Y-time
   int    ytfitfailed,ytndof;
   double ytintersect,ytslope,yterrinter,yterrslope,yterrcovar,ytchisquare;
   double ytext[nlayer],ytdev[nlayer],ytexter[nlayer];




   return kTRUE;
}

void Selector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
   if (fFile) {
     Bool_t cleanup = kFALSE;
     TDirectory *savedir = gDirectory;
     fFile->cd();
     //Writing all histograms
     h1->Write(0, TObject::kOverwrite);
     for(unsigned int ij=0;ij<nlayer;ij++){
       xlayer_occu[ij]->Write(0, TObject::kOverwrite);
       ylayer_occu[ij]->Write(0, TObject::kOverwrite);
     }
     for(unsigned int ij=0;ij<nlayer;ij++){
       xlayer_mult[ij]->Write(0, TObject::kOverwrite);
       ylayer_mult[ij]->Write(0, TObject::kOverwrite);
     }
     for(unsigned int ij=0;ij<nlayer;ij++){
       rawhits_corr_xymul[ij]->Write(0, TObject::kOverwrite);
     }
     for(unsigned int ij=0;ij<nlayer;ij++){
       raw_occu[ij]->Write(0, TObject::kOverwrite);
     }
     for(unsigned int ij=0;ij<nlayer;ij++){
       for(unsigned int jk=ij+1;jk<nlayer;jk++){
         rawhits_xlay_corr_mul[ij][jk]->Write(0, TObject::kOverwrite);
         rawhits_ylay_corr_mul[ij][jk]->Write(0, TObject::kOverwrite);
       }
     }
     for(unsigned int ij=0;ij<nlayer;ij++){
       xlayer_reso[ij]->Write(0, TObject::kOverwrite);
       ylayer_reso[ij]->Write(0, TObject::kOverwrite);
     }
     for(unsigned int ij=0;ij<nlayer;ij++){
       totalentry[ij]->Write(0, TObject::kOverwrite);
     }

     fProofFile->Print();
     fOutput->Add(fProofFile);

     h1->SetDirectory(0);
     gDirectory = savedir;
     fFile->Close();
     // Cleanup, if needed
     if (cleanup) {
       TUrl uf(*(fFile->GetEndpointUrl()));
       SafeDelete(fFile);
       gSystem->Unlink(uf.GetFile());
       SafeDelete(fProofFile);
     }
   }
}

void Selector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   TPostScript ps("output.ps",111);
   ps.Range(20,30); //ps.Range(10,20);
   h1 = dynamic_cast<TH1F*>(fOutput->FindObject("h1"));
   TCanvas *c1 = new TCanvas("c1","h1 canvas",200,10,700,700);
   c1->cd();
   //ps.NewPage();
   gStyle->SetOptStat(0);
   gStyle->SetPadBottomMargin(0.11);
   gStyle->SetPadTopMargin(0.07);
   gStyle->SetPadLeftMargin(0.10);
   gStyle->SetPadRightMargin(0.11);
   gStyle->SetOptTitle(1);
   h1->SetMarkerColor(2);
   h1->GetXaxis()->SetTitle("Time");
   h1->GetYaxis()->SetTitle("Muon Rate in sec");
   h1->GetXaxis()->CenterTitle();
   h1->GetYaxis()->CenterTitle();
   //h1->GetXaxis()->SetTitleSize(0.055);
   //h1->GetYaxis()->SetLabelSize(0.065);
   h1->GetXaxis()->SetTitleOffset(1.4);
   h1->GetXaxis()->SetTitleFont(62);
   h1->GetYaxis()->SetTitleFont(62);
   h1->GetXaxis()->SetLabelFont(62);
   h1->GetYaxis()->SetLabelFont(62);
   h1->GetXaxis()->SetTitleOffset(1.25);
   h1->GetXaxis()->SetTickLength(0.05);
   h1->GetYaxis()->SetTickLength(0.03);
   h1->GetXaxis()->SetLabelOffset(0.002);
   h1->GetYaxis()->SetLabelOffset(0.01);
   h1->Draw();
   c1->Update();
   ps.Close();
}


Double_t Selector::cal_slope2(Double_t x, Double_t* par) {
  if (x<nstrip/2. -0.5) {
    return par[0] + par[1]*(x - nstrip/4. +0.5);
  } else {
    double par3 = (par[2]-par[0]-par[1]*nstrip/4.)/(nstrip/4.);
    return par[2] + par3*(x - 3*nstrip/4. +0.5);
  }
}

void Selector::GetPosInStrip(int ixy, double* ext, double* otherext, double* off, double* pos, double* local) {
  for (int ij=0; ij<nlayer; ij++) {
    local[ij] = ext[ij]; //+off[ij];
    if (ixy==0) {
      local[ij] +=cal_slope2(otherext[ij], &align_ystr_xdev[ij][0]);
      local[ij] +=cal_slope2(local[ij], &align_xstr_xdev[ij][0]);
    } else {
      local[ij] +=cal_slope2(otherext[ij], &align_xstr_ydev[ij][0]);
      local[ij] +=cal_slope2(local[ij], &align_ystr_ydev[ij][0]);
    }

    local[ij] +=off[ij];
    int istr = int(local[ij]);
    if (local[ij]<0.0) {
      pos[ij] = local[ij] - istr + 0.5;
    } else {
      pos[ij] = local[ij] - istr - 0.5;
    }

    local[ij] -= 0.5;
  }
}


Selector::~Selector()
{
   // Destructor
   SafeDelete(fFile);
}
