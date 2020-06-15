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
#include "TMath.h"
#include "TH2F.h"

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
     //Timing
     sprintf(title,"xtime_layer_1hit_%i",ij);
     xtime_layer_1hit[ij] = new TH1F(title,title,8000,-0.5,799.5);
     sprintf(title,"ytime_layer_1hit_%i",ij);
     ytime_layer_1hit[ij] = new TH1F(title,title,8000,-0.5,799.5);
     sprintf(title,"xtime_layer_hough_%i",ij);
     xtime_layer_hough[ij] = new TH1F(title,title,8000,-0.5,799.5);
     sprintf(title,"ytime_layer_hough_%i",ij);
     ytime_layer_hough[ij] = new TH1F(title,title,8000,-0.5,799.5);
     sprintf(title,"xtime_layer_%i",ij);
     xtime_layer[ij] = new TH1F(title,title,8000,-0.5,799.5);
     sprintf(title,"ytime_layer_%i",ij);
     ytime_layer[ij] = new TH1F(title,title,8000,-0.5,799.5);
     sprintf(title,"distance_from_pos_l%i",ij);
     distance_from_pos_l[ij] = new TH1F(title,title,200,-0.5,199.5);
     sprintf(title,"distance_from_ext_l%i",ij);
     distance_from_ext_l[ij] = new TH1F(title,title,200,-0.5,199.5);
     sprintf(title,"tdc_mul_twohit_x_l%i",ij);
     tdc_mul_twohit_x[ij] = new TH1F(title,title,40,-0.5,39.5);
     sprintf(title,"tdc_mul_twohit_y_l%i",ij);
     tdc_mul_twohit_y[ij] = new TH1F(title,title,40,-0.5,39.5);
   }

   sprintf(title,"distance_from_pos");
   distance_from_pos = new TH1F(title,title,200,-0.5,199.5);
   sprintf(title,"distance_from_ext");
   distance_from_ext = new TH1F(title,title,200,-0.5,199.5);

   sprintf(title,"tdc_mul_mustopped_x");
   tdc_mul_mustopped[0] = new TH1F(title,title,40,-0.5,39.5);
   sprintf(title,"tdc_mul_mustopped_y");
   tdc_mul_mustopped[1] = new TH1F(title,title,40,-0.5,39.5);

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


   //////////////////////////////////////////////
   //                                          //
   //           Position Analysis              //
   //                                          //
   //////////////////////////////////////////////

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
   vector<double> xptsall[nlayer]; //For trigger criteria
   vector<double> yptsall[nlayer];
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
         // Only layers with one hit or one cluster is used
         if (ix<xhits[jk]-1 && abs(xpts[jk][ix]-xpts[jk][ix+1])>1)
         {xpos[jk]=-100; break;}
         xpos[jk] += xpts[jk][ix];
       }
       int tempx=0;
       int mul=1;
       for(int ix=0; ix<xhits[jk]; ix++) {
         //Find clusters and save to xptsall[jk]
         if(mul==1){tempx += xpts[jk][ix];}
         if(ix<xhits[jk]-1 &&  abs(xpts[jk][ix]-xpts[jk][ix+1])==1)
         {tempx += xpts[jk][ix+1]; mul++;}
         else {
           double val = (double)tempx/(double)mul + 0.5 - xoff[jk];
           val -= cal_slope2(val, &align_xstr_xdev[jk][0]);
           xptsall[jk].push_back(val);
           mul=1;
           tempx=0;
         }
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
         //Only layers with one hit or one cluster is used.
         if (ix<yhits[jk]-1 && abs(ypts[jk][ix]-ypts[jk][ix+1])>1)
         {ypos[jk]=-100; break;}
         ypos[jk] += ypts[jk][ix];
       }
       int tempy=0;
       int mul=1;
       for(int ix=0; ix<yhits[jk]; ix++) {
         //Find clusters and save to yptsall[jk]
         if(mul==1){tempy += ypts[jk][ix];}
         if(ix<yhits[jk]-1 &&  abs(ypts[jk][ix]-ypts[jk][ix+1])==1)
         {tempy += ypts[jk][ix+1]; mul++;}
         else {
           double val = (double)tempy/(double)mul + 0.5 - yoff[jk];
           val -= cal_slope2(val, &align_ystr_ydev[jk][0]);
           yptsall[jk].push_back(val);
           mul=1;
           tempy=0;
         }
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

   if(occulyr<nlayer)
   {
     if(xfitfailed==0)
     {
       double minx=10000.0;
       for(unsigned int ix=0; ix<xptsall[occulyr].size(); ix++)
       {
         if(abs(xext[occulyr]-xptsall[occulyr][ix])<minx)
         {
           tempxpos[occulyr] = xptsall[occulyr][ix];
           minx = abs(xext[occulyr]-xptsall[occulyr][ix]);
         }
       }
     }
     if(yfitfailed==0)
     {
       double miny=10000.0;
       for(unsigned int ix=0; ix<yptsall[occulyr].size(); ix++)
       {
         if(abs(yext[occulyr]-yptsall[occulyr][ix])<miny)
         {
           tempypos[occulyr] = yptsall[occulyr][ix];
           miny = abs(yext[occulyr]-yptsall[occulyr][ix]);
         }
       }
     }
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


   //////////////////////////////////////////////
   //                                          //
   //           Timing Analysis                //
   //                                          //
   //////////////////////////////////////////////

   //X-time
   int    xtfitfailed,xtndof;
   double xtime[nlayer];
   double xtintersect,xtslope,xterrinter,xterrslope,xterrcovar,xtchisquare;
   double xtext[nlayer],xtdev[nlayer],xtexter[nlayer];
   //Y-time
   int    ytfitfailed,ytndof;
   double ytime[nlayer];
   double ytintersect,ytslope,yterrinter,yterrslope,yterrcovar,ytchisquare;
   double ytext[nlayer],ytdev[nlayer],ytexter[nlayer];

   //removing noise hits using hough transform______________________________
   #ifdef hough_space
   // double rmin,rmax,thetamin,thetamax;
   // int nbinr,nbintheta;
   // rmin=-10000.;
   // rmax=10000.;
   // thetamin=0.;
   // thetamax=180.;
   // nbinr=1000;
   // nbintheta=1000;
   //
   // TH2F hough_space_x("hough_space_x","rx vs thetax",nbinr,rmin,rmax,nbintheta,thetamin,thetamax);
   // TH2F hough_space_y("hough_space_y","ry vs thetay",nbinr,rmin,rmax,nbintheta,thetamin,thetamax);
   #endif

   for(int jk=0;jk<nlayer;jk++)
   {
     for (int itdc=0; itdc<nTDCpLayer; itdc++)
     {
       double tmpxtime = -10000;
       if(vxtdc_l[jk][itdc]->size())
       {
         tmpxtime=600.+ 0.1*int(vxtdc_l[jk][itdc]->at(0)-tdc_ref_l[jk]);
         xtime_layer_1hit[jk]->Fill(tmpxtime);
       }
       double tmpytime = -10000;
       if(vytdc_l[jk][itdc]->size())
       {
         tmpytime=600.+ 0.1*int(vytdc_l[jk][itdc]->at(0)-tdc_ref_l[jk]);
         ytime_layer_1hit[jk]->Fill(tmpytime);

       }
       #ifdef hough_space
       // for(unsigned int multitdc=0;multitdc<vxtdc_l[jk][itdc]->size();multitdc++)
       // {
       //   double z = (double)jk;
       //   double xt = 600.+ 0.1*(double)(vxtdc_l[jk][itdc]->at(multitdc)-tdc_ref_l[jk]);
       //   for(double theta=thetamin;theta<=thetamax;theta=theta+180./1000.)
       //   {
       //     double r=xt*cos(TMath::DegToRad()*theta)+z*sin(TMath::DegToRad()*theta);
       //     hough_space_x.Fill(r,theta);
       //   }
       // }
       //
       // for(unsigned int multitdc=0;multitdc<vytdc_l[jk][itdc]->size();multitdc++)
       // {
       //   double z = (double)jk;
       //   double yt = 600.+ 0.1*(double)(vytdc_l[jk][itdc]->at(multitdc)-tdc_ref_l[jk]);
       //   for(double theta=thetamin;theta<=thetamax;theta=theta+180./1000.)
       //   {
       //     double r=yt*cos(TMath::DegToRad()*theta)+z*sin(TMath::DegToRad()*theta);
       //     hough_space_y.Fill(r,theta);
       //   }
       // }
       #endif
     }//nTDCpLayer
   }//nlayer

   #ifdef hough_space
   // Int_t rlocx,thetalocx,zlocx;
   // hough_space_x.GetBinXYZ(hough_space_x.GetMaximumBin(), rlocx, thetalocx, zlocx);
   // double rvalx=rmin+(rmax-rmin)*((double)rlocx-0.5)/nbinr;
   // double thetavalx=thetamin+(thetamax-thetamin)*((double)thetalocx-0.5)/nbintheta;
   // Int_t rlocy,thetalocy,zlocy;
   // hough_space_y.GetBinXYZ(hough_space_y.GetMaximumBin(), rlocy, thetalocy, zlocy);
   // double rvaly=rmin+(rmax-rmin)*((double)rlocy-0.5)/nbinr;
   // double thetavaly=thetamin+(thetamax-thetamin)*((double)thetalocy-0.5)/nbintheta;
   #endif


   double tmpxtime[nlayer];
   double tmpytime[nlayer];
   for(int jk=0;jk<nlayer;jk++)
   {
     tmpxtime[jk]=-10000.0;
     tmpytime[jk]=-10000.0;
     double z=jk;
     double xt=xtime_layer_1hit[jk]->GetMean();
     double yt=xtime_layer_1hit[jk]->GetMean();
     #ifdef hough_space
     //double xt=(rvalx-(z*sin(TMath::DegToRad()*thetavalx)))/cos(TMath::DegToRad()*thetavalx);
     //double yt=(rvaly-(z*sin(TMath::DegToRad()*thetavaly)))/cos(TMath::DegToRad()*thetavaly);
     #endif
     int minx=100000;
     int miny=100000;
     for (int itdc=0; itdc<nTDCpLayer; itdc++)
     {
       for(unsigned int multitdc=0;multitdc<vxtdc_l[jk][itdc]->size();multitdc++)
       {
         if(abs(xt - (600. + 0.1*int(vxtdc_l[jk][itdc]->at(multitdc) - tdc_ref_l[jk])))<minx)
         {
           minx=abs(xt - (600. + 0.1*int(vxtdc_l[jk][itdc]->at(multitdc) - tdc_ref_l[jk])));
           tmpxtime[jk]=600. + 0.1*int(vxtdc_l[jk][itdc]->at(multitdc) - tdc_ref_l[jk]);
         }
       }
       for(unsigned int multitdc=0;multitdc<vytdc_l[jk][itdc]->size();multitdc++)
       {
         if(abs(yt - (600. + 0.1*int(vytdc_l[jk][itdc]->at(multitdc) -tdc_ref_l[jk])))<miny)
         {
           miny=abs(yt - (600. + 0.1*int(vytdc_l[jk][itdc]->at(multitdc) -tdc_ref_l[jk])));
           tmpytime[jk]=600. + 0.1*int(vytdc_l[jk][itdc]->at(multitdc) -tdc_ref_l[jk]);
         }
       }

     }
     xtime_layer_hough[jk]->Fill(tmpxtime[jk]);
     ytime_layer_hough[jk]->Fill(tmpytime[jk]);
   }

   //Timing Corrections start

   double initxpos, initypos, initzpos;
   double dist[nlayer];
   initxpos = -100.0;
   initypos = -100.0;
   initzpos = 0;
   int init = -1;
   int istr;
   int istrxtime[nlayer], istrytime[nlayer]; //strip with early timing

   for(int jk=0; jk<nlayer; jk++)
   {
     istrxtime[jk] = -1;
     istrytime[jk] = -1;
     xtime[jk] = -10000.0;
     ytime[jk] = -10000.0;
     if (init<0) {
 		  initxpos = xextloc[jk];
 		  initypos = yextloc[jk];
 		  dist[jk] = 0.0;
 		  init = jk;
 		  initzpos = layerzpos[jk];
 		  } else {
        //Distance form bottom layer(along the partical trajectory)
 		     dist[jk] = sqrt( pow((xextloc[jk] - initxpos)*stripwidth, 2.) +
 				   pow((yextloc[jk] - initypos)*stripwidth, 2.) +
 				   pow(layerzpos[jk] - initzpos, 2.));
 		  }
      if(abs(xdev[jk]) < 2.0 && abs(ydev[jk]) < 2.0 &&
       xpts[jk].size()>0 && xpts[jk].size() <=nmxhits && ypts[jk].size()>0 && ypts[jk].size() <=nmxhits)
       {//if condition for timing correction starts here
        int channelnum[nstrip];
        double inittime = 11000;
        double tshft=0;
        //X-Side
        //Find the strip with the first time
        for (int ix=0; ix<min(nstrip,int(xpts[jk].size())); ix++)
        {
          channelnum[ix] = xpts[jk][ix]%8;
          if (vxtdc_l[jk][channelnum[ix]]->size()) {
            double tmpoxtime = 600.0 + 0.1*(int(vxtdc_l[jk][channelnum[ix]]->at(0)-tdc_ref_l[jk]));
            if (tmpoxtime < inittime) {
    			    inittime = tmpoxtime; tshft =xtoffset[jk][xpts[jk][ix]]; istr = istrxtime[jk] = xpts[jk][ix];
    			  }
          }
        }
        //Save the time to xtime
        if(inittime <10000){xtime[jk] = inittime;}
        //Apply corrections
        xtime[jk] -= timeoffsetx[jk];
        xtime[jk] -= slope_path*yext[jk];
        xtime[jk] -= tshft;

        istr = int(yextloc[jk]+0.5);
        if (istr<0) istr=0;
        if (istr>=nstrip) istr = nstrip-1;
        double dx = yextloc[jk]-istr;

        // Linear extrapolation using only two points
  		  if ((istr==0 && dx<=0.0) || (istr==nstrip-1 && dx>=0.0)) {
   		    xtime[jk] -=xt_slope_cor[jk][istrxtime[jk]][istr];
   		  } else if (dx>0) {
   		    xtime[jk] -=(1-dx)*xt_slope_cor[jk][istrxtime[jk]][istr]+dx*xt_slope_cor[jk][istrxtime[jk]][istr+1];
   		  } else {
   		    xtime[jk] -=abs(dx)*xt_slope_cor[jk][istrxtime[jk]][istr-1]+(1-abs(dx))*xt_slope_cor[jk][istrxtime[jk]][istr];
   		  }

        inittime = 11000;
        tshft=0;
        //Y-Side
        //Find the strip with the first time
        for (int ix=0; ix<min(nstrip,int(ypts[jk].size())); ix++)
        {
          channelnum[ix] = ypts[jk][ix]%8;
          if (vytdc_l[jk][channelnum[ix]]->size()) {
            double tmpoytime = 600.0 + 0.1*(int(vytdc_l[jk][channelnum[ix]]->at(0)-tdc_ref_l[jk]));
            if (tmpoytime < inittime) {
    			    inittime = tmpoytime; tshft =ytoffset[jk][ypts[jk][ix]]; istr = istrytime[jk] = ypts[jk][ix];
    			  }
          }
        }
        //Save the time to xtime
        if(inittime <10000){ytime[jk] = inittime;}
        //Apply corrections
        ytime[jk] -= timeoffsety[jk];
        ytime[jk] -= ytimeshift; // shift y-time to match with x-time
        ytime[jk] -= slope_path*xext[jk];
        ytime[jk] -= tshft;

        istr = int(xextloc[jk]+0.5);
        if (istr<0) istr=0;
        if (istr>=nstrip) istr = nstrip-1;
        dx = xextloc[jk]-istr;

        // Linear extrapolation using only two points
  		  if ((istr==0 && dx<=0.0) || (istr==nstrip-1 && dx>=0.0)) {
   		    ytime[jk] -=yt_slope_cor[jk][istrytime[jk]][istr];
   		  } else if (dx>0) {
   		    ytime[jk] -=(1-dx)*yt_slope_cor[jk][istrytime[jk]][istr]+dx*yt_slope_cor[jk][istrytime[jk]][istr+1];
   		  } else {
   		    ytime[jk] -=abs(dx)*yt_slope_cor[jk][istrytime[jk]][istr-1]+(1-abs(dx))*yt_slope_cor[jk][istrytime[jk]][istr];
   		  }
   //Timing Corrections end

        xtime_layer[jk]->Fill(xtime[jk]);
        ytime_layer[jk]->Fill(ytime[jk]);

      }  //if condition for timing correction ends here
    }//nlayer for loop


   //////////////////////////////////////////////
   //                                          //
   //          Dead Space Analysis             //
   //                                          //
   //////////////////////////////////////////////

   //Trigger layers - 6,7,8,9
   //We are intrested layers 3,4,5(where muon stopped).
   //means 4,5,6 are the layers we look for efficency after a hit
   //Bottom layers are excluded for our study because of chance coincidence that
   //bottom layers may be inefficient

   bool mustopped;
   for(int laystudy=4;laystudy<=6;laystudy++)
   {
     //The layer under study should have hit, if not check above layers
     //Reasons for layer under study dont have a hit are inefficiency or
     //the muon is stopped in the above layer
     if(!xusedpos[laystudy])
     {continue;}
     //Three layers just below the layer under study should have the expected location
     if(xext[laystudy-1]>59 || xext[laystudy-1]<4 || yext[laystudy-1]>59 ||yext[laystudy-1]<4
      ||xext[laystudy-2]>59 || xext[laystudy-2]<4 || yext[laystudy-2]>59 ||yext[laystudy-2]<4
      ||xext[laystudy-3]>59 || xext[laystudy-3]<4 || yext[laystudy-3]>59 ||yext[laystudy-3]<4)
     {
       break;
     }
     mustopped=true;
     //bottom layers should not have a latch hit within two strip width of the expected location
     for(int botlay=laystudy-1;botlay>=0;botlay--)
     {
       if(abs(xpos[botlay]-xext[botlay]) <= 2)
       {mustopped=false;break;}
     }
     if(mustopped==false)
     {break;}
     //Now we are sure that the muon stopped. But electron can go isotropically
     //So we dont know the laystudy layer hit is due to a muon or electron
     //To confirm it is due to a muon, we impose the time should be within 10 ns
     //and position to be within one stripwidth.

     //if(xtime[laystudy])                  ---incorporate later
     double muhitx=-100;
     double muhity=-100;
     double elhitx=-100;
     double elhity=-100;

     if(xptsall[laystudy].size()==2  && yptsall[laystudy].size()==2)
     {
        if(abs(xext[laystudy]-xptsall[laystudy][0])<abs(xext[laystudy]-xptsall[laystudy][1])) {
          muhitx=xptsall[laystudy][0];
          elhitx=xptsall[laystudy][1]; }
        else {
          muhitx=xptsall[laystudy][1];
          elhitx=xptsall[laystudy][0]; }

        if(abs(yext[laystudy]-yptsall[laystudy][0])<abs(yext[laystudy]-yptsall[laystudy][1])) {
            muhity=yptsall[laystudy][0];
            elhity=yptsall[laystudy][1]; }
          else {
            muhity=yptsall[laystudy][1];
            elhity=yptsall[laystudy][0]; }

      double distancepos = sqrt(pow((muhitx - elhitx)*stripwidth, 2.) +
          pow((muhity - elhity)*stripwidth, 2.));
        double distanceext = sqrt(pow((xext[laystudy] - elhitx)*stripwidth, 2.) +
          pow((yext[laystudy] - elhity)*stripwidth, 2.));

        distance_from_pos->Fill(distancepos);
        distance_from_ext->Fill(distanceext);

     }

     int mul=0;
     for (int itdc=0; itdc<nTDCpLayer; itdc++)
     {
       mul += vxtdc_l[laystudy][itdc]->size();
     }
     tdc_mul_mustopped[0]->Fill(mul);

     mul=0;
     for (int itdc=0; itdc<nTDCpLayer; itdc++)
     {
       mul += vytdc_l[laystudy][itdc]->size();
     }
     tdc_mul_mustopped[1]->Fill(mul);
   }

   for(int jk=0; jk<nlayer; jk++)
   {
     if(xfitfailed==1 || yfitfailed==1)
     {
       break;
     }
     //if(xtime[laystudy])                  ---incorporate later
     double muhitx=-100;
     double muhity=-100;
     double elhitx=-100;
     double elhity=-100;

     if(xptsall[jk].size()==2  && yptsall[jk].size()==2)
     {
       if(abs(xext[jk]-xptsall[jk][0])<abs(xext[jk]-xptsall[jk][1])) {
         muhitx=xptsall[jk][0];
         elhitx=xptsall[jk][1]; }
       else {
         muhitx=xptsall[jk][1];
         elhitx=xptsall[jk][0]; }
       if(abs(yext[jk]-yptsall[jk][0])<abs(yext[jk]-yptsall[jk][1])) {
           muhity=yptsall[jk][0];
           elhity=yptsall[jk][1]; }
         else {
           muhity=yptsall[jk][1];
           elhity=yptsall[jk][0]; }

        double distancepos = sqrt(pow((muhitx - elhitx)*stripwidth, 2.) +
          pow((muhity - elhity)*stripwidth, 2.));
        double distanceext = sqrt(pow((xext[jk] - elhitx)*stripwidth, 2.) +
          pow((yext[jk] - elhity)*stripwidth, 2.));

        distance_from_pos_l[jk]->Fill(distancepos);
        distance_from_ext_l[jk]->Fill(distanceext);

        int mul=0;
        for (int itdc=0; itdc<nTDCpLayer; itdc++)
        {
          mul += vxtdc_l[jk][itdc]->size();
        }
        tdc_mul_twohit_x[jk]->Fill(mul);

        mul=0;
        for (int itdc=0; itdc<nTDCpLayer; itdc++)
        {
          mul += vytdc_l[jk][itdc]->size();
        }
        tdc_mul_twohit_y[jk]->Fill(mul);

     }
   }




   //////////////////////////////////////////////
   //                                          //
   //          Dead Time Analysis              //
   //                                          //
   //////////////////////////////////////////////













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
     for(unsigned int ij=0;ij<nlayer;ij++){
       xtime_layer_1hit[ij]->Write(0, TObject::kOverwrite);
       ytime_layer_1hit[ij]->Write(0, TObject::kOverwrite);
     }
     for(unsigned int ij=0;ij<nlayer;ij++){
       xtime_layer_hough[ij]->Write(0, TObject::kOverwrite);
       ytime_layer_hough[ij]->Write(0, TObject::kOverwrite);
     }
     for(unsigned int ij=0;ij<nlayer;ij++){
       xtime_layer[ij]->Write(0, TObject::kOverwrite);
       ytime_layer[ij]->Write(0, TObject::kOverwrite);
     }
     distance_from_pos->Write(0, TObject::kOverwrite);
     distance_from_ext->Write(0, TObject::kOverwrite);
     for(unsigned int ij=0;ij<nlayer;ij++){
       distance_from_pos_l[ij]->Write(0, TObject::kOverwrite);
       distance_from_ext_l[ij]->Write(0, TObject::kOverwrite);
     }
     tdc_mul_mustopped[0]->Write(0, TObject::kOverwrite);
     tdc_mul_mustopped[1]->Write(0, TObject::kOverwrite);
     for(unsigned int ij=0;ij<nlayer;ij++){
       tdc_mul_twohit_x[ij]->Write(0, TObject::kOverwrite);
       tdc_mul_twohit_y[ij]->Write(0, TObject::kOverwrite);
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

//After second fitting if we called GetPosInStrip two times alignment correction will be done. Is this correct?
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
