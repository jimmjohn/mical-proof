//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr 23 18:17:11 2020 by ROOT version 6.18/04
// from TTree evetree/event tree
// found on file: RPCv4t_evtraw-02092017-122902.rre
//////////////////////////////////////////////////////////

#ifndef Selector_h
#define Selector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

#include <vector>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TGraph2DErrors.h"
#include "TTree.h"
#include "TProfile.h"
#include "TString.h"
#include <TProofOutputFile.h>
#include "Constants.h"

// Headers needed by this particular selector
#include "TTimeStamp.h"
#include "TBits.h"



class   TH1F;

class Selector : public TSelector {
public :
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           ENum[nlayer];
   Int_t           REnum[nlayer];
   ULong64_t       CEnum;
   TTimeStamp *EveTS[nlayer];
   TBits *xLayer[nlayer];
   TBits *yLayer[nlayer];

   Int_t           tdc_ref_l[nlayer];
   Int_t           tdc_ref_t[nlayer];
   Int_t           trigCntDiff[nlayer];

   std::vector<unsigned int> *vxtdc_l[nlayer][nchannel];
   std::vector<unsigned int> *vytdc_l[nlayer][nchannel];

   std::vector<unsigned int> *vxtdc_t[nlayer][nchannel];
   std::vector<unsigned int> *vytdc_t[nlayer][nchannel];


   // List of branches
   TBranch        *b_ENum;   //!
   TBranch        *b_REnum;   //!
   TBranch        *b_CEnum;   //!
   TBranch        *b_Evetime[nlayer];   //!
   TBranch        *b_xstriphitsL[nlayer];   //!
   TBranch        *b_ystriphitsL[nlayer];   //!

   TBranch        *b_tdc_ref_l;   //!
   TBranch        *b_tdc_ref_t;   //!
   TBranch        *b_trigCntDiff;   //!

   TBranch        *b_xtdc_l[nlayer][nchannel];   //!
   TBranch        *b_ytdc_l[nlayer][nchannel];   //!

   TBranch        *b_xtdc_t[nlayer][nchannel];   //!
   TBranch        *b_ytdc_t[nlayer][nchannel];   //!

   // Specific members
   TFile            *fFile;
   TProofOutputFile *fProofFile; // For optimized merging of the ntuple

   TH1F* xlayer_occu[nlayer];
   TH1F* ylayer_occu[nlayer];
   TH1F* xlayer_mult[nlayer];
   TH1F* ylayer_mult[nlayer];
   TH2F* rawhits_corr_xymul[nlayer];
   TH2F* raw_occu[nlayer];
   TH2F* rawhits_xlay_corr_mul[nlayer][nlayer];
   TH2F* rawhits_ylay_corr_mul[nlayer][nlayer];
   TH1F* xlayer_reso[nlayer];
   TH1F* ylayer_reso[nlayer];
   TH2D* totalentry[nlayer];

   TH1F* xtime_layer_1hit[nlayer];
   TH1F* ytime_layer_1hit[nlayer];
   TH1F* xtime_layer_hough[nlayer];
   TH1F* ytime_layer_hough[nlayer];
   TH1F* xtime_layer[nlayer];
   TH1F* ytime_layer[nlayer];

   TH1F* distance_from_pos;
   TH1F* distance_from_ext;

   TH1F* distance_from_pos_l[nlayer];
   TH1F* distance_from_ext_l[nlayer];

   TH1F *h1;//!

   Selector(TTree *tree =0):fFile(0), fProofFile(0){ }
   virtual ~Selector();
   void    Reset();

   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();
   //Custom Defined Functions
   Double_t cal_slope2(Double_t x, Double_t* par);
   void GetPosInStrip(int ixy, double* ext, double* otherext, double* off, double* pos, double* local);

   ClassDef(Selector,0);


private:

  vector<int>deadchannel[2][nlayer];

};

#endif

#ifdef Selector_cxx
void Selector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   const char *sideMark[2] = {"x","y"};
   for(int i=0;i<nlayer;i++){
     for(int j=0;j<nchannel;j++){
       vxtdc_l[i][j]=0;
       vytdc_l[i][j]=0;
       vxtdc_t[i][j]=0;
       vytdc_t[i][j]=0;
     }
   }

   for (int ixy=0; ixy<2; ixy++) {
     for (int il=0; il<nlayer; il++) {
       deadchannel[ixy][il].clear();
     }
   }

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ENum", ENum, &b_ENum);
   fChain->SetBranchAddress("REnum", REnum, &b_REnum);
   fChain->SetBranchAddress("CEnum", &CEnum, &b_CEnum);

   fChain->SetBranchAddress("tdc_ref_l", tdc_ref_l, &b_tdc_ref_l);
   fChain->SetBranchAddress("tdc_ref_t", tdc_ref_t, &b_tdc_ref_t);
   fChain->SetBranchAddress("trigCntDiff", trigCntDiff, &b_trigCntDiff);
   for(int ij=0;ij<nlayer;ij++) {
     EveTS[ij] = 0;
     xLayer[ij] = 0;
     yLayer[ij] = 0;
     fChain->SetBranchAddress(TString::Format("Evetime_%i",ij), &EveTS[ij], &b_Evetime[ij]);
     fChain->SetBranchAddress(TString::Format("xstriphitsL%i",ij), &xLayer[ij], &b_xstriphitsL[ij]);
     fChain->SetBranchAddress(TString::Format("ystriphitsL%i",ij), &yLayer[ij], &b_ystriphitsL[ij]);
     for(int jk=0;jk<nchannel;jk++) {
       fChain->SetBranchAddress(TString::Format("xtdc_l_%i_%i",ij,jk), &vxtdc_l[ij][jk], &b_xtdc_l[ij][jk]);
       fChain->SetBranchAddress(TString::Format("ytdc_l_%i_%i",ij,jk), &vytdc_l[ij][jk], &b_ytdc_l[ij][jk]);

       fChain->SetBranchAddress(TString::Format("xtdc_t_%i_%i",ij,jk), &vxtdc_t[ij][jk], &b_xtdc_t[ij][jk]);
       fChain->SetBranchAddress(TString::Format("ytdc_t_%i_%i",ij,jk), &vytdc_t[ij][jk], &b_ytdc_t[ij][jk]);
     }
   }

}

Bool_t Selector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Selector::Reset()
{
   // Reset the data members to theit initial value

   h1 = 0;
   //fChain = 0;
   //elist = 0;
   //fillList = kFALSE;
   //useList  = kFALSE;
   //fProcessed = 0;
}

#endif // #ifdef Selector_cxx
