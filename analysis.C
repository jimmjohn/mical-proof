/// This code uses an existing PROOF session or starts one at the indicated URL.
/// In the case non existing PROOF session is found and no URL is given, the
/// macro tries to start a local PROOF session.


#include <iostream>
#include <functional>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <new>
#include <climits>
#include <vector>

//Root classes
#include "TCanvas.h"
#include "TChain.h"
#include "TDSet.h"
#include "TEnv.h"
#include "TEntryList.h"
#include "TFile.h"
#include "TFileCollection.h"
#include "TFrame.h"
#include "THashList.h"
#include "TList.h"
#include "TMath.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TProof.h"
#include "TProofDebug.h"
#include "TString.h"
#include <TSystem.h>

//Custom made classes
#include "Constants.h"        // All the constants are put inside this file
#include "getProof.C"

void SavePerfTree(TProof *proof, const char *fn);

using namespace std;

//void analysis(const char *what = "simple",
//              const char *masterurl = "proof://localhost:40000",
//              Int_t nwrks = -1, TList *ins = 0)
void analysis()
{
  const char *what = "simple(nevt=1000000)";
  const char *masterurl = "lite://";
  Int_t nwrks = 6;
  TList *ins = 0;
  gEnv->SetValue("Proof.StatsHist",1);

  if (nmnhits <=2) { //GMA
    cerr<<"nmnhits must be more than two, but value is ,"<<nmnhits <<",Change this value and rerun the code"<<endl;
    return;
  }

  TString u(masterurl);
  // Determine locality of this session
  Bool_t isProofLocal = kFALSE;
  if (!u.IsNull() && u != "lite://") {
     TUrl uu(masterurl);
     TString uopts(uu.GetOptions());
     if ((!strcmp(uu.GetHost(), "localhost") && !uopts.Contains("external")) ||
        !strcmp(uu.GetHostFQDN(), TUrl(gSystem->HostName()).GetHostFQDN())) {
        isProofLocal = kTRUE;
     }
     // Adjust URL
     if (!u.BeginsWith(uu.GetProtocol())) uu.SetProtocol("proof");
     uopts.ReplaceAll("external", "");
     uu.SetOptions(uopts.Data());
     u = uu.GetUrl();
  }
  const char *url = u.Data();

  // Temp dir for PROOF work
  // Force "/tmp/<user>" whenever possible to avoid length problems on MacOsX
  TString tmpdir("/tmp");
  if (gSystem->AccessPathName(tmpdir, kWritePermission)) tmpdir = gSystem->TempDirectory();
  TString us;
  UserGroup_t *ug = gSystem->GetUserInfo(gSystem->GetUid());
  if (!ug) {
     Printf("runProof: could not get user info");
     return 0;
  }
  us.Form("/%s", ug->fUser.Data());
  if (!tmpdir.EndsWith(us.Data())) tmpdir += us;
  gSystem->mkdir(tmpdir.Data(), kTRUE);
  if (gSystem->AccessPathName(tmpdir, kWritePermission)) {
     Printf("runProof: unable to get a writable tutorial directory (tried: %s)"
            " - cannot continue", tmpdir.Data());
     return 0;
  }
  TString workdir = Form("%s/.proof-work", tmpdir.Data());
  if (gSystem->AccessPathName(workdir)) {
     Printf("runProof: creating the temporary directory"
               " for the work (%s) ... ", workdir.Data());
     if (gSystem->mkdir(workdir, kTRUE) != 0) {
        Printf("runProof: could not assert / create the temporary directory"
               " for the work (%s)", workdir.Data());
        return 0;
     }
  }

  Printf("work dir:\t%s", workdir.Data());

  // Get the PROOF Session
  getProof *goProof = new getProof();
  TProof *proof = goProof->getProofSection(url, nwrks, workdir.Data(), "ask",kFALSE,kFALSE);
  if (!proof) {
     Printf("runProof: could not start/attach a PROOF session");
     return 0;
  }

  // Refine locality (PROOF-Lite always local)
  if (proof->IsLite()) isProofLocal = kTRUE;

  TString proofsessions(Form("%s/sessions",workdir.Data()));
  // Save tag of the used session
  FILE *fs = fopen(proofsessions.Data(), "a");
  if (!fs) {
     Printf("runProof: could not create files for sessions tags");
  } else {
     fprintf(fs,"session-%s\n", proof->GetSessionTag());
     fclose(fs);
  }
  if (!proof) {
     Printf("runProof: could not start/attach a PROOF session");
     return 0;
  }

  proof->Load("Constants.h");
  proof->Load("StraightLineFit.C+");


  // Set the number of workers (may only reduce the number of active workers
  // in the session)
  if (nwrks > 0)
     proof->SetParallel(nwrks);

  // Parse 'what'; it is in the form 'analysis(arg1,arg2,...)'
  TString args(what);
  args.ReplaceAll("("," ");
  args.ReplaceAll(")"," ");
  args.ReplaceAll(","," ");
  Ssiz_t from = 0;
  TString act, tok;
  if (!args.Tokenize(act, from, " ")) {
    // Cannot continue
    Printf("runProof: action not found: check your arguments (%s)", what);
    return 0;
  }
  // Extract ACLiC mode
  TString aMode = "+";
  if (act.EndsWith("+")) {
     aMode += "+";
     while (act.EndsWith("+")) { act.Remove(TString::kTrailing,'+'); }
  }
  Printf("runProof: %s: ACLiC mode: '%s'", act.Data(), aMode.Data());

  // Parse out number of events and  'asyn' option, used almost by every test
  TString aNevt, aFirst, aNwrk, opt, sel, punzip("on"), aCache, aOutFile,
          aDebug, aDebugEnum, aRateEst, aPerfTree("perftree.root"),
          aFeedback("fb=stats");
  Long64_t suf = 1;
  Int_t aSubMg = -1;
  Bool_t makePerfTree = kFALSE;
  while (args.Tokenize(tok, from, " ")) {
     // Debug controllers
     if (tok.BeginsWith("debug=")) {
        aDebug = tok;
        aDebug.ReplaceAll("debug=","");
        Int_t icol = kNPOS;
        if ((icol = aDebug.Index(":")) != kNPOS) {
           aDebugEnum = aDebug(0, icol);
           aDebug.Remove(0, icol+1);
        }
        if (!aDebug.IsDigit()) {
           Printf("runProof: %s: error parsing the 'debug=' option (%s) - ignoring", act.Data(), tok.Data());
           aDebug = "";
           aDebugEnum = "";
        }
     }
     // Number of events
     if (tok.BeginsWith("nevt=")) {
        aNevt = tok;
        aNevt.ReplaceAll("nevt=","");
        if (!aNevt.IsDigit()) {
           Printf("runProof: %s: error parsing the 'nevt=' option (%s) - ignoring", act.Data(), tok.Data());
           aNevt = "";
        }
     }
     // First event
     if (tok.BeginsWith("first=")) {
        aFirst = tok;
        aFirst.ReplaceAll("first=","");
        if (!aFirst.IsDigit()) {
           Printf("runProof: %s: error parsing the 'first=' option (%s) - ignoring", act.Data(), tok.Data());
           aFirst = "";
        }
     }
     // Sync or async ?
     if (tok.BeginsWith("asyn"))
        opt = "ASYN";
     // Number of workers
     if (tok.BeginsWith("nwrk=")) {
        aNwrk = tok;
        aNwrk.ReplaceAll("nwrk=","");
        if (!aNwrk.IsDigit()) {
           Printf("runProof: %s: error parsing the 'nwrk=' option (%s) - ignoring", act.Data(), tok.Data());
           aNwrk = "";
        }
     }
     // Parallel unzipping ?
     if (tok.BeginsWith("punzip"))
        punzip = "on";
     // Number of workers
     if (tok.BeginsWith("cache=")) {
        aCache = tok;
        aCache.ReplaceAll("cache=","");
        if (aCache.EndsWith("k")) { aCache.Remove(TString::kTrailing, 'k'); suf = 1024; }
        if (aCache.EndsWith("K")) { aCache.Remove(TString::kTrailing, 'K'); suf = 1024; }
        if (aCache.EndsWith("M")) { aCache.Remove(TString::kTrailing, 'M'); suf = 1024*1024; }
        if (!aCache.IsDigit()) {
           Printf("runProof: %s: error parsing the 'cache=' option (%s) - ignoring", act.Data(), tok.Data());
           aCache = "";
        }
     }
     // Use submergers?
     if (tok.BeginsWith("submergers")) {
        tok.ReplaceAll("submergers","");
        aSubMg = 0;
        if (tok.BeginsWith("=")) {
           tok.ReplaceAll("=","");
           if (tok.IsDigit()) aSubMg = tok.Atoi();
        }
     }

     // Rate estimation technique
     if (tok.BeginsWith("rateest=")) {
        tok.ReplaceAll("rateest=","");
        if (!(tok.IsNull())) aRateEst = tok;
        Printf("runProof: %s: progress-bar rate estimation option: '%s'", act.Data(), aRateEst.Data());
     }
     // Create and save the preformance tree?
     if (tok.BeginsWith("perftree")) {
        makePerfTree = kTRUE;
        if (tok.BeginsWith("perftree=")) {
           tok.ReplaceAll("perftree=","");
           if (!(tok.IsNull())) aPerfTree = tok;
        }
        Printf("runProof: %s: saving performance tree to '%s'", act.Data(), aPerfTree.Data());
     }
     // Location of the output file, if any
     if (tok.BeginsWith("outfile")) {
        if (tok.BeginsWith("outfile=")) {
           tok.ReplaceAll("outfile=","");
           if (!(tok.IsNull())) aOutFile = tok;
        }
        Printf("runProof: %s: output file: '%s'", act.Data(), aOutFile.Data());
     }
     // Feedback
     if (tok.BeginsWith("feedback=")) {
        tok.ReplaceAll("feedback=","");
        if (tok == "off" || tok == "OFF" || tok == "0") {
           aFeedback = "";
        } else if (!(tok.IsNull())) {
           if (tok.BeginsWith("+")) {
              tok[0] = ',';
              aFeedback += tok;
           } else {
              aFeedback.Form("fb=%s", tok.Data());
           }
        }
        Printf("runProof: %s: feedback: '%s'", act.Data(), aFeedback.Data());
     }
  }
  Long64_t nevt = (aNevt.IsNull()) ? -1 : aNevt.Atoi();
  Long64_t first = (aFirst.IsNull()) ? 0 : aFirst.Atoi();
  Long64_t nwrk = (aNwrk.IsNull()) ? -1 : aNwrk.Atoi();
  from = 0;

  // Set number workers
  if (nwrk > 0) {
     if (proof->GetParallel() < nwrk) {
        Printf("runProof: %s: request for a number of workers larger then available - ignored", act.Data());
     } else {
        proof->SetParallel(nwrk);
     }
  }

  // Debug controllers
  if (!aDebug.IsNull()) {
     Int_t dbg = aDebug.Atoi();
     Int_t scope = TProofDebug::kAll;
     if (!aDebugEnum.IsNull()) scope = goProof->getDebugEnum(aDebugEnum.Data());
     proof->SetLogLevel(dbg, scope);
     Printf("runProof: %s: verbose mode for '%s'; level: %d", act.Data(), aDebugEnum.Data(), dbg);
  }

  // Have constant progress reporting based on estimated info
  // (NB: may screw up the progress bar in some cases)
  if (aRateEst == "average")
     proof->SetParameter("PROOF_RateEstimation", aRateEst);

  // Parallel unzip
  if (punzip == "on") {
     proof->SetParameter("PROOF_UseParallelUnzip", (Int_t)1);
     Printf("runProof: %s: parallel unzip enabled", act.Data());
  } else {
     proof->SetParameter("PROOF_UseParallelUnzip", (Int_t)0);
  }

  // Tree cache
  if (!aCache.IsNull()) {
     Long64_t cachesz = aCache.Atoi() * suf;
     if (cachesz <= 0) {
        proof->SetParameter("PROOF_UseTreeCache", (Int_t)0);
        Printf("runProof: %s: disabling tree cache", act.Data());
     } else {
        proof->SetParameter("PROOF_UseTreeCache", (Int_t)1);
        proof->SetParameter("PROOF_CacheSize", cachesz);
        Printf("runProof: %s: setting cache size to %lld", act.Data(), cachesz);
     }
  } else {
     // Use defaults
     proof->DeleteParameters("PROOF_UseTreeCache");
     proof->DeleteParameters("PROOF_CacheSize");
  }

  // Enable submergers, if required
  if (aSubMg >= 0) {
     proof->SetParameter("PROOF_UseMergers", aSubMg);
     if (aSubMg > 0) {
        Printf("runProof: %s: enabling merging via %d sub-mergers", act.Data(), aSubMg);
     } else {
        Printf("runProof: %s: enabling merging via sub-mergers (optimal number)", act.Data());
     }
  } else {
     proof->DeleteParameters("PROOF_UseMergers");
  }

  // The performance tree
  if (makePerfTree) {
     proof->SetParameter("PROOF_StatsHist", "");
     proof->SetParameter("PROOF_StatsTrace", "");
     proof->SetParameter("PROOF_SlaveStatsTrace", "");
  }

  // Additional inputs from the argument 'ins'
  if (ins && ins->GetSize() > 0) {
     TObject *oin = 0;
     TIter nxo(ins);
     while ((oin = nxo())) { proof->AddInput(oin); }
  }
  // Full lits of inputs so far
  proof->GetInputList()->Print();

  // Action
  if (act == "simple") {
    TString rootfiles("data.log");
    TString datafile;
    TString outfile;
    TString outfilx;
    TString infile;
    ifstream file_db;
    file_db.open(rootfiles);
    // Create the chain
    TChain *chain = new TChain("evetree");
    while(!(file_db.eof())) {
      file_db >> datafile;           //Each root file name
      if(file_db.eof()) break;
      Printf("datafile : %s",datafile.Data());
      outfilx = datafile(0, datafile.Length()-4);
      outfile.Form("%s.root",outfilx.Data());
      infile.Form("%s/data/%s",gSystem->WorkingDirectory(),datafile.Data());
      chain->Add(TString::Format("%s", infile.Data()));
      Printf("infile : %s",infile.Data());
    }

    // Output file
    TString fout(aOutFile);
    if (fout.IsNull()) {
       fout.Form("%s/ProofOutput.root", gSystem->WorkingDirectory());
       // Cleanup any existing instance of the output file
       gSystem->Unlink(fout);

       if (!isProofLocal) {
          // Setup a local basic xrootd to receive the file
          Bool_t xrdok = kFALSE;
          Int_t port = 9000;
          while (port < 9010) {
             if (goProof->checkXrootdAt(port) != 1) {
                if (goProof->startXrootdAt(port, gSystem->WorkingDirectory(), kTRUE) == 0) {
                   xrdok = kTRUE;
                   break;
                }
             }
             port++;
          }
          if (!xrdok) {
             Printf("runProof: could not start basic xrootd on ports 9000-9009 - cannot continue");
             return;
          }
          fout.Insert(0, TString::Format("root://%s:%d/", TUrl(gSystem->HostName()).GetHostFQDN(), port));
          // Make a copy of the files on the master before merging
          proof->AddInput(new TNamed("PROOF_OUTPUTFILE_LOCATION", "LOCAL"));
       }
    }
    proof->AddInput(new TNamed("PROOF_OUTPUTFILE", fout.Data()));


    chain->ls();
    // We run on Proof
    chain->SetProof();
    // The selector
    sel.Form("Selector.C%s", aMode.Data());
    // Run it
    Printf("\nrunProof: running \"Chain\"\n");
    TString xopt = aFeedback; if (!opt.IsNull()) xopt += TString::Format(" %s", opt.Data());
    chain->Process(sel.Data(),xopt,nevt,first);
  }

/*
  [out=std::ref(std::cout << "Pressure")](){out.get() << " Coefficient Project\n"; }();
  RootFileStruct *pAnalysis = new RootFileStruct();
  RunAction *runAct = new RunAction();
  runAct->BeginOfRunAction();
  runAct->EndOfRunAction();
*/
// Save the performance tree
if (makePerfTree) {
   SavePerfTree(proof, aPerfTree.Data());
   // Cleanup parameters
   gProof->DeleteParameters("PROOF_StatsHist");
   gProof->DeleteParameters("PROOF_StatsTrace");
   gProof->DeleteParameters("PROOF_SlaveStatsTrace");
}

  return 0;
}

//______________________________________________________________________________
void SavePerfTree(TProof *proof, const char *fn)
{
   // Save PROOF timing information from TPerfStats to file 'fn'

   if (!proof) {
      Printf("PROOF must be run to save output performance information");;
      return;
   }
   if (!proof->GetOutputList() || proof->GetOutputList()->GetSize() <= 0) {
      Printf("PROOF outputlist undefined or empty");;
      return;
   }

   TFile f(fn, "RECREATE");
   if (f.IsZombie()) {
      Printf("ERROR: could not open file '%s' for writing", fn);;
   } else {
      f.cd();
      TIter nxo(proof->GetOutputList());
      TObject* obj = 0;
      while ((obj = nxo())) {
         TString objname(obj->GetName());
         if (objname.BeginsWith("PROOF_")) {
            // Must list the objects since other PROOF_ objects exist
            // besides timing objects
            if (objname == "PROOF_PerfStats" ||
                objname == "PROOF_PacketsHist" ||
                objname == "PROOF_EventsHist" ||
                objname == "PROOF_NodeHist" ||
                objname == "PROOF_LatencyHist" ||
                objname == "PROOF_ProcTimeHist" ||
                objname == "PROOF_CpuTimeHist")
               obj->Write();
         }
      }
      f.Close();
   }

}
