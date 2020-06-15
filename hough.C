
#include "TInterpreter.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TFrame.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include <TH2F.h>
#include <TMath.h>

#include "vector"

void hough()
{
  printf("Program to find straightline fit using hough transform\n");
  int nlayer=12;
  int nchannel=8;
  //Layer Info
  double z[nlayer];
  for(int i=0;i<nlayer;i++)
  {
    z[i]=double(i);
  }
  std::vector<unsigned int> vxtdc_l[nlayer][nchannel];
  // for(int i=0;i<nlayer;i++){
  //   for(int j=0;j<nchannel;j++){
  //     vxtdc_l[i][j]=0;
  //   }
  // }

  //TDC values

  vxtdc_l[0][2].push_back(100);
  vxtdc_l[0][2].push_back(20);
  vxtdc_l[1][1].push_back(30);
  vxtdc_l[2][2].push_back(40);
  vxtdc_l[3][6].push_back(50);
  vxtdc_l[4][0].push_back(60);
  vxtdc_l[5][0].push_back(70);
  vxtdc_l[6][0].push_back(80);
  vxtdc_l[7][0].push_back(90);
  vxtdc_l[8][0].push_back(100);
  vxtdc_l[9][0].push_back(110);


  double rmin,rmax,thetamin,thetamax;
  int nbinr,nbintheta;

  rmin=-100.;
  rmax=100.;
  thetamin=0.;
  thetamax=180.;
  nbinr=1000;
  nbintheta=1000;

  gStyle->SetOptFit();
  TCanvas *c1 = new TCanvas("c1","multigraph",700,500);
  c1->SetGrid();
  c1->Divide(2,1);

   // draw a frame to define the range
   TMultiGraph *mg = new TMultiGraph();

  TH2F *hough_space = new TH2F("hough_space","r vs theta",nbinr,rmin,rmax,nbintheta,thetamin,thetamax);
  TGraph *graph1 = new TGraph();
  TGraph *graph2 = new TGraph();

  double r;
  double x,y;
  int i=0;
  for(int ij=0;ij<nlayer;ij++)
  {
    for(int jk=0;jk<nchannel;jk++)
    {
      for(int tsize=0;tsize<vxtdc_l[ij][jk].size();tsize++)
      {
        y=(double)ij;
        x=(double)vxtdc_l[ij][jk].at(tsize);
        graph1->SetPoint(i,x,y);
        i=i+1;
        printf("x=%f, y=%f\n",x,y);
        for(double theta=thetamin;theta<=thetamax;theta=theta+180./1000.)
          {
            r=x*cos(TMath::DegToRad()*theta)+y*sin(TMath::DegToRad()*theta);
            //printf("r=%f, theta=%f\n",r,theta);
            hough_space->Fill(r,theta);
          }
      }
    }
  }
  c1->cd(1);
  printf("%d\n", hough_space->GetMaximumBin());
  Int_t rloc,thetaloc,zloc;
  hough_space->GetBinXYZ(hough_space->GetMaximumBin(), rloc, thetaloc, zloc);
  printf("The bin having the maximum value is (%d,%d)\n",rloc,thetaloc);
  double rval=rmin+(rmax-rmin)*((double)rloc-0.5)/nbinr;
  double thetaval=thetamin+(thetamax-thetamin)*((double)thetaloc-0.5)/nbintheta;
  printf("The r,theta having the maximum value is (%f,%f)\n",rval,thetaval);

  //hough_space->Reset("ICESM");
  hough_space->Draw("colz");

  graph1->GetXaxis()->SetTitle("time");
  graph1->GetYaxis()->SetTitle("z");
  graph1->GetXaxis()->CenterTitle();
  graph1->GetYaxis()->CenterTitle();
  graph1->GetXaxis()->SetTitleFont(62);
  graph1->GetYaxis()->SetTitleFont(62);
  graph1->GetXaxis()->SetLabelFont(62);
  graph1->GetYaxis()->SetLabelFont(62);
  graph1->SetTitle();
  graph1->SetMarkerColor(kBlue);
  graph1->SetMarkerSize(1);
  graph1->SetMarkerStyle(20);
  //graph1->Draw("AP");

  i=0;
  for(int ij=0;ij<nlayer;ij++)
  {
    y=ij;
    x=(rval-(y*sin(TMath::DegToRad()*thetaval)))/cos(TMath::DegToRad()*thetaval);
    graph2->SetPoint(i,x,y);
    i=i+1;
  }

  graph2->GetXaxis()->SetTitle("time");
  graph2->GetYaxis()->SetTitle("z");
  graph2->GetXaxis()->CenterTitle();
  graph2->GetYaxis()->CenterTitle();
  graph2->GetXaxis()->SetTitleFont(62);
  graph2->GetYaxis()->SetTitleFont(62);
  graph2->GetXaxis()->SetLabelFont(62);
  graph2->GetYaxis()->SetLabelFont(62);
  graph2->SetTitle();
  graph2->SetMarkerColor(kGreen);
  graph2->SetMarkerSize(1);
  graph2->SetMarkerStyle(20);
  //graph2->Draw("AP");

  c1->cd(2);
  mg->Add(graph1);
  mg->Add(graph2);
  mg->Draw("ap");

}
