//Generalllllll

gStyle->SetPalette(1,0);
gStyle->SetFillColor(0);
gStyle->SetCanvasBorderMode(0);
gStyle->SetPadBorderMode(0);
gStyle->SetStatBorderSize(1);
gStyle->SetStatStyle(1001);
gStyle->SetCanvasColor(10);
gStyle->SetPadColor(10);
gStyle->SetStatColor(10);
gStyle->SetTitleFillColor(10);
gStyle->SetTitleBorderSize(1);
gStyle->SetTitleYOffset(0.05);

gStyle->SetStatFont(22);        // Times New Roman
gStyle->SetTextFont(22);        // Times New Roman
gStyle->SetTitleFont(22,"XYZ"); // Times New Roman
gStyle->SetLabelFont(22,"XYZ"); // Times New Roman
gStyle->SetLabelSize(0.07, "XYZ"); // Times New Roman
gStyle->SetNdivisions(606, "XYZ");
gStyle->SetPaintTextFormat("6.4f");


gStyle->SetOptTitle(0);
gStyle->SetFuncWidth(1);
gStyle->SetFuncColor(2);
gStyle->SetOptStat(1110);
gStyle->SetOptFit(101);
gStyle->SetOptLogy(0);
gStyle->SetStatW(.18);
gStyle->SetStatH(.08);
gStyle->SetPadTopMargin(.001); //0.09
gStyle->SetPadBottomMargin(0.08);
gStyle->SetPadLeftMargin(0.001);
gStyle->SetPadRightMargin(0.001);


TLatex latex;
latex.SetNDC();
latex.SetTextSize(0.12);
latex.SetTextFont(42);
latex.SetTextAlign(1); //(31); // align right
