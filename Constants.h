#ifndef Constants_h
#define Constants_h

const double cval = 29.979; // velocity of light in cm/ns
static const int  nlayer =12; //Maximum RPC layers
static const int  nchannel=8;  //Number of TDC Channels
const int  nstrip = 64; // Strips in RPC chamber
const int  nusedstrip =64; // Strips connected with Amplifier/ADC in RPC chamber
const int  lastXstrip =58; // last strip in X
const int  lastYstrip =61; // last strip in Y

bool isOnlyCom  = true; //Combined fit only, no iteration, not proper efficiency
bool isTimeCorOrReso = true; // 10th Aug 2016, must be true always //false; //true; // true: Time offset Correction and
                                    // false: Time resolution, they should come one by one


const double  mxchisq =2.0;//15;//2.0;   //Maximum value of Normalised chi^2 (position fit);
const double  mxtimechisq=2.0; // Maximum value of Normalised chi^2 (time fit);
const double  accprng =1.0; //Acceptance region with respect tot strip centre
const double  effirng =0.25;//0.25;//0.25; //Additional region for efficeincy calculation
                           //Extrapolation error < 0.2 of strip width
const double pival=acos(-1);

const int  nTDCpLayer=8;
const int nmxhits=4;
const int nmnhits =4;//Minimum number of layers required for position and time measurements--Must greater than two
const int nmxusedhits=3;
const int nmxtimehit=4;
const int nmxusedtimehit=3;
const int nPosAlignPar=3;
//Used for fitting
const double layerzpos[nlayer] = {0.0,10.4,20.2,30.2,40.7,51.0,61.1,71.3,81.3,91.9,101.9,102.9};//9.6,19.2,28.8,38.4,48.0,57.6,67.2,76.8,86.4,96.0,105.6};
const int  layfirst =0;  //Used first layer in track fitting
const int  laylast =11;  //Used last layer in track fitting
const float xyPosDev = 7.0; // seven sigma 2.0; //maximum deviation of points from fit line

//Trigger layers
const int trigly1 =6;
const int trigly2 =7;
const int trigly3 =8;
const int trigly4 =9;

//Physical information
const double stripwidth = 3.0; // cm
const double stripgap = 0.2; // cm
const double layergap = 16.0; // cm

//corrections from RPCv4t_evtraw-20181226_192924_miical_jj14_3.root
double align_xstr_ydev[nlayer][nPosAlignPar] ={{0.893067, -0.001879, 0.900119},
{0.562203, -0.00614427, 0.439102},
{0.600602, -0.000126263, 0.651504},
{0.568652, -0.00510885, 0.410989},
{-0.00312035, 0.00782, 0.258904},
{-0.0288475, 0.00535717, 0.170647},
{-0.14004, 0.000312532, -0.113553},
{-0.148819, 0.00241532, -0.0323629},
{0.0373388, 0.00676629, 0.246137},
{0.0321182, -0.00204897, 0.00499599},
{0.117565, 0.00409309, 0.235423},
{-0.0716595, 0.000813741, -0.0348974}};
double align_ystr_xdev[nlayer][nPosAlignPar] ={{-0.371309, -0.00241894, -0.481752},
{-0.158913, 0.000610099, -0.175186},
{-0.389387, -0.00418929, -0.449643},
{-0.0828907, 0.00487446, 0.0366989},
{-0.104944, -0.00410707, -0.278469},
{0.131809, -0.00417544, -0.0118498},
{0.164401, -0.00209133, 0.0363426},
{0.192258, -0.00121329, 0.122072},
{0.159859, -0.00544231, -0.0387917},
{0.183074, 0.0015277, 0.133429},
{-0.0991423, -0.00246185, -0.205385},
{-0.215174, -0.000523797, -0.279876}};
double align_xstr_xdev[nlayer][nPosAlignPar] ={{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0}};
double align_ystr_ydev[nlayer][nPosAlignPar] ={{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0},
{0, 0, 0}};


double xoff[nlayer] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
double yoff[nlayer] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

double xposerrsq[nmxhits][nlayer] = {{0.115948, 0.075825, 0.0848383, 0.0646547, 0.0670733, 0.0701558, 0.0646064, 0.0591128, 0.0727697, 0.114561, 0.02, 0.02},
{0.0911147, 0.0399934, 0.0318722, 0.0358208, 0.0394841, 0.0301098, 0.045686, 0.02, 0.0595045, 0.10909, 0.02, 0.02},
{0.168243, 0.132579, 0.192114, 0.135142, 0.113249, 0.195536, 0.139515, 0.140938, 0.1433, 0.195018, 0.02, 0.02},
{0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02}};
double yposerrsq[nmxhits][nlayer] = {{0.110846, 0.0767004, 0.0640187, 0.0615692, 0.0579787, 0.0634886, 0.056111, 0.0518296, 0.0561367, 0.0791761, 0.02, 0.02},
{0.0865317, 0.0467293, 0.0378272, 0.0389136, 0.0484361, 0.0425223, 0.0435002, 0.0398323, 0.0396466, 0.0704588, 0.02, 0.02},
{0.200208, 0.147808, 0.158143, 0.150396, 0.121554, 0.165004, 0.115626, 0.14542, 0.13671, 0.176862, 0.02, 0.02},
{0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02}};


float errxco[nlayer]={0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675};
float erryco[nlayer]={0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675};




#endif
