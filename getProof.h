#ifndef getProof_h
#define getProof_h 1


#include "Bytes.h"
#include "Getline.h"
#include "TEnv.h"
#include "TProof.h"
#include "TSocket.h"
#include "TString.h"
#include "TSystem.h"


class getProof
{
public:
  getProof();
  ~getProof();

public:
  TProof* getProofSection(const char *url,Int_t nwrks,const char *dir,const char *opt,Bool_t dyn,Bool_t tutords);

  // Auxilliary functions
  int getDebugEnum(const char *what);
  Int_t getXrootdPid(Int_t port, const char *subdir = "xpdtut");
  Int_t checkXrootdAt(Int_t port, const char *host = "localhost");
  Int_t checkXproofdAt(Int_t port, const char *host = "localhost");
  Int_t startXrootdAt(Int_t port, const char *exportdirs = 0, Bool_t force = kFALSE);
  Int_t killXrootdAt(Int_t port, const char *id = 0);

  // Auxilliary structures for Xrootd/Xproofd pinging ...
  // The client request
  typedef struct {
    int first;
    int second;
    int third;
    int fourth;
    int fifth;
  } clnt_HS_t;
  // The body received after the first handshake's header
  typedef struct {
    int msglen;
    int protover;
    int msgval;
  } srv_HS_t;

  // By default we start a cluster on the local machine
  const char *refloc = "proof://localhost:40000";

private:

};

#endif
