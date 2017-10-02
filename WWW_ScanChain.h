//  .
// ..: P. Chang, philip@physics.ucsd.edu

// RooUtil tool
#include "rooutil/looper.h"
#include "rooutil/autohist.h"
#include "rooutil/eventlist.h"
#include "WWW_CORE/WWWTree.h"
#include "WWW_CORE/WWWTools.h"

#include "CORE/Tools/dorky/dorky.h"

using namespace std;
using namespace RooUtil;

#define MAXOBJ 3
#define NSYST 3

void ScanChain( TChain* chain, TString output_name, TString optstr, int nevents = -1 ); // the default nevents=-1 option means loop over all events.

bool doAnalysis( RooUtil::AutoHist&, bool doskim );

void fillHistograms( RooUtil::AutoHist&, TString, int regionid, int isyst );
void fillHistogramsFull( RooUtil::AutoHist&, TString, TString, TString, int regionid, int isyst );
void fillLepHistograms( RooUtil::AutoHist&, TString, TString, TString, int );
void fillJetHistograms( RooUtil::AutoHist&, TString, TString, TString, int );
void fillWWWHistograms( RooUtil::AutoHist&, TString, int );

void printevent( TString );

//void doTmpAnalysis( RooUtil::AutoHist& hists );
