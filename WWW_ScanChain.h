//  .
// ..: P. Chang, philip@physics.ucsd.edu

// RooUtil tool
#include "rooutil/looper.cc"
#include "rooutil/autohist.cc"
#include "rooutil/eventlist.cc"
#include "WWW_CORE/WWWTree.h"
#include "WWW_CORE/WWWTools.h"

#include "CORE/Tools/dorky/dorky.h"

using namespace std;
using namespace RooUtil;

#define MAXOBJ 3
#define NSYST 3

void ScanChain( TChain* chain, TString output_name, TString optstr, int nevents = -1 ); // the default nevents=-1 option means loop over all events.

bool doAnalysis( RooUtil::AutoHist& );

void fillHistograms( bool& passed, RooUtil::AutoHist&, TString, int regionid, int isyst );
void fillHistogramsFull( RooUtil::AutoHist&, TString, TString, TString, int regionid, int isyst );
void fillLepHistograms( RooUtil::AutoHist&, TString, TString, TString );
void fillJetHistograms( RooUtil::AutoHist&, TString, TString, TString );
void fillWWWHistograms( RooUtil::AutoHist&, TString );

void printevent( TString );

//void doTmpAnalysis( RooUtil::AutoHist& hists );
