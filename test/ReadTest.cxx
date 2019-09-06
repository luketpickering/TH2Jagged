
#define TH2JAGGED_DEBUG
#include "TH2Jagged.h"

#include "TFile.h"
#include "TRandom3.h"
#include "TH2Poly.h"

#include <iostream>

int main() {

  TH2JaggedD *jag = 0;
  TFile fr("TH2JagWriteTest.root", "READ");
  fr.GetObject("jag", jag);
  assert(jag);

  assert(jag->fBinMappingNonFlowToWithFlowFlat.size());

  TH2JaggedD *jag2 = dynamic_cast<TH2JaggedD *>(jag->Clone("clone"));

  assert(jag2->fBinMappingNonFlowToWithFlowFlat.size());
  assert(jag->fBinMappingNonFlowToWithFlowFlat.size() == jag2->fBinMappingNonFlowToWithFlowFlat.size());

  TFile fw("TH2JagReadTest.root", "RECREATE");
  TH2Poly *p = jag->ToTH2Poly();
  fw.WriteTObject(p, "poly");
  TH2Poly *pc = jag2->ToTH2Poly();
  fw.WriteTObject(pc, "polyclone");
  fw.Close();
}
