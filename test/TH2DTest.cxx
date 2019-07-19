
#define TH2JAGGED_DEBUG
#include "TH2Jagged.h"

#include "TFile.h"
#include "TH2Poly.h"
#include "TRandom3.h"

#include <iostream>

int main() {

  std::vector<Int_t> NXBins = {1, 2, 3, 4, 5};
  std::vector<double> XLowEdges(5, 0);
  std::vector<double> XUpEdges(5, 5);

  TH2JaggedD jag("testjag", "some title; some x axis; some y axis; some z axis",
                 NXBins.data(), XLowEdges.data(), XUpEdges.data(), 5, 0, 5);

  for (Int_t i = 0; i < jag.GetNbins(); ++i) {
    if (jag.IsFlowBin(i)) {
      continue;
    }
    jag.SetBinContent(i, i);
  }

  TFile f("TH2DTest.root", "RECREATE");
  TH2Poly *p = jag.ToTH2Poly();
  f.WriteTObject(p, "poly");
  TH2D *d = jag.ToUniformTH2("width");
  f.WriteTObject(d, "th2d");
  f.WriteTObject(&jag, "jag");

  f.Close();
}
