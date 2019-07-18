
#define TH2JAGGED_DEBUG
#include "TH2Jagged.h"

#include "TFile.h"
#include "TRandom3.h"
#include "TH2Poly.h"

int main() {

  std::vector<Int_t> NXBins = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<double> XLowEdges = {-1, -2, -3, -4, -5, -5, -4, -3, -2, -1};
  std::vector<double> XUpEdges = {1, 2, 3, 4, 5, 5, 4, 3, 2, 1};

  TH2JaggedD jag("testjag", "some title; some x axis; some y axis; some z axis",
                 NXBins.data(), XLowEdges.data(), XUpEdges.data(), 10, -5, 5);

  NXBins.clear();
  XLowEdges.clear();
  XUpEdges.clear();
  for (size_t i = 0; i < 100; ++i) {
    NXBins.push_back(100);
    XLowEdges.push_back(-5);
    XUpEdges.push_back(5);
  }

  for (size_t i = 0; i < 1000000; ++i) {
    double x = gRandom->Gaus(0, 2);
    double y = gRandom->Gaus(0, 2);
    jag.Fill(x, y);
  }

  TFile f("TH2JagWriteTest.root", "RECREATE");
  f.WriteTObject(&jag, "jag");
  f.Close();
}
