
#define TH2JAGGED_DEBUG
#include "TH2Jagged.h"

#include "TFile.h"
#include "TRandom3.h"

int main() {

  std::vector<Int_t> NXBins = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<double> XLowEdges = {-1, -2, -3, -4, -5, -5, -4, -3, -2, -1};
  std::vector<double> XUpEdges = {1, 2, 3, 4, 5, 5, 4, 3, 2, 1};

  TH2D uni("uniform", "", 100, -5, 5, 100, -5, 5);

  TH2DJagged jag("testjag", "some title; some x axis; some y axis; some z axis",
                 NXBins.data(), XLowEdges.data(), XUpEdges.data(), 10, -5, 5);
  TH2DJagged jag2("testjag2",
                  "some title; some x axis; some y axis; some z axis",
                  NXBins.data(), XLowEdges.data(), XUpEdges.data(), 10, -5, 5);

  NXBins.clear();
  XLowEdges.clear();
  XUpEdges.clear();
  for (size_t i = 0; i < 100; ++i) {
    NXBins.push_back(100);
    XLowEdges.push_back(-5);
    XUpEdges.push_back(5);
  }

  TH2DJagged jaguni(
      "testjaguni", "some title; some x axis; some y axis; some z axis",
      NXBins.data(), XLowEdges.data(), XUpEdges.data(), 100, -5, 5);

  for (size_t i = 0; i < 1000000; ++i) {
    double x = gRandom->Gaus(0, 2);
    double y = gRandom->Gaus(0, 2);
    jag.Fill(x, y);
    jaguni.Fill(x, y);
    uni.Fill(x, y);
  }

  TFile f("TH2JagTestPoly.root", "RECREATE");
  TH2Poly *p = jag.ToTH2Poly();
  f.WriteTObject(p, "poly");
  jag.Scale(0.1);
  TH2Poly *psd = jag.ToTH2Poly();
  f.WriteTObject(psd, "polyscaledown");
  jag.Scale(1, "width");
  TH2Poly *pws = jag.ToTH2Poly();
  f.WriteTObject(pws, "polywidthscale");
  f.WriteTObject(&uni, "uni");
  TH2Poly *pu = jaguni.ToTH2Poly();
  f.WriteTObject(pu, "jaguni");

  jag2.Add(&jag);
  TH2Poly *p2 = jag2.ToTH2Poly();
  f.WriteTObject(p2, "poly2");
  jag2.Add(&jag);
  TH2Poly *p4 = jag2.ToTH2Poly();
  f.WriteTObject(p4, "poly2poly");

  TH1D *t1 = jag2.ToFlatTH1();
  f.WriteTObject(t1, "t1");
  jag.SetBinContentFromFlatTH1(t1);
  TH2Poly *p4reset = jag.ToTH2Poly();
  f.WriteTObject(p4reset, "p4reset");

  bool didthrow = false;
  try {
    jag2.Add(&jaguni);
  } catch (...) {
    didthrow = true;
    std::cout
        << "Cannot add TH2Jagged with differing binnings!\nTest successful."
        << std::endl;
  }

  assert(didthrow);

  f.Close();
}