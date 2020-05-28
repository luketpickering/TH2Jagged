#define TH2JAGGED_DEBUG
#include "TH2Jagged.h"

#include "TFile.h"
#include "TH2Poly.h"
#include "TRandom3.h"

#include <iostream>

int main() {

  std::vector<Int_t> NXBins = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::vector<double> XLowEdges = {-1, -2, -3, -4, -5, -4, -3, -2, -1};
  std::vector<double> XUpEdges = {1, 2, 3, 4, 5, 4, 3, 2, 1};

  TH2D uni("uniform", "", 100, -5, 5, 100, -2.25, 2.25);

  TH2JaggedD jag("testjag", "some title; some x axis; some y axis; some z axis",
                 NXBins.data(), XLowEdges.data(), XUpEdges.data(), 9, -2.25,
                 2.25);
  TH2JaggedD jag2(
      "testjag2", "some title; some x axis; some y axis; some z axis",
      NXBins.data(), XLowEdges.data(), XUpEdges.data(), 9, -2.25, 2.25);

  NXBins.clear();
  XLowEdges.clear();
  XUpEdges.clear();
  for (size_t i = 0; i < 100; ++i) {
    NXBins.push_back(100);
    XLowEdges.push_back(-5);
    XUpEdges.push_back(5);
  }

  TH2JaggedD jaguni(
      "testjaguni", "some title; some x axis; some y axis; some z axis",
      NXBins.data(), XLowEdges.data(), XUpEdges.data(), 100, -2.25, 2.25);

  for (size_t i = 0; i < 1000000; ++i) {
    double x = gRandom->Gaus(0, 2);
    double y = gRandom->Gaus(0, 1);
    jag.Fill(x, y);
    jaguni.Fill(x, y);
    uni.Fill(x, y);
  }

  auto jag_slice_whole = jag.UniformRange("jag_slice_whole", 0, 11);
  auto jag_slice_keep_low = jag.UniformRange("jag_slice_keep_low", 1, 9);
  auto jag_slice_keep_up = jag.UniformRange("jag_slice_keep_up", 3, 10);
  auto jag_slice = jag.UniformRange("jag_slice", 3, 9);

  TFile f("TH2JagTestPoly.root", "RECREATE");

  f.WriteTObject(jag_slice_whole, jag_slice_whole->GetName());
  f.WriteTObject(jag_slice_keep_low, jag_slice_keep_low->GetName());
  f.WriteTObject(jag_slice_keep_up, jag_slice_keep_up->GetName());
  f.WriteTObject(jag_slice, jag_slice->GetName());

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
  TH2D *th2jag_nw = jag.ToUniformTH2();
  f.WriteTObject(th2jag_nw, "th2jag_nw");
  TH2D *th2jag_w = jag.ToUniformTH2("width");
  f.WriteTObject(th2jag_w, "th2jag_w");

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
  TH2JaggedF *jagf = ToTHJaggedF(&jag);
  f.WriteTObject(jagf, "jagf");

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

  NXBins = {3, 3, 3, 5, 5, 5, 10, 10, 10};
  XLowEdges = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  XUpEdges = {1, 1, 1, 1, 1, 1, 1, 1, 1};
  TH2JaggedD jag_projslice(
      "jag_projslice", "some title; some x axis; some y axis; some z axis",
      NXBins.data(), XLowEdges.data(), XUpEdges.data(), 9, -2.25, 2.25);

  for (size_t i = 0; i < 1000000; ++i) {
    double x = gRandom->Gaus(0, 1);
    double y = gRandom->Gaus(0, 2);
    jag_projslice.Fill(x, y);
  }

  f.WriteTObject(&jag_projslice, "projslice");

  TH2JaggedD *projslice_low = jag_projslice.UniformRange("projslice_low", 1, 4);
  TH2JaggedD *projslice_mid = jag_projslice.UniformRange("projslice_mid", 4, 7);
  TH2JaggedD *projslice_up = jag_projslice.UniformRange("projslice_up", 7, 10);

  f.WriteTObject(projslice_low, "projslice_low");
  f.WriteTObject(projslice_mid, "projslice_mid");
  f.WriteTObject(projslice_up, "projslice_up");

  TH2JaggedD *projslice_low_merge =
      jag_projslice.UniformRange("projslice_low_merge", 1, 4, true);
  TH2JaggedD *projslice_mid_merge =
      jag_projslice.UniformRange("projslice_mid_merge", 4, 7, true);
  TH2JaggedD *projslice_up_merge =
      jag_projslice.UniformRange("projslice_up_merge", 7, 10, true);

  f.WriteTObject(projslice_low_merge, "projslice_low_merge");
  f.WriteTObject(projslice_mid_merge, "projslice_mid_merge");
  f.WriteTObject(projslice_up_merge, "projslice_up_merge");

  TH2JaggedD *fail =
      jag_projslice.UniformRange("projslice_fail_merge", 3, 8, true);
  assert(!fail);

  f.Close();
}
