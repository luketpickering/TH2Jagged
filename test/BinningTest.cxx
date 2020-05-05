#include "TH2Jagged.h"

#include <cassert>
#include <iostream>

#define say_assert(l, op, r)                                                   \
  {                                                                            \
    std::cout << "[ASSERT]: " #l << "(" << l << ") " << #op << " " << #r       \
              << "(" << r << ")\n";                                            \
    assert(l op r);                                                            \
  }

int main() {
  std::vector<Int_t> NXBins(5, 5);
  std::vector<double> XLowEdges(5, 0);
  std::vector<double> XUpEdges(5, 5);

  TH2JaggedD jag("test", "", NXBins.data(), XLowEdges.data(), XUpEdges.data(),
                 5, 0, 5);

  int gbin_11 = jag.GetBin(1, 1);
  say_assert(gbin_11, ==, 8); // 7 x underflow + 1 y underflow

  int gbin_12 = jag.GetBin(1, 2); // 7 x underflow + 7 x_1 + 1 y underflow
  say_assert(gbin_12, ==, 15);

  int gbin_21 = jag.GetBin(2, 1); // 7 x underflow + 1 y underflow + 1 y_1
  say_assert(gbin_21, ==, 9);

  int gbin_OF = jag.GetBin(100, 100); // 7 * 7 - 1
  say_assert(gbin_OF, ==, 48);

  int gbin_neg = jag.GetBin(-100, -100); // 7 * 7 - 1
  say_assert(gbin_neg, ==, 0);

  int gbin_find_00 = jag.FindFixBin(0, 0);
  say_assert(gbin_find_00, ==, 8);

  int gbin_find_12 = jag.FindFixBin(0.5, 1.5);
  say_assert(gbin_find_12, ==, 15);

  int gbin_find_21 = jag.FindFixBin(1.5, 0.5);
  say_assert(gbin_find_21, ==, 9);

  int gbin_find_OF = jag.FindFixBin(10, 10);
  say_assert(gbin_find_OF, ==, 48);

  int gbin_find_OFx = jag.FindFixBin(10, 0.5);
  say_assert(gbin_find_OFx, ==, 13);

  int gbin_find_UF = jag.FindFixBin(-1, -1);
  say_assert(gbin_find_UF, ==, 0);

  int gbin_find_UFy = jag.FindFixBin(4.5, -1);
  say_assert(gbin_find_UFy, ==, 5);

  std::vector<Int_t> NXBins_ylong(20, 5);
  std::vector<double> XLowEdges_ylong(20, 0);
  std::vector<double> XUpEdges_ylong(20, 5);

  for (int i = 0; i < 20; ++i) {
    NXBins_ylong[i] = i / 2 + 2;
    XLowEdges_ylong[i] = i;
    XUpEdges_ylong[i] = i + 1 + (20 * i) / 2;
  }

  TH2JaggedD jag_ylong("test", "", NXBins_ylong.data(), XLowEdges_ylong.data(),
                       XUpEdges_ylong.data(), 20, 0, 20);
  std::cout << jag_ylong.BinningString() << std::endl;
  for (int ybin_it = 0; ybin_it < jag_ylong.GetYaxis()->GetNbins(); ++ybin_it) {
    std::cout << "ybin_it: " << ybin_it << " ["
              << jag_ylong.GetYaxis()->GetBinLowEdge(ybin_it + 1) << "-"
              << jag_ylong.GetYaxis()->GetBinUpEdge(ybin_it + 1) << "]"
              << std::endl;
    TAxis const *xax = jag_ylong.GetNonUniformAxis_UniformAxisBin(ybin_it + 1);
    std::cout << "nbinsx: " << xax->GetNbins() << std::endl;

    for (int xbin_it = 0; xbin_it < xax->GetNbins(); ++xbin_it) {
      std::cout << "\txbin_it: " << xbin_it << " ["
                << xax->GetBinLowEdge(xbin_it + 1) << "-"
                << xax->GetBinUpEdge(xbin_it + 1) << "]" << std::endl;

      double ybc = jag_ylong.GetYaxis()->GetBinCenter(ybin_it + 1);
      double xbc = xax->GetBinCenter(xbin_it + 1);

      int gb = jag_ylong.GetBin(xbin_it + 1, ybin_it + 1);

      std::cout << "\t--[xbc]: " << xbc << std::endl;
      std::cout << "\t--[ybc]: " << ybc << std::endl;
      int ffb = jag_ylong.FindFixBin(xbc, ybc);

      std::cout << "\t--[gbin]: " << gb << std::endl;
      std::cout << "\t--[FFB]: " << ffb << std::endl;

      say_assert(gb, ==, ffb);
    }
  }

  std::vector<Int_t> NXBins_xlong(20, 100);
  std::vector<double> XLowEdges_xlong(20, 0);
  std::vector<double> XUpEdges_xlong(20, 5);
  for (int i = 0; i < 20; ++i) {
    NXBins_xlong[i] = ((i / 2) + 10);
    XLowEdges_xlong[i] = i;
    XUpEdges_xlong[i] = i + 1 + (20 * i) / 2;
  }

  TH2JaggedD jag_xlong("test", "", NXBins_xlong.data(), XLowEdges_xlong.data(),
                       XUpEdges_xlong.data(), 20, 0, 20);
  std::cout << jag_xlong.BinningString() << std::endl;
  for (int ybin_it = 0; ybin_it < jag_xlong.GetYaxis()->GetNbins(); ++ybin_it) {
    std::cout << "ybin_it: " << ybin_it << " ["
              << jag_xlong.GetYaxis()->GetBinLowEdge(ybin_it + 1) << "-"
              << jag_xlong.GetYaxis()->GetBinUpEdge(ybin_it + 1) << "]"
              << std::endl;
    TAxis const *xax = jag_xlong.GetNonUniformAxis_UniformAxisBin(ybin_it + 1);
    std::cout << "nbinsx: " << xax->GetNbins() << std::endl;

    for (int xbin_it = 0; xbin_it < xax->GetNbins(); ++xbin_it) {
      std::cout << "\txbin_it: " << xbin_it << " ["
                << xax->GetBinLowEdge(xbin_it + 1) << "-"
                << xax->GetBinUpEdge(xbin_it + 1) << "]" << std::endl;

      double ybc = jag_xlong.GetYaxis()->GetBinCenter(ybin_it + 1);
      double xbc = xax->GetBinCenter(xbin_it + 1);

      int gb = jag_xlong.GetBin(xbin_it + 1, ybin_it + 1);

      std::cout << "\t--[xbc]: " << xbc << std::endl;
      std::cout << "\t--[ybc]: " << ybc << std::endl;
      int ffb = jag_xlong.FindFixBin(xbc, ybc);

      std::cout << "\t--[gbin]: " << gb << std::endl;
      std::cout << "\t--[FFB]: " << ffb << std::endl;
      say_assert(gb, ==, ffb);
    }
  }
}
