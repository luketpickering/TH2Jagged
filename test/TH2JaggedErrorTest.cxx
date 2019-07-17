#include "TH2Jagged.h"
#include "TRandom3.h"

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

  TH2JaggedD jag_uy("test_uy", "", NXBins.data(), XLowEdges.data(),
                    XUpEdges.data(), 5, 0, 5);

  TH2JaggedD jag_ux("test_ux", "", 5, 0, 5, NXBins.data(), XLowEdges.data(),
                    XUpEdges.data());

  TH2D uni("test2", "", 5, 0, 5, 5, 0, 5);

  for (size_t i = 0; i < 1000000; ++i) {
    double x = gRandom->Gaus(0, 2);
    double y = gRandom->Gaus(0, 1);
    double w = gRandom->Uniform(0, 1);
    jag_uy.Fill(x, y, w);
    jag_ux.Fill(x, y, w);
    uni.Fill(x, y, w);
  }

  for (Int_t i = 0; i < 5; ++i) {
    for (Int_t j = 0; j < 5; ++j) {
      double jbc_uy = jag_uy.GetBinContent(i + 1, j + 1);
      double jbc_ux = jag_ux.GetBinContent(i + 1, j + 1);
      double ubc = uni.GetBinContent(i + 1, j + 1);
      double jbe_uy = jag_uy.GetBinError(i + 1, j + 1);
      double jbe_ux = jag_ux.GetBinError(i + 1, j + 1);
      double ube = uni.GetBinError(i + 1, j + 1);

      say_assert(fabs(jbc_uy-jbc_ux), <, 1E-8);
      say_assert(fabs(jbc_uy-ubc), <, 1E-8);

      say_assert(fabs(jbe_uy-jbe_ux), <, 1E-8);
      say_assert(fabs(jbe_uy-ube), <, 1E-8);

    }
  }
}
