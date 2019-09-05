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

  int gbin_find_00 = jag.FindFixBin(0, 0);
  say_assert(gbin_find_00, ==, 8);

  int nfbin_find_00 = jag.FindFixBinNoFlow(0, 0);
  say_assert(nfbin_find_00, ==, 0);

  int gbin_find_12 = jag.FindFixBin(0.5, 1.5);
  say_assert(gbin_find_12, ==, 15);

  int nfbin_find_12 = jag.FindFixBinNoFlow(0.5, 1.5);
  say_assert(nfbin_find_12, ==, 5);

  int gbin_find_21 = jag.FindFixBin(1.5, 0.5);
  say_assert(gbin_find_21, ==, 9);

  int nfbin_find_21 = jag.FindFixBinNoFlow(1.5, 0.5);
  say_assert(nfbin_find_21, ==, 1);

  int gbin_find_OF = jag.FindFixBin(10, 10);
  say_assert(gbin_find_OF, ==, 48);

  int nfbin_find_OF = jag.FindFixBinNoFlow(10, 10);
  say_assert(nfbin_find_OF, ==, -1);

  int gbin_find_OFx = jag.FindFixBin(10, 0.5);
  say_assert(gbin_find_OFx, ==, 13);

  int nfbin_find_OFx = jag.FindFixBinNoFlow(10, 0.5);
  say_assert(nfbin_find_OFx, ==, -1);

  int gbin_find_UF = jag.FindFixBin(-1, -1);
  say_assert(gbin_find_UF, ==, 0);

  int nfbin_find_UF = jag.FindFixBinNoFlow(-1, -1);
  say_assert(nfbin_find_UF, ==, -1);

  int gbin_find_UFy = jag.FindFixBin(4.5, -1);
  say_assert(gbin_find_UFy, ==, 5);

  int nfbin_find_UFy = jag.FindFixBinNoFlow(4.5, -1);
  say_assert(nfbin_find_UFy, ==, -1);

  int nfbin_find_top = jag.FindFixBinNoFlow(4.5, 4.5);
  say_assert(nfbin_find_top, ==, 24);
}
