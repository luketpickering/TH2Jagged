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
}
