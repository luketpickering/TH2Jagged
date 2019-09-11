#include "TH2Jagged_impl.h"

using TH2JaggedD = TH2Jagged<Double_t>;
using TH2JaggedF = TH2Jagged<Float_t>;

TH2JaggedF *ToTHJaggedF(TH2JaggedD const *d, char const *newname = 0);
