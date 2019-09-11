
#include "TH2Jagged.h"

#include "TH2Poly.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>

bool operator<(JBinId const &l, JBinId const &r) {
  if (l.UniBin < r.UniBin) {
    return true;
  }
  if (l.UniBin > r.UniBin) {
    return false;
  }
  return (l.NonUniBin < r.NonUniBin);
}

template <typename ST> void TH2Jagged<ST>::BuildBinMappings() {
  Int_t GBin = 0;
  Int_t GNoFlowBin = 0;
  fBinMappingNonFlowToWithFlowFlat.clear();
  fBinMappingToNonFlowFlat.clear();
  fBinMappingToFlat.clear();

  double minU = fUniformAxis.GetBinLowEdge(1);
  double maxU = fUniformAxis.GetBinUpEdge(fUniformAxis.GetNbins());

  double minNU = std::numeric_limits<double>::max();
  double maxNU = -std::numeric_limits<double>::max();

  Int_t NNUAxes = fNonUniformAxes.size();
  for (Int_t ubin = 0; ubin < NNUAxes; ++ubin) {
    minNU = std::min(minNU, fNonUniformAxes[ubin].GetBinLowEdge(1));
    Int_t NNUBins = fNonUniformAxes[ubin].GetNbins();
    for (Int_t nubin = 0; nubin < (NNUBins + 2); ++nubin) {

      // Builds a mapping of 2D binnings to a linear binning [0,NBinsNonFlow)
      // which ignores flow bins
      if (ubin && ((ubin + 1) < NNUAxes) && nubin && (nubin < (NNUBins + 1))) {
        fBinMappingNonFlowToWithFlowFlat.push_back(GBin);
        fBinMappingToNonFlowFlat[JBinId{ubin, nubin}] = GNoFlowBin++;
      }

      fBinMappingToFlat[JBinId{ubin, nubin}] = GBin;
      GBin++;
    }
    maxNU = std::max(maxNU, fNonUniformAxes[ubin].GetBinUpEdge(NNUBins));
  }

  fMinX = GetXAxisT(minU, minNU);
  fMaxX = GetXAxisT(maxU, maxNU);
  fMinY = GetYAxisT(minU, minNU);
  fMaxY = GetYAxisT(maxU, maxNU);

  if (fBinContent.size() != size_t(GBin)) {
    fBinContent.clear();
    std::fill_n(std::back_inserter(fBinContent), GBin, 0);
  }
  if (fBinError.size() != size_t(GBin)) {
    fBinError.clear();
    std::fill_n(std::back_inserter(fBinError), GBin, 0);
  }
  if (fBinSumW2.size() != size_t(GBin)) {
    fBinSumW2.clear();
    std::fill_n(std::back_inserter(fBinSumW2), GBin, 0);
  }
}

template <typename ST>
TH2Jagged<ST>::TH2Jagged(char const *name, char const *title, Int_t NUBins,
                         Double_t UMin, Double_t UMax, Int_t const *NNUbins,
                         Double_t const *NUMin, Double_t const *NUMax,
                         bool XIsUniform) {
  TH2::SetName(name);
  TH2::SetTitle(title);
  fOTitle = title;

  fXIsUniform = XIsUniform;
  fUniformAxis = TAxis(NUBins, UMin, UMax);
  (fXIsUniform ? fXaxis : fYaxis) = fUniformAxis;

  // Underflow axis is the same as the first bin
  fNonUniformAxes.push_back(TAxis(NNUbins[0], NUMin[0], NUMax[0]));
  for (Int_t ubin = 0; ubin < NUBins; ++ubin) {
    fNonUniformAxes.push_back(TAxis(NNUbins[ubin], NUMin[ubin], NUMax[ubin]));
  }
  // Overflow axis is the same as the last bin
  fNonUniformAxes.push_back(
      TAxis(NNUbins[NUBins - 1], NUMin[NUBins - 1], NUMax[NUBins - 1]));

  BuildBinMappings();
}

template <typename ST>
TH2Jagged<ST>::TH2Jagged(char const *name, char const *title, Int_t NUBins,
                         Double_t const *UBinEdges, Int_t const *NNUbins,
                         Double_t const **NUBinEdges, bool XIsUniform) {
  TH2::SetName(name);
  TH2::SetTitle(title);
  fOTitle = title;

  fXIsUniform = XIsUniform;
  fUniformAxis = TAxis(NUBins, UBinEdges);
  (fXIsUniform ? fXaxis : fYaxis) = fUniformAxis;

  // Underflow axis is the same as the first bin
  fNonUniformAxes.push_back(TAxis(NNUbins[0], NUBinEdges[0]));
  for (Int_t ubin = 0; ubin < NUBins; ++ubin) {
    fNonUniformAxes.push_back(TAxis(NNUbins[ubin], NUBinEdges[ubin]));
  }
  // Overflow axis is the same as the last bin
  fNonUniformAxes.push_back(TAxis(NNUbins[NUBins - 1], NUBinEdges[NUBins - 1]));

  BuildBinMappings();
}

template <typename ST>
bool TH2Jagged<ST>::CheckConsistency(const TH2Jagged *h) {
  if (!TH1::CheckEqualAxes(&fUniformAxis, &h->fUniformAxis)) {
    return false;
  }

  Int_t NNUAxes = fNonUniformAxes.size();
  for (Int_t ubin = 0; ubin < NNUAxes; ++ubin) {
    if (!TH1::CheckEqualAxes(&fNonUniformAxes[ubin],
                             &h->fNonUniformAxes[ubin])) {
      return false;
    }
  }
  return true;
}

template <typename ST> TH2Jagged<ST>::TH2Jagged() {
  ResetUniformAxis();
  // Added to elide a bug where the copy of the bin mappings was missed, could
  // be removed in future versions
  BuildBinMappings();
}
template <typename ST>
TH2Jagged<ST>::TH2Jagged(char const *name, char const *title, Int_t NXbins,
                         Double_t XMin, Double_t XMax, Int_t const *NYbins,
                         Double_t const *YMin, Double_t const *YMax)
    : TH2Jagged(name, title, NXbins, XMin, XMax, NYbins, YMin, YMax, true) {}

template <typename ST>
TH2Jagged<ST>::TH2Jagged(char const *name, char const *title,
                         Int_t const *NXbins, Double_t const *XMin,
                         Double_t const *XMax, Int_t NYbins, Double_t YMin,
                         Double_t YMax)
    : TH2Jagged(name, title, NYbins, YMin, YMax, NXbins, XMin, XMax, false) {}

template <typename ST>
TH2Jagged<ST>::TH2Jagged(char const *name, char const *title, Int_t NXbins,
                         Double_t const *XBinEdges, Int_t const *NYbins,
                         Double_t const **YBinEdges)
    : TH2Jagged(name, title, NXbins, XBinEdges, NYbins, YBinEdges, true) {}

template <typename ST>
TH2Jagged<ST>::TH2Jagged(char const *name, char const *title,
                         Int_t const *NXbins, Double_t const **XBinEdges,
                         Int_t NYbins, Double_t const *YBinEdges)
    : TH2Jagged(name, title, NYbins, YBinEdges, NXbins, XBinEdges, false) {}

template <typename ST> Int_t TH2Jagged<ST>::Fill(Double_t x, Double_t y) {
  return Fill(x, y, 1);
}
template <typename ST>
Int_t TH2Jagged<ST>::Fill(Double_t x, Double_t y, Double_t w) {
  return FillKnownBin(TH2Jagged<ST>::FindFixBin(x, y), w);
}

template <typename ST>
Int_t TH2Jagged<ST>::FillKnownBin(Int_t gbin, Double_t w) {
  fBinContent[gbin] += w;
  // From here:
  // https://www.pp.rhul.ac.uk/~cowan/stat/notes/errors_with_weights.pdf
  fBinSumW2[gbin] += pow(w, 2);
  fBinError[gbin] = sqrt(fBinSumW2[gbin]);

  return fBinContent[gbin];
}

template <typename ST> Int_t TH2Jagged<ST>::GetNbins() const {
  // This allows loops that look like TH1 loops
  return (fBinContent.size() - 2);
}

template <typename ST> Int_t TH2Jagged<ST>::GetNbinsNonFlow() const {
  return fBinMappingNonFlowToWithFlowFlat.size();
}

template <typename ST>
Int_t TH2Jagged<ST>::GetGBinFromNonFlowIndex(Int_t nfidx) {
  if ((nfidx >= 0) &&
      (size_t(nfidx) < fBinMappingNonFlowToWithFlowFlat.size())) {
    return fBinMappingNonFlowToWithFlowFlat[nfidx];
  }
  return -1;
}

template <typename ST>
Int_t TH2Jagged<ST>::GetBin(Int_t binx, Int_t biny, Int_t) const {
  binx = std::max(0, binx);
  biny = std::max(0, biny);

  Int_t ubin = GetUniformAxisT(binx, biny);
  ubin = std::min(ubin, Int_t(fNonUniformAxes.size() - 1));

  Int_t nubin = GetNonUniformAxisT(binx, biny);
  nubin = std::min(nubin, fNonUniformAxes[ubin].GetNbins() + 1);

  return fBinMappingToFlat.at({ubin, nubin});
}
template <typename ST>
void TH2Jagged<ST>::GetBinXYZ(Int_t gbin, Int_t &binx, Int_t &biny,
                              Int_t &) const {
  gbin = std::max(0, gbin);
  // Default top right bin
  Int_t maxubin = fUniformAxis.GetNbins() + 1;
  Int_t maxnubin = fNonUniformAxes.back().GetNbins() + 1;
  binx = GetXAxisT(maxubin, maxnubin);
  biny = GetYAxisT(maxubin, maxnubin);

  for (std::map<JBinId, Int_t>::const_iterator b_it = fBinMappingToFlat.begin();
       b_it != fBinMappingToFlat.end(); ++b_it) {
    if (b_it->second == gbin) {
      binx = GetXAxisT(b_it->first.UniBin, b_it->first.NonUniBin);
      biny = GetYAxisT(b_it->first.UniBin, b_it->first.NonUniBin);
      break;
    }
  }
}

template <typename ST>
TAxis const *TH2Jagged<ST>::GetNonUniformAxis(Int_t gbin) const {
  Int_t x, y, z;
  GetBinXYZ(gbin, x, y, z);
  Int_t ubin = GetUniformAxisT(x, y);
  return &fNonUniformAxes[ubin];
}

template <typename ST>
Int_t TH2Jagged<ST>::FindFixBin(Double_t x, Double_t y, Double_t) const {
  Double_t u = GetUniformAxisT(x, y);
  Double_t nu = GetNonUniformAxisT(x, y);

  Int_t ubin = fUniformAxis.FindFixBin(u);
  Int_t nubin = fNonUniformAxes[ubin].FindFixBin(nu);
  return fBinMappingToFlat.at({ubin, nubin});
}
template <typename ST>
Int_t TH2Jagged<ST>::FindFixBinNoFlow(Double_t x, Double_t y, Double_t) const {
  Double_t u = GetUniformAxisT(x, y);
  Double_t nu = GetNonUniformAxisT(x, y);

  Int_t ubin = fUniformAxis.FindFixBin(u);
  Int_t nubin = fNonUniformAxes[ubin].FindFixBin(nu);
  if ((ubin == 0) || (ubin >= (fUniformAxis.GetNbins() + 1))) {
    return -1;
  }
  if ((nubin == 0) || (nubin >= (fNonUniformAxes[ubin].GetNbins() + 1))) {
    return -1;
  }
  return fBinMappingToNonFlowFlat.at({ubin, nubin});
}

template <typename ST> void TH2Jagged<ST>::Scale(Double_t c, Option_t *option) {
  std::string opt = option;
  std::transform(opt.begin(), opt.end(), opt.begin(), ::tolower);

  bool doWidth = false;
  if (opt.find("width") != std::string::npos) {
    doWidth = true;
  }

  for (Int_t ubin = 1; ubin < (fUniformAxis.GetNbins() + 1); ++ubin) {
    std::pair<double, double> ubinedges{fUniformAxis.GetBinLowEdge(ubin),
                                        fUniformAxis.GetBinUpEdge(ubin)};
    for (Int_t nubin = 1; nubin < (fNonUniformAxes[ubin].GetNbins() + 1);
         ++nubin) {
      std::pair<double, double> nubinedges{
          fNonUniformAxes[ubin].GetBinLowEdge(nubin),
          fNonUniformAxes[ubin].GetBinUpEdge(nubin)};

      Int_t gbin = fBinMappingToFlat.at({ubin, nubin});

      double bc = fBinContent[gbin];
      double be = fBinError[gbin];

      double ws = doWidth ? (ubinedges.second - ubinedges.first) *
                                (nubinedges.second - nubinedges.first)
                          : 1;

      fBinContent[gbin] = (bc * c / ws);
      fBinError[gbin] = (be * c / ws);
    }
  }
}

template <typename ST> Double_t TH2Jagged<ST>::GetBinContent(Int_t gbin) const {
  gbin = std::max(0, gbin);
  gbin = std::min(gbin, Int_t(fBinContent.size() - 1));
  return fBinContent[gbin];
}
template <typename ST>
Double_t TH2Jagged<ST>::GetBinContent(Int_t binx, Int_t biny) const {
  return GetBinContent(GetBin(binx, biny));
}

template <typename ST> void TH2Jagged<ST>::SetBinContent(Int_t gbin, ST c) {
  gbin = std::max(0, gbin);
  gbin = std::min(gbin, Int_t(fBinContent.size() - 1));
  fBinContent[gbin] = c;
}
template <typename ST>
void TH2Jagged<ST>::SetBinContent(Int_t binx, Int_t biny, ST c) {
  SetBinContent(GetBin(binx, biny), c);
}

template <typename ST> Double_t TH2Jagged<ST>::GetBinError(Int_t gbin) const {
  gbin = std::max(0, gbin);
  gbin = std::min(gbin, Int_t(fBinError.size() - 1));
  return fBinError[gbin];
}
template <typename ST>
Double_t TH2Jagged<ST>::GetBinError(Int_t binx, Int_t biny) const {
  return GetBinError(GetBin(binx, biny));
}

template <typename ST> void TH2Jagged<ST>::SetBinError(Int_t gbin, ST c) {
  gbin = std::max(0, gbin);
  gbin = std::min(gbin, Int_t(fBinError.size() - 1));
  fBinError[gbin] = c;
}
template <typename ST>
void TH2Jagged<ST>::SetBinError(Int_t binx, Int_t biny, ST c) {
  SetBinError(GetBin(binx, biny), c);
}

template <typename ST> TH2Poly *TH2Jagged<ST>::ToTH2Poly() const {
  std::stringstream ss("");
  ss << TH2::GetName() << "_poly";
  TH2Poly *pol = new TH2Poly(ss.str().c_str(), fOTitle.c_str(), fMinX, fMaxX,
                             fMinY, fMaxY);

  for (Int_t ubin = 1; ubin < fUniformAxis.GetNbins() + 1; ++ubin) {
    std::pair<double, double> ubinedges{fUniformAxis.GetBinLowEdge(ubin),
                                        fUniformAxis.GetBinUpEdge(ubin)};

    for (Int_t nubin = 1; nubin < fNonUniformAxes[ubin].GetNbins() + 1;
         ++nubin) {
      std::pair<double, double> nubinedges{
          fNonUniformAxes[ubin].GetBinLowEdge(nubin),
          fNonUniformAxes[ubin].GetBinUpEdge(nubin)};
      std::pair<double, double> XBinEdges = GetXAxisT(ubinedges, nubinedges);
      std::pair<double, double> YBinEdges = GetYAxisT(ubinedges, nubinedges);

      Int_t b = pol->AddBin(XBinEdges.first, YBinEdges.first, XBinEdges.second,
                            YBinEdges.second);

      pol->SetBinContent(b, GetBinContent(fBinMappingToFlat.at({ubin, nubin})));
      pol->SetBinError(b, GetBinError(fBinMappingToFlat.at({ubin, nubin})));
    }
  }
  return pol;
}

template <typename ST> void TH2Jagged<ST>::Reset(Option_t *) {
  std::fill_n(fBinContent.begin(), fBinContent.size(), 0);
  std::fill_n(fBinError.begin(), fBinError.size(), 0);
  std::fill_n(fBinSumW2.begin(), fBinSumW2.size(), 0);
}

template <typename T> struct sqrtfunct {
  T operator()(T const &o) { return std::sqrt(o); }
};

template <typename ST> void TH2Jagged<ST>::RecalculateErrors(Option_t *option) {
  std::string opt(option);
  std::transform(opt.begin(), opt.end(), opt.begin(), ::tolower);

  if (!opt.size() || (opt.find("wpoisson") != std::string::npos)) {
    std::transform(fBinSumW2.begin(), fBinSumW2.end(), fBinError.begin(),
                   sqrtfunct<ST>());
  } else if (opt.find("poisson") != std::string::npos) {
    std::transform(fBinContent.begin(), fBinContent.end(), fBinError.begin(),
                   sqrtfunct<ST>());
  }
}

namespace {
template <typename T> struct quad {
  double fFactL, fFactR;
  quad(double factl = 1, double factr = 1) {
    fFactL = factl;
    fFactR = factr;
  }
  T operator()(T const &l, T const &r) {
    return sqrt(fFactL * l * fFactL * l + fFactR * r * fFactR * r);
  }
};

template <typename T> struct plus {
  double fFactL, fFactR;
  plus(double factl = 1, double factr = 1) {
    fFactL = factl;
    fFactR = factr;
  }
  T operator()(T const &l, T const &r) { return (fFactL * l + fFactR * r); }
};
} // namespace

template <typename ST>
Bool_t TH2Jagged<ST>::Add(const TH2Jagged<ST> *h1, Double_t c1) {
  if (!CheckConsistency(h1)) {
    throw "[ERROR]: Inconsistent axes in TH2Jagged<ST>::Add";
  }
  std::transform(fBinContent.begin(), fBinContent.end(),
                 h1->fBinContent.begin(), fBinContent.begin(), plus<ST>(1, c1));

  std::transform(fBinError.begin(), fBinError.end(), h1->fBinError.begin(),
                 fBinError.begin(), quad<ST>(1, c1));

  std::transform(fBinSumW2.begin(), fBinSumW2.end(), h1->fBinSumW2.begin(),
                 fBinSumW2.begin(), plus<ST>(1, c1));
  return true;
}

template <typename ST>
Bool_t TH2Jagged<ST>::Add(const TH2Jagged<ST> *h1, const TH2Jagged<ST> *h2,
                          Double_t c1, Double_t c2) {
  if (!CheckConsistency(h1)) {
    throw "[ERROR]: Inconsistent axes in TH2Jagged<ST>::Add";
  }
  if (!CheckConsistency(h2)) {
    throw "[ERROR]: Inconsistent axes in TH2Jagged<ST>::Add";
  }
  std::transform(h1->fBinContent.begin(), h1->fBinContent.end(),
                 h2->fBinContent.begin(), fBinContent.begin(),
                 plus<ST>(c1, c2));

  std::transform(h1->fBinError.begin(), h1->fBinError.end(),
                 h2->fBinError.begin(), fBinError.begin(), quad<ST>(c1, c2));

  std::transform(h1->fBinSumW2.begin(), h1->fBinSumW2.end(),
                 h2->fBinSumW2.begin(), fBinSumW2.begin(), plus<ST>(c1, c2));
  return true;
}

template <typename ST> Bool_t TH2Jagged<ST>::Add(const TH1 *h1, Double_t c1) {
  if (dynamic_cast<const TH2Jagged<ST> *>(h1)) {
    return Add(dynamic_cast<const TH2Jagged<ST> *>(h1), c1);
  } else {
    throw std::string("[ERROR]: Cannot add ") + h1->ClassName() +
        " to TH2Jagged<ST>.";
  }
}

template <typename ST>
Bool_t TH2Jagged<ST>::Add(const TH1 *h1, const TH1 *h2, Double_t c1,
                          Double_t c2) {
  if (dynamic_cast<const TH2Jagged<ST> *>(h1) &&
      dynamic_cast<const TH2Jagged<ST> *>(h2)) {
    return Add(dynamic_cast<const TH2Jagged<ST> *>(h1),
               dynamic_cast<const TH2Jagged<ST> *>(h2), c1, c2);
  } else {
    throw std::string("[ERROR]: Cannot add ") + h1->ClassName() + " to a " +
        h2->ClassName() + " with a TH2Jagged<ST>";
  }
}

template <typename ST> bool TH2Jagged<ST>::IsFlowBin(Int_t gbin) const {
  Int_t x, y, z;
  GetBinXYZ(gbin, x, y, z);
  Int_t ubin = GetUniformAxisT(x, y);
  Int_t nubin = GetNonUniformAxisT(x, y);

  return ((ubin == 0) || (ubin == (fUniformAxis.GetNbins() + 1)) ||
          (nubin == 0) || (nubin == (fNonUniformAxes[ubin].GetNbins() + 1)));
}

template <typename ST>
Double_t TH2Jagged<ST>::Integral(Option_t *option) const {
  std::string opt = option;
  std::transform(opt.begin(), opt.end(), opt.begin(), ::tolower);

  bool doWidth = false;
  if (opt.find("width") != std::string::npos) {
    doWidth = true;
  }

  Double_t integ = 0;

  for (Int_t ubin = 1; ubin < (fUniformAxis.GetNbins() + 1); ++ubin) {
    std::pair<double, double> ubinedges{fUniformAxis.GetBinLowEdge(ubin),
                                        fUniformAxis.GetBinUpEdge(ubin)};
    for (Int_t nubin = 1; nubin < (fNonUniformAxes[ubin].GetNbins() + 1);
         ++nubin) {
      std::pair<double, double> nubinedges{
          fNonUniformAxes[ubin].GetBinLowEdge(nubin),
          fNonUniformAxes[ubin].GetBinUpEdge(nubin)};

      Int_t gbin = fBinMappingToFlat.at({ubin, nubin});

      double bc = fBinContent[gbin];
      double ws = doWidth ? (ubinedges.second - ubinedges.first) *
                                (nubinedges.second - nubinedges.first)
                          : 1;

      integ += (bc * ws);
    }
  }
  return integ;
}

template <typename ST>
typename TH2TypeTraits<ST>::TH2Type *
TH2Jagged<ST>::ToUniformTH2(Option_t *option) const {
  std::string opt = option;
  std::transform(opt.begin(), opt.end(), opt.begin(), ::tolower);

  // BY default should spread events between sub bins, but if it is a width
  // histogram, then the ratio of bin widths has already be accounted for
  bool doWidth = true;
  if (opt.find("width") != std::string::npos) {
    doWidth = false;
  }

  std::vector<double> UBins;
  for (int ubin = 1; ubin < (fUniformAxis.GetNbins() + 1); ++ubin) {
    UBins.push_back(fUniformAxis.GetBinLowEdge(ubin));
  }
  UBins.push_back(fUniformAxis.GetBinUpEdge(fUniformAxis.GetNbins()));

  std::vector<double> NUBins;
  NUBins.push_back(GetNonUniformAxisT(fMinX, fMinY));
  double cedge = NUBins.back();
  while (true) {
    for (int ubin = 1; ubin < (fUniformAxis.GetNbins() + 1); ++ubin) {
      for (int nubin = 1; nubin < (fNonUniformAxes[ubin].GetNbins() + 1);
           ++nubin) {
        // std::cout << "x: " << GetXAxisT(ubin, nubin)
        //           << ", y: " << GetYAxisT(ubin, nubin) << " low edge = "
        //           << fNonUniformAxes[ubin].GetBinLowEdge(nubin) << std::endl;
        if (fNonUniformAxes[ubin].GetBinLowEdge(nubin) > NUBins.back()) {

          if (cedge != NUBins.back()) {
            cedge = std::min(cedge, fNonUniformAxes[ubin].GetBinLowEdge(nubin));
          } else {
            cedge = fNonUniformAxes[ubin].GetBinLowEdge(nubin);
          }

          // std::cout << "  --NCE " << cedge << " (" << NUBins.back() << ")"
          //           << std::endl;
          break;
        }
      }
      // std::cout << "x: " << GetXAxisT(ubin, fNonUniformAxes[ubin].GetNbins())
      //           << ", y: " << GetYAxisT(ubin,
      //           fNonUniformAxes[ubin].GetNbins())
      //           << " up edge = "
      //           << fNonUniformAxes[ubin].GetBinUpEdge(
      //                  fNonUniformAxes[ubin].GetNbins())
      //           << std::endl;
      if (fNonUniformAxes[ubin].GetBinUpEdge(fNonUniformAxes[ubin].GetNbins()) >
          NUBins.back()) {

        if (cedge != NUBins.back()) {
          cedge = std::min(cedge, fNonUniformAxes[ubin].GetBinUpEdge(
                                      fNonUniformAxes[ubin].GetNbins()));
        } else {
          cedge = fNonUniformAxes[ubin].GetBinUpEdge(
              fNonUniformAxes[ubin].GetNbins());
        }

        // std::cout << "  --NCE " << cedge << " (" << NUBins.back() << ")"
        //           << std::endl;
      }
    }
    if (cedge != NUBins.back()) {
      NUBins.push_back(cedge);
      std::cout << "--NE: " << cedge << std::endl;
    } else {
      break;
    }
  }

  std::vector<double> *XBins = GetXAxisT(&UBins, &NUBins);
  std::vector<double> *YBins = GetYAxisT(&UBins, &NUBins);

  std::stringstream ss("");
  ss << TH2::GetName() << "_th2";

  typename TH2TypeTraits<ST>::TH2Type *ret = new
      typename TH2TypeTraits<ST>::TH2Type(ss.str().c_str(), fOTitle.c_str(),
                                          (XBins->size() - 1), XBins->data(),
                                          (YBins->size() - 1), YBins->data());

  for (Int_t xb_it = 1; xb_it < (ret->GetXaxis()->GetNbins() + 1); ++xb_it) {
    for (Int_t yb_it = 1; yb_it < (ret->GetYaxis()->GetNbins() + 1); ++yb_it) {
      double xc = ret->GetXaxis()->GetBinCenter(xb_it);
      double yc = ret->GetYaxis()->GetBinCenter(yb_it);

      Int_t gbin = TH2Jagged<ST>::FindFixBin(xc, yc);
      if (IsFlowBin(gbin)) {
        continue;
      }

      double ws = 1;
      if (doWidth) {
        Int_t x, y, z;
        GetBinXYZ(gbin, x, y, z);
        Int_t ubin = GetUniformAxisT(x, y);
        Int_t nubin = GetNonUniformAxisT(x, y);

        double jag_bw = (fNonUniformAxes[ubin].GetBinUpEdge(nubin) -
                         fNonUniformAxes[ubin].GetBinLowEdge(nubin)) *
                        (fUniformAxis.GetBinUpEdge(ubin) -
                         fUniformAxis.GetBinLowEdge(ubin));

        double uni_bw = (ret->GetXaxis()->GetBinUpEdge(xb_it) -
                         ret->GetXaxis()->GetBinLowEdge(xb_it)) *
                        (ret->GetYaxis()->GetBinUpEdge(yb_it) -
                         ret->GetYaxis()->GetBinLowEdge(yb_it));

        ws = uni_bw / jag_bw;
      }

      Int_t rgbin = ret->GetBin(xb_it, yb_it);
      ret->SetBinContent(rgbin, fBinContent[gbin] * ws);
      ret->SetBinError(rgbin, 0);
    }
  }
  return ret;
}

template <typename ST>
typename TH2TypeTraits<ST>::TH1Type *
TH2Jagged<ST>::NonUniformSlice(Int_t ubin) const {
  if (ubin >= fUniformAxis.GetNbins()) {
    std::cout << "[ERROR]: Requested non uniform slice " << ubin
              << " but only have " << fUniformAxis.GetNbins()
              << " bins along the uniform axis." << std::endl;
    throw;
  }
  ubin++;

  std::vector<double> NUBins;
  std::vector<std::pair<double, double>> NUBinContent;

  NUBins.push_back(fNonUniformAxes[ubin].GetBinLowEdge(1));

  for (int nubin = 1; nubin < fNonUniformAxes[ubin].GetNbins() + 1; ++nubin) {
    NUBins.push_back(NUBins.back() + fNonUniformAxes[ubin].GetBinWidth(nubin));
    Int_t gbin = fBinMappingToFlat.at({ubin, nubin});
    NUBinContent.push_back(
        std::make_pair<double, double>(GetBinContent(gbin), GetBinError(gbin)));
  }

  std::stringstream ss("");
  ss << TH2::GetName() << "_ubin" << ubin;
  T1T *t1 = new T1T(ss.str().c_str(), "", (NUBins.size() - 1), NUBins.data());
  for (size_t i = 0; i < NUBinContent.size(); ++i) {
    t1->SetBinContent(i + 1, NUBinContent[i].first);
    t1->SetBinError(i + 1, NUBinContent[i].second);
  }
  return t1;
}

template <typename ST>
typename TH2TypeTraits<ST>::TH1Type *TH2Jagged<ST>::ToFlatTH1() const {
  std::vector<double> NUBins;
  std::vector<std::pair<double, double>> NUBinContent;

  NUBins.push_back(fNonUniformAxes[0].GetBinLowEdge(1));

  for (int ubin = 1; ubin < fUniformAxis.GetNbins() + 1; ++ubin) {
    for (int nubin = 1; nubin < fNonUniformAxes[ubin].GetNbins() + 1; ++nubin) {
      NUBins.push_back(NUBins.back() +
                       fNonUniformAxes[ubin].GetBinWidth(nubin));
      Int_t gbin = fBinMappingToFlat.at({ubin, nubin});
      NUBinContent.push_back(std::make_pair<double, double>(GetBinContent(gbin),
                                                            GetBinError(gbin)));
    }
  }

  std::stringstream ss("");
  ss << TH2::GetName() << "_1D";
  T1T *t1 = new T1T(ss.str().c_str(), "", (NUBins.size() - 1), NUBins.data());
  for (size_t i = 0; i < NUBinContent.size(); ++i) {
    t1->SetBinContent(i + 1, NUBinContent[i].first);
    t1->SetBinError(i + 1, NUBinContent[i].second);
  }
  return t1;
}
template <typename ST>
void TH2Jagged<ST>::SetBinContentFromFlatTH1(T1T const *h) {
  std::vector<std::pair<double, double>> NUBinContent;

  for (Int_t i = 0; i < h->GetXaxis()->GetNbins(); ++i) {
    NUBinContent.push_back(std::make_pair<double, double>(
        h->GetBinContent(i + 1), h->GetBinError(i + 1)));
  }
  size_t ctr = 0;
  for (int ubin = 1; ubin < (fUniformAxis.GetNbins() + 1); ++ubin) {
    for (int nubin = 1; nubin < (fNonUniformAxes[ubin].GetNbins() + 1);
         ++nubin) {
      Int_t gbin = fBinMappingToFlat.at({ubin, nubin});

      SetBinContent(gbin, NUBinContent[ctr].first);
      SetBinError(gbin, NUBinContent[ctr].second);
      ctr++;
    }
  }
}

template <typename ST> void TH2Jagged<ST>::ResetUniformAxis() {
  (fXIsUniform ? fXaxis : fYaxis) = fUniformAxis;
}

template <typename ST>
TObject *TH2Jagged<ST>::Clone(char const *newname) const {
  TH2Jagged<ST> *n = new TH2Jagged<ST>();

  n->fBinMappingToFlat = fBinMappingToFlat;
  n->fBinMappingToNonFlowFlat = fBinMappingToNonFlowFlat;
  n->fBinMappingNonFlowToWithFlowFlat = fBinMappingNonFlowToWithFlowFlat;
  n->fUniformAxis = fUniformAxis;
  n->fNonUniformAxes = fNonUniformAxes;
  n->fBinContent = fBinContent;
  n->fBinError = fBinError;
  n->fBinSumW2 = fBinSumW2;
  n->fXIsUniform = fXIsUniform;
  n->fMinX = fMinX;
  n->fMaxX = fMaxX;
  n->fMinY = fMinY;
  n->fMaxY = fMaxY;
  n->fOTitle = fOTitle;
  if (newname) {
    n->SetName(newname);
  } else {
    n->SetName(GetName());
  }
  n->SetTitle(GetTitle());

  // Added to elide a bug where the copy of the bin mappings was missed, could
  // be removed in future versions
  n->BuildBinMappings();

  n->ResetUniformAxis();

  return n;
}

template <typename ST> void TH2Jagged<ST>::Draw(Option_t *option) {
  TH2Poly *p = ToTH2Poly();
  p->Draw(option);
}

TH2JaggedF *ToTHJaggedF(TH2JaggedD const *d, char const *newname) {
  TH2JaggedF *f = new TH2JaggedF();

  f->fBinMappingToFlat = d->fBinMappingToFlat;
  f->fUniformAxis = d->fUniformAxis;
  f->fNonUniformAxes = d->fNonUniformAxes;

  std::fill_n(std::back_inserter(f->fBinContent), d->GetNbins() + 2, 0);
  std::fill_n(std::back_inserter(f->fBinError), d->GetNbins() + 2, 0);
  std::fill_n(std::back_inserter(f->fBinSumW2), d->GetNbins() + 2, 0);

  for (Int_t i = 0; i < d->GetNbins() + 2; ++i) {
    f->fBinContent[i] = float(d->fBinContent[i]);
    f->fBinError[i] = float(d->fBinError[i]);
    f->fBinSumW2[i] = float(d->fBinSumW2[i]);
  }
  f->fXIsUniform = d->fXIsUniform;
  f->fMinX = d->fMinX;
  f->fMaxX = d->fMaxX;
  f->fMinY = d->fMinY;
  f->fMaxY = d->fMaxY;
  f->fOTitle = d->fOTitle;
  if (newname) {
    f->SetName(newname);
  } else {
    f->SetName(d->GetName());
  }
  f->SetTitle(d->GetTitle());

  return f;
}

template <typename ST>
TH1D *TH2Jagged<ST>::DoProjection(bool onX, const char *name, Int_t firstbin,
                                  Int_t lastbin, Option_t *) const {

  bool IsU = (onX == fXIsUniform);
  if (!IsU && (lastbin == -1)) {
    lastbin = firstbin;
  }
  if (!IsU && (firstbin != lastbin)) {
    throw "Cannot project multiple bins for non uniform axis";
  }
  if (IsU) {
    throw "Projecting to uniform axis is not yet implemented.";
  }

  Int_t ubin = firstbin - 1;
  if (ubin >= fUniformAxis.GetNbins()) {
    std::cout << "[ERROR]: Requested non uniform slice " << ubin
              << " but only have " << fUniformAxis.GetNbins()
              << " bins along the uniform axis." << std::endl;
    throw;
  }
  ubin++;

  std::vector<double> NUBins;
  std::vector<std::pair<double, double>> NUBinContent;

  NUBins.push_back(fNonUniformAxes[ubin].GetBinLowEdge(1));

  for (int nubin = 1; nubin < fNonUniformAxes[ubin].GetNbins() + 1; ++nubin) {
    NUBins.push_back(NUBins.back() + fNonUniformAxes[ubin].GetBinWidth(nubin));
    Int_t gbin = fBinMappingToFlat.at({ubin, nubin});
    NUBinContent.push_back(
        std::make_pair<double, double>(GetBinContent(gbin), GetBinError(gbin)));
  }

  std::stringstream ss("");
  if (name) {
    ss << name;
  } else {
    ss << TH2::GetName() << "_ubin" << ubin;
  }
  TH1D *t1 = new TH1D(ss.str().c_str(), "", (NUBins.size() - 1), NUBins.data());
  for (size_t i = 0; i < NUBinContent.size(); ++i) {
    t1->SetBinContent(i + 1, NUBinContent[i].first);
    t1->SetBinError(i + 1, NUBinContent[i].second);
  }
  return t1;
}
