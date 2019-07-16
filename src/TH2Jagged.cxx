
#include "TH2Jagged.h"

#include "TH2Poly.h"

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

template <class TH2T> void TH2Jagged<TH2T>::BuildBinMappings() {
  Int_t GBin = 0;

  double minU = fUniformAxis.GetBinLowEdge(1);
  double maxU = fUniformAxis.GetBinUpEdge(fUniformAxis.GetNbins());

  double minNU = std::numeric_limits<double>::max();
  double maxNU = -std::numeric_limits<double>::max();

  Int_t NNUAxes = fNonUniformAxes.size();
  for (Int_t ubin = 0; ubin < NNUAxes; ++ubin) {
    minNU = std::min(minNU, fNonUniformAxes[ubin].GetBinLowEdge(1));
    Int_t NNUBins = fNonUniformAxes[ubin].GetNbins();
    for (Int_t nubin = 0; nubin < (NNUBins + 2); ++nubin) {
      JBinId binId{ubin, nubin};
      std::cout << "GBin = " << GBin << " = {X: " << GetXAxisT(ubin, nubin)
                << ", Y: " << GetYAxisT(ubin, nubin) << "}, {U: " << ubin
                << ", NU: " << nubin << "}\n";

      fBinMappingToFlat[binId] = GBin;
      GBin++;
    }
    maxNU = std::max(maxNU, fNonUniformAxes[ubin].GetBinUpEdge(NNUBins));
  }

  fMinX = GetXAxisT(minU, minNU);
  fMaxX = GetXAxisT(maxU, maxNU);
  fMinY = GetYAxisT(minU, minNU);
  fMaxY = GetYAxisT(maxU, maxNU);

  std::fill_n(std::back_inserter(fBinContent), GBin, 0);
  std::fill_n(std::back_inserter(fBinError), GBin, 0);
  std::fill_n(std::back_inserter(fBinSumW2), GBin, 0);
}

template <class TH2T>
TH2Jagged<TH2T>::TH2Jagged(const char *name, const char *title, Int_t NUbins,
                           Double_t UMin, Double_t UMax, Int_t *NNUbins,
                           Double_t *NUMin, Double_t *NUMax, bool XIsUniform) {
  TH2::SetName(name);
  TH2::SetTitle(title);
  fOTitle = title;

  fXIsUniform = XIsUniform;
  fUniformAxis = TAxis(NUbins, UMin, UMax);

  // Underflow axis is the same as the first bin
  fNonUniformAxes.push_back(TAxis(NNUbins[0], NUMin[0], NUMax[0]));
  for (Int_t ubin = 0; ubin < NUbins; ++ubin) {
    fNonUniformAxes.push_back(TAxis(NNUbins[ubin], NUMin[ubin], NUMax[ubin]));
  }
  // Overflow axis is the same as the last bin
  fNonUniformAxes.push_back(
      TAxis(NNUbins[NUbins - 1], NUMin[NUbins - 1], NUMax[NUbins - 1]));

  BuildBinMappings();
}

template <class TH2T>
TH2Jagged<TH2T>::TH2Jagged(const char *name, const char *title, Int_t NUbins,
                           Double_t *UBinEdges, Int_t *NNUbins,
                           Double_t **NUBinEdges, bool XIsUniform) {
  TH2::SetName(name);
  TH2::SetTitle(title);
  fOTitle = title;

  fXIsUniform = XIsUniform;
  fUniformAxis = TAxis(NUbins, UBinEdges);

  // Underflow axis is the same as the first bin
  fNonUniformAxes.push_back(TAxis(NNUbins[0], NUBinEdges[0]));
  for (Int_t ubin = 0; ubin < NUbins; ++ubin) {
    fNonUniformAxes.push_back(TAxis(NNUbins[NUbins], NUBinEdges[NUbins]));
  }
  // Overflow axis is the same as the last bin
  fNonUniformAxes.push_back(TAxis(NNUbins[NUbins - 1], NUBinEdges[NUbins - 1]));

  BuildBinMappings();
}

template <class TH2T>
bool TH2Jagged<TH2T>::CheckConsistency(const TH2Jagged *h) {
  if (!TH1::CheckEqualAxes(&fUniformAxis, &h->fUniformAxis)) {
    return false;
  }

  Int_t NNUAxes = fNonUniformAxes.size();
  for (size_t ubin = 0; ubin < NNUAxes; ++ubin) {
    if (!TH1::CheckEqualAxes(&fNonUniformAxes[ubin],
                             &h->fNonUniformAxes[ubin])) {
      return false;
    }
  }
  return true;
}

template <class TH2T> TH2Jagged<TH2T>::TH2Jagged() {}
template <class TH2T>
TH2Jagged<TH2T>::TH2Jagged(const char *name, const char *title, Int_t NXbins,
                           Double_t XMin, Double_t XMax, Int_t *NYbins,
                           Double_t *YMin, Double_t *YMax)
    : TH2Jagged(name, title, NXbins, XMin, XMax, NYbins, YMin, YMax, true) {}

template <class TH2T>
TH2Jagged<TH2T>::TH2Jagged(const char *name, const char *title, Int_t *NXbins,
                           Double_t *XMin, Double_t *XMax, Int_t NYbins,
                           Double_t YMin, Double_t YMax)
    : TH2Jagged(name, title, NYbins, YMin, YMax, NXbins, XMin, XMax, false) {}

template <class TH2T>
TH2Jagged<TH2T>::TH2Jagged(const char *name, const char *title, Int_t NXbins,
                           Double_t *XBinEdges, Int_t *NYbins,
                           Double_t **YBinEdges)
    : TH2Jagged(name, title, NXbins, XBinEdges, NYbins, YBinEdges, true) {}

template <class TH2T>
TH2Jagged<TH2T>::TH2Jagged(const char *name, const char *title, Int_t *NXbins,
                           Double_t **XBinEdges, Int_t NYbins,
                           Double_t *YBinEdges)
    : TH2Jagged(name, title, NYbins, YBinEdges, NXbins, XBinEdges, false) {}

template <class TH2T> Int_t TH2Jagged<TH2T>::Fill(Double_t x, Double_t y) {
  return Fill(x, y, 1);
}
template <class TH2T>
Int_t TH2Jagged<TH2T>::Fill(Double_t x, Double_t y, Double_t w) {
  Int_t gbin = FindFixBin(x, y);
  fBinContent[gbin] += w;
  fBinSumW2[gbin] += pow(w, 2);
  fBinError[gbin] = sqrt(fBinSumW2[gbin]);

  return fBinContent[gbin];
}
template <class TH2T>
Int_t TH2Jagged<TH2T>::GetBin(Int_t binx, Int_t biny, Int_t) const {
  binx = std::max(0, binx);
  biny = std::max(0, biny);

  Int_t ubin = GetUniformAxisT(binx, biny);
  ubin = std::min(ubin, Int_t(fNonUniformAxes.size() - 1));

  Int_t nubin = GetNonUniformAxisT(binx, biny);
  nubin = std::min(nubin, fNonUniformAxes[ubin].GetNbins() + 1);

  return fBinMappingToFlat.at({ubin, nubin});
}
template <class TH2T>
void TH2Jagged<TH2T>::GetBinXYZ(Int_t gbin, Int_t &binx, Int_t &biny,
                                Int_t &binz) const {
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

template <class TH2T>
TAxis const *TH2Jagged<TH2T>::GetNonUniformAxis(Int_t ubin) const {
  ubin = std::max(0, ubin);
  ubin = std::min(ubin, Int_t(fNonUniformAxes.size() - 1));
  return &fNonUniformAxes[ubin];
}

template <class TH2T>
Int_t TH2Jagged<TH2T>::FindFixBin(Double_t x, Double_t y) const {
  Double_t u = GetUniformAxisT(x, y);
  Double_t nu = GetNonUniformAxisT(x, y);

  Int_t ubin = fUniformAxis.FindFixBin(u);
  Int_t nubin = fNonUniformAxes[ubin].FindFixBin(nu);
  return fBinMappingToFlat.at({ubin, nubin});
}

template <class TH2T>
void TH2Jagged<TH2T>::Scale(Double_t c, Option_t *option) {
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

template <class TH2T>
Double_t TH2Jagged<TH2T>::GetBinContent(Int_t gbin) const {
  gbin = std::max(0, gbin);
  gbin = std::min(gbin, Int_t(fBinContent.size() - 1));
  return fBinContent[gbin];
}
template <class TH2T>
Double_t TH2Jagged<TH2T>::GetBinContent(Int_t binx, Int_t biny) const {
  return GetBinContent(GetBin(binx, biny));
}

template <class TH2T> void TH2Jagged<TH2T>::SetBinContent(Int_t gbin, ST c) {
  gbin = std::max(0, gbin);
  gbin = std::min(gbin, Int_t(fBinContent.size() - 1));
  fBinContent[gbin] = c;
}
template <class TH2T>
void TH2Jagged<TH2T>::SetBinContent(Int_t binx, Int_t biny, ST c) {
  SetBinContent(GetBin(binx, biny), c);
}

template <class TH2T> Double_t TH2Jagged<TH2T>::GetBinError(Int_t gbin) const {
  gbin = std::max(0, gbin);
  gbin = std::min(gbin, Int_t(fBinError.size() - 1));
  return fBinError[gbin];
}
template <class TH2T>
Double_t TH2Jagged<TH2T>::GetBinError(Int_t binx, Int_t biny) const {
  return GetBinError(GetBin(binx, biny));
}

template <class TH2T> void TH2Jagged<TH2T>::SetBinError(Int_t gbin, ST c) {
  gbin = std::max(0, gbin);
  gbin = std::min(gbin, Int_t(fBinError.size() - 1));
  fBinError[gbin] = c;
}
template <class TH2T>
void TH2Jagged<TH2T>::SetBinError(Int_t binx, Int_t biny, ST c) {
  SetBinError(GetBin(binx, biny), c);
}

template <class TH2T> TH2Poly *TH2Jagged<TH2T>::ToTH2Poly() const {
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

template <class TH2T> void TH2Jagged<TH2T>::Reset(Option_t *option) {
  std::fill_n(fBinContent.begin(), 0, fBinContent.size());
  std::fill_n(fBinError.begin(), 0, fBinError.size());
  std::fill_n(fBinSumW2.begin(), 0, fBinSumW2.size());
}

template <typename T> struct quad {
  T operator()(T const &l, T const &r) { return sqrt(l * l + r * r); }
};

template <class TH2T>
Bool_t TH2Jagged<TH2T>::Add(const TH2Jagged<TH2T> *h1, Double_t c1) {
  if (!CheckConsistency(h1)) {
    throw "[ERROR]: Inconsistent axes in TH2Jagged<TH2T>::Add";
  }
  std::transform(fBinContent.begin(), fBinContent.end(),
                 h1->fBinContent.begin(), fBinContent.begin(),
                 std::plus<TH2Jagged<TH2T>::ST>());

  std::transform(fBinError.begin(), fBinError.end(), h1->fBinError.begin(),
                 fBinError.begin(), quad<TH2Jagged<TH2T>::ST>());

  std::transform(fBinSumW2.begin(), fBinSumW2.end(), h1->fBinSumW2.begin(),
                 fBinSumW2.begin(), std::plus<TH2Jagged<TH2T>::ST>());
  return true;
}

template <class TH2T>
typename TH2TypeTraits<TH2T>::TH1Type *TH2Jagged<TH2T>::ToFlatTH1() const {
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
template <class TH2T>
void TH2Jagged<TH2T>::SetBinContentFromFlatTH1(T1T const *h) {
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

template <class TH2T> TObject *TH2Jagged<TH2T>::Clone(const char *newname) {
  TH2Jagged<TH2T> *n = new TH2Jagged<TH2T>();

  n->fBinMappingToFlat = fBinMappingToFlat;
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

  return n;
}

template <class TH2T> void TH2Jagged<TH2T>::Draw(Option_t *option) {
  TH2Poly *p = ToTH2Poly();
  p->Draw(option);
}
