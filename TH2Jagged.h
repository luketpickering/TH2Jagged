#ifndef TH2JAGGED_HXX_SEEN
#define TH2JAGGED_HXX_SEEN

#include "TH1.h"
#include "TH2.h"
#include "TH2Poly.h"

#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#ifdef TH2JAGGED_DEBUG
#include <iostream>
#endif

template <typename THT> struct TH2TypeTraits {};

template <> struct TH2TypeTraits<TH2D> {
  using TH1Type = TH1D;
  using TH2Type = TH2D;
  using StorageType = Double_t;
};

template <> struct TH2TypeTraits<TH2F> {
  using TH1Type = TH1F;
  using TH2Type = TH2F;
  using StorageType = float;
};

struct JBinId {
  Int_t UniBin;
  Int_t NonUniBin;
};

bool operator<(JBinId const &l, JBinId const &r) {
  if (l.UniBin < r.UniBin) {
    return true;
  }
  if (l.UniBin > r.UniBin) {
    return false;
  }
  return (l.NonUniBin < r.NonUniBin);
}

template <typename TH2T> class TH2Jagged : public TH2TypeTraits<TH2T>::TH2Type {

  using T1T = typename TH2TypeTraits<TH2T>::TH1Type;
  using T2T = typename TH2TypeTraits<TH2T>::TH2Type;
  using ST = typename TH2TypeTraits<TH2T>::StorageType;

  std::vector<T1T> fHistos;

  std::map<JBinId, Int_t> fBinMappingToFlat;
  std::map<Int_t, JBinId> fBinMappingFromFlat;
  Int_t fNAllBins;
  TAxis fUniformAxis;
  bool fXIsUniform;
  double fMinX, fMaxX;
  double fMinY, fMaxY;
  std::string fOTitle;

  void BuildBinMappings() {
    Int_t GBin = 0;

    double minU = fUniformAxis.GetBinLowEdge(1);
    double maxU = fUniformAxis.GetBinUpEdge(fUniformAxis.GetNbins());

    double minNU = std::numeric_limits<double>::max();
    double maxNU = -std::numeric_limits<double>::max();

    Int_t NHistos = fHistos.size();
    for (Int_t NUIdx = 0; NUIdx < NHistos; ++NUIdx) {
      minNU = std::min(minNU, fHistos[NUIdx].GetXaxis()->GetBinLowEdge(1));
      Int_t NBins = fHistos[NUIdx].GetXaxis()->GetNbins();
      for (Int_t bi_it = 0; bi_it < (NBins + 2); ++bi_it) {
        JBinId binId{NUIdx + 1, bi_it};
        fBinMappingToFlat[binId] = GBin;
        fBinMappingFromFlat[GBin] = binId;
        GBin++;
      }
      maxNU = std::max(maxNU, fHistos[NUIdx].GetXaxis()->GetBinUpEdge(NBins));
    }

    fMinX = GetXAxisT(minU, minNU);
    fMaxX = GetXAxisT(maxU, maxNU);
    fMinY = GetYAxisT(minU, minNU);
    fMaxY = GetYAxisT(maxU, maxNU);

    fNAllBins = GBin;
  }

  template <typename T> T GetUniformAxisT(T x, T y) const {
    return fXIsUniform ? x : y;
  }

  template <typename T> T GetNonUniformAxisT(T x, T y) const {
    return fXIsUniform ? y : x;
  }

  template <typename T> T GetXAxisT(T u, T nu) const {
    return fXIsUniform ? u : nu;
  }

  template <typename T> T GetYAxisT(T u, T nu) const {
    return fXIsUniform ? nu : u;
  }

  TH2Jagged(const char *name, const char *title, Int_t NUbins, Double_t UMin,
            Double_t UMax, Int_t *NNUbins, Double_t *NUMin, Double_t *NUMax,
            bool XIsUniform) {
    T2T::SetName(name);
    T2T::SetTitle(title);
    fOTitle = title;

    fXIsUniform = XIsUniform;
    fUniformAxis = TAxis(NUbins, UMin, UMax);

    std::stringstream name_ss("");
    name_ss << name << "_UF";
    fHistos.push_back(
        T1T(name_ss.str().c_str(), title, NNUbins[0], NUMin[0], NUMax[0]));
    for (Int_t ub = 0; ub < NUbins; ++ub) {
      name_ss.str("");

      name_ss << name << "_" << ub;
      fHistos.push_back(
          T1T(name_ss.str().c_str(), title, NNUbins[ub], NUMin[ub], NUMax[ub]));
    }
    name_ss.str("");
    name_ss << name << "_OF";
    fHistos.push_back(T1T(name_ss.str().c_str(), title, NNUbins[NUbins - 1],
                          NUMin[NUbins - 1], NUMax[NUbins - 1]));

    BuildBinMappings();
  }

  TH2Jagged(const char *name, const char *title, Int_t NUbins,
            Double_t *UBinEdges, Int_t *NNUbins, Double_t **NUBinEdges,
            bool XIsUniform) {
    T2T::SetName(name);
    T2T::SetTitle(title);
    fOTitle = title;

    fXIsUniform = XIsUniform;
    fUniformAxis = TAxis(NUbins, UBinEdges);

    std::stringstream name_ss("");
    name_ss << name << "_UF";
    fHistos.push_back(
        T1T(name_ss.str().c_str(), title, NNUbins[0], NUBinEdges[0]));
    for (Int_t ub = 0; ub < NUbins; ++ub) {
      name_ss.str("");
      name_ss << name << "_" << ub;
      fHistos.push_back(
          T1T(name_ss.str().c_str(), title, NNUbins[ub], NUBinEdges[ub]));
    }
    name_ss.str("");
    name_ss << name << "_OF";
    fHistos.push_back(T1T(name_ss.str().c_str(), title, NNUbins[NUbins - 1],
                          NUBinEdges[NUbins - 1]));

    BuildBinMappings();
  }

  bool CheckConsistency(const TH2Jagged<TH2T> *h) {
    if (!TH1::CheckEqualAxes(&fUniformAxis, &h->fUniformAxis)) {
      return false;
    }
    for (size_t ubin = 0; ubin < fHistos.size(); ++ubin) {

      if (!TH1::CheckEqualAxes(fHistos[ubin].GetXaxis(),
                               h->fHistos[ubin].GetXaxis())) {
        return false;
      }
    }
    return true;
  }

public:
  TH2Jagged(const char *name, const char *title, Int_t NXbins, Double_t XMin,
            Double_t XMax, Int_t *NYbins, Double_t *YMin, Double_t *YMax)
      : TH2Jagged(name, title, NXbins, XMin, XMax, NYbins, YMin, YMax, true) {}

  TH2Jagged(const char *name, const char *title, Int_t *NXbins, Double_t *XMin,
            Double_t *XMax, Int_t NYbins, Double_t YMin, Double_t YMax)
      : TH2Jagged(name, title, NYbins, YMin, YMax, NXbins, XMin, XMax, false) {}

  TH2Jagged(const char *name, const char *title, Int_t NXbins,
            Double_t *XBinEdges, Int_t *NYbins, Double_t **YBinEdges)
      : TH2Jagged(name, title, NXbins, XBinEdges, NYbins, YBinEdges, true) {}

  TH2Jagged(const char *name, const char *title, Int_t *NXbins,
            Double_t **XBinEdges, Int_t NYbins, Double_t *YBinEdges)
      : TH2Jagged(name, title, NYbins, YBinEdges, NXbins, XBinEdges, false) {}

  Int_t Fill(Double_t x, Double_t y) { return Fill(x, y, 1); }
  Int_t Fill(Double_t x, Double_t y, Double_t w) {
    Double_t u = GetUniformAxisT(x, y);
    Double_t nu = GetNonUniformAxisT(x, y);
    Int_t ubin = fUniformAxis.FindFixBin(u);
    fHistos[ubin].Fill(nu, w);
  }
  Int_t GetBin(Int_t binx, Int_t biny, Int_t) const {
    Int_t binu = GetUniformAxisT(binx, biny);

    if (binu >= (fHistos.size() + 2)) {
      return (fNAllBins - 1);
    }

    Int_t binnu = GetNonUniformAxisT(binx, biny);
    if (binnu >= (fHistos[binu].GetXaxis()->GetNbins() + 2)) {
      return (fNAllBins - 1);
    }

    return fBinMappingToFlat.at({binu, binnu});
  }

  void Scale(Double_t c = 1, Option_t *option = "") {
    std::string opt = option;
    std::transform(opt.begin(), opt.end(), opt.begin(), ::tolower);

    bool doWidth = false;
    if (opt.find("width") != std::string::npos) {
      doWidth = true;
    }

    for (Int_t ubin = 0; ubin < fUniformAxis.GetNbins(); ++ubin) {
      std::pair<double, double> ubinedges{fUniformAxis.GetBinLowEdge(ubin + 1),
                                          fUniformAxis.GetBinUpEdge(ubin + 1)};
      for (Int_t nubin = 0; nubin < fHistos[ubin + 1].GetXaxis()->GetNbins();
           ++nubin) {
        std::pair<double, double> nubinedges{
            fHistos[ubin + 1].GetXaxis()->GetBinLowEdge(nubin + 1),
            fHistos[ubin + 1].GetXaxis()->GetBinUpEdge(nubin + 1)};

        double bc = fHistos[ubin + 1].GetBinContent(nubin + 1);
        double be = fHistos[ubin + 1].GetBinError(nubin + 1);

        double ws = doWidth ? (ubinedges.second - ubinedges.first) *
                                  (nubinedges.second - nubinedges.first)
                            : 1;

        fHistos[ubin + 1].SetBinContent(nubin + 1, bc * c / ws);
        fHistos[ubin + 1].SetBinError(nubin + 1, be * c / ws);
      }
    }
  }

  Double_t GetBinContent(Int_t bin) const {
    JBinId const &jbin = fBinMappingFromFlat.at(bin);
    Int_t binx = GetXAxisT(jbin.UniBin, jbin.NonUniBin);
    Int_t biny = GetYAxisT(jbin.UniBin, jbin.NonUniBin);
    return GetBinContent(binx, biny);
  }
  Double_t GetBinContent(Int_t binx, Int_t biny) const {
    Int_t binu = GetUniformAxisT(binx, biny);
    binu = std::min(binu, Int_t(fHistos.size() + 2));

    Int_t binnu = GetNonUniformAxisT(binx, biny);
    binnu = std::min(binnu, Int_t(fHistos[binu].GetXaxis()->GetNbins() + 2));

    return fHistos[binu].GetBinContent(binnu);
  }

  void SetBinContent(Int_t bin, ST c) {
    JBinId const &jbin = fBinMappingFromFlat.at(bin);
    Int_t binx = GetXAxisT(jbin.UniBin, jbin.NonUniBin);
    Int_t biny = GetYAxisT(jbin.UniBin, jbin.NonUniBin);
    return SetBinContent(binx, biny, c);
  }
  void SetBinContent(Int_t binx, Int_t biny, ST c) {
    Int_t binu = GetUniformAxisT(binx, biny);
    if (binu >= (fHistos.size() + 2)) {
      return;
    }

    Int_t binnu = GetNonUniformAxisT(binx, biny);
    if (binnu >= (fHistos[binu].GetXaxis()->GetNbins() + 2)) {
      return;
    }

    return fHistos[binu].SetBinContent(binnu, c);
  }

  Double_t GetBinError(Int_t bin) const {
    JBinId const &jbin = fBinMappingFromFlat.at(bin);
    Int_t binx = GetXAxisT(jbin.UniBin, jbin.NonUniBin);
    Int_t biny = GetYAxisT(jbin.UniBin, jbin.NonUniBin);
    return GetBinError(binx, biny);
  }
  Double_t GetBinError(Int_t binx, Int_t biny) const {
    Int_t binu = GetUniformAxisT(binx, biny);
    binu = std::min(binu, Int_t(fHistos.size() + 2));

    Int_t binnu = GetNonUniformAxisT(binx, biny);
    binnu = std::min(binnu, Int_t(fHistos[binu].GetXaxis()->GetNbins() + 2));

    return fHistos[binu].GetBinError(binnu);
  }

  void SetBinError(Int_t bin, ST c) {
    JBinId const &jbin = fBinMappingFromFlat.at(bin);
    Int_t binx = GetXAxisT(jbin.UniBin, jbin.NonUniBin);
    Int_t biny = GetYAxisT(jbin.UniBin, jbin.NonUniBin);
    return SetBinError(binx, biny, c);
  }
  void SetBinError(Int_t binx, Int_t biny, ST c) {
    Int_t binu = GetUniformAxisT(binx, biny);
    if (binu >= (fHistos.size() + 2)) {
      return;
    }

    Int_t binnu = GetNonUniformAxisT(binx, biny);
    if (binnu >= (fHistos[binu].GetXaxis()->GetNbins() + 2)) {
      return;
    }

    return fHistos[binu].SetBinError(binnu, c);
  }

  TH2Poly *ToTH2Poly() const {
    std::stringstream ss("");
    ss << T2T::GetName() << "_poly";
    TH2Poly *pol = new TH2Poly(ss.str().c_str(), fOTitle.c_str(), fMinX, fMaxX,
                               fMinY, fMaxY);

    for (size_t ubin = 0; ubin < fUniformAxis.GetNbins(); ++ubin) {
      std::pair<double, double> ubinedges{fUniformAxis.GetBinLowEdge(ubin + 1),
                                          fUniformAxis.GetBinUpEdge(ubin + 1)};

      for (Int_t nubin = 0; nubin < fHistos[ubin + 1].GetXaxis()->GetNbins();
           ++nubin) {
        std::pair<double, double> nubinedges{
            fHistos[ubin + 1].GetXaxis()->GetBinLowEdge(nubin + 1),
            fHistos[ubin + 1].GetXaxis()->GetBinUpEdge(nubin + 1)};
        std::pair<double, double> XBinEdges = GetXAxisT(ubinedges, nubinedges);
        std::pair<double, double> YBinEdges = GetYAxisT(ubinedges, nubinedges);

        Int_t b = pol->AddBin(XBinEdges.first, YBinEdges.first,
                              XBinEdges.second, YBinEdges.second);
        pol->SetBinContent(b, fHistos[ubin + 1].GetBinContent(nubin + 1));
        pol->SetBinError(b, fHistos[ubin + 1].GetBinError(nubin + 1));
      }
    }
    return pol;
  }
  void Reset(Option_t *option = "") {
    for (size_t ubin = 0; ubin < fHistos.size(); ++ubin) {
      fHistos[ubin].Reset(option);
    }
  }
  void SetDirectory(TDirectory *dir) {
    for (size_t ubin = 0; ubin < fHistos.size(); ++ubin) {
      fHistos[ubin].SetDirectory(dir);
    }
  }

  Bool_t Add(const TH2Jagged<TH2T> *h1, Double_t c1 = 1) {
    if (!CheckConsistency(h1)) {
      throw "Inconsistent axes in TH2Jagged::Add";
    }
    for (Int_t ubin = 0; ubin < fUniformAxis.GetNbins(); ++ubin) {
      fHistos[ubin + 1].Add(&h1->fHistos.at(ubin + 1), c1);
    }
    return false;
  }

  T1T *ToFlatTH1() const {
    std::vector<double> NUBins;
    std::vector<std::pair<double, double>> NUBinContent;
    NUBins.push_back(fHistos[1].GetXaxis()->GetBinLowEdge(1));

    for (int ub = 0; ub < fUniformAxis.GetNbins(); ++ub) {
      for (int nub = 0; nub < fHistos[ub + 1].GetXaxis()->GetNbins(); ++nub) {
        NUBins.push_back(NUBins.back() +
                         fHistos[ub + 1].GetXaxis()->GetBinWidth(nub + 1));
        NUBinContent.push_back(std::make_pair<double, double>(
            fHistos[ub + 1].GetBinContent(nub + 1),
            fHistos[ub + 1].GetBinError(nub + 1)));
      }
    }

    std::stringstream ss("");
    ss << T2T::GetName() << "_1D";
    T1T *t1 = new T1T(ss.str().c_str(), "", (NUBins.size() - 1), NUBins.data());
    for (size_t i = 0; i < NUBinContent.size(); ++i) {
      t1->SetBinContent(i + 1, NUBinContent[i].first);
      t1->SetBinError(i + 1, NUBinContent[i].second);
    }
    return t1;
  }
  void SetBinContentFromFlatTH1(T1T const *h) {
    std::vector<std::pair<double, double>> NUBinContent;

    for (size_t i = 0; i < h->GetXaxis()->GetNbins(); ++i) {
      NUBinContent.push_back(std::make_pair<double, double>(
          h->GetBinContent(i + 1), h->GetBinError(i + 1)));
    }
    size_t ctr = 0;
    for (int ub = 0; ub < fUniformAxis.GetNbins(); ++ub) {
      for (int nub = 0; nub < fHistos[ub + 1].GetXaxis()->GetNbins(); ++nub) {
        fHistos[ub + 1].SetBinContent(nub + 1, NUBinContent[ctr].first);
        fHistos[ub + 1].SetBinError(nub + 1, NUBinContent[ctr].second);
        ctr++;
      }
    }
  }
};

using TH2DJagged = TH2Jagged<TH2D>;
using TH2FJagged = TH2Jagged<TH2F>;

#endif
