#ifndef TH2JAGGED_HXX_SEEN
#define TH2JAGGED_HXX_SEEN

#include "TH1.h"
#include "TH2.h"

#include <map>
#include <string>
#include <vector>

class TH2Poly;

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
bool operator<(JBinId const &l, JBinId const &r);

template <class TH2T> class TH2Jagged : public TH2 {

  using T1T = typename TH2TypeTraits<TH2T>::TH1Type;
  using ST = typename TH2TypeTraits<TH2T>::StorageType;

  std::map<JBinId, Int_t> fBinMappingToFlat;

  TAxis fUniformAxis;
  std::vector<TAxis> fNonUniformAxes;
  std::vector<ST> fBinContent;
  std::vector<ST> fBinError;
  std::vector<ST> fBinSumW2;

  bool fXIsUniform;
  double fMinX, fMaxX;
  double fMinY, fMaxY;
  std::string fOTitle;

  // Builds mappings between 'flat' and X/Y bin numbers
  void BuildBinMappings();

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
            bool XIsUniform);

  TH2Jagged(const char *name, const char *title, Int_t NUbins,
            Double_t *UBinEdges, Int_t *NNUbins, Double_t **NUBinEdges,
            bool XIsUniform);

  bool CheckConsistency(const TH2Jagged *h);

public:
  TH2Jagged();
  TH2Jagged(const char *name, const char *title, Int_t NXbins, Double_t XMin,
            Double_t XMax, Int_t *NYbins, Double_t *YMin, Double_t *YMax);

  TH2Jagged(const char *name, const char *title, Int_t *NXbins, Double_t *XMin,
            Double_t *XMax, Int_t NYbins, Double_t YMin, Double_t YMax);

  TH2Jagged(const char *name, const char *title, Int_t NXbins,
            Double_t *XBinEdges, Int_t *NYbins, Double_t **YBinEdges);

  TH2Jagged(const char *name, const char *title, Int_t *NXbins,
            Double_t **XBinEdges, Int_t NYbins, Double_t *YBinEdges);

  // Shuts up compiler warning but pulls methods that will do nothing into scope
  using TH2::Fill;
  Int_t Fill(Double_t x, Double_t y);
  Int_t Fill(Double_t x, Double_t y, Double_t w);

  Int_t GetBin(Int_t binx, Int_t biny, Int_t binz = 0) const;
  void GetBinXYZ(Int_t gbin, Int_t &binx, Int_t &biny, Int_t &binz) const;

  TAxis const *GetNonUniformAxis(Int_t gbin) const;

  Int_t FindFixBin(Double_t x, Double_t y) const;

  void Scale(Double_t c = 1, Option_t *option = "");

  // Shuts up compiler warning but pulls methods that will do nothing into scope
  using TH2::GetBinContent;
  Double_t GetBinContent(Int_t gbin) const;
  Double_t GetBinContent(Int_t binx, Int_t biny) const;

  // Shuts up compiler warning but pulls methods that will do nothing into scope
  using TH2::SetBinContent;
  void SetBinContent(Int_t gbin, ST c);
  void SetBinContent(Int_t binx, Int_t biny, ST c);

  // Shuts up compiler warning but pulls methods that will do nothing into scope
  using TH2::GetBinError;
  Double_t GetBinError(Int_t gbin) const;
  Double_t GetBinError(Int_t binx, Int_t biny) const;

  // Shuts up compiler warning but pulls methods that will do nothing into scope
  using TH2::SetBinError;
  void SetBinError(Int_t gbin, ST c);
  void SetBinError(Int_t binx, Int_t biny, ST c);

  TH2Poly *ToTH2Poly() const;

  void Reset(Option_t *option = "");

  // Shuts up compiler warning but pulls methods that will do nothing into scope
  using TH2::Add;
  Bool_t Add(const TH2Jagged<TH2T> *h1, Double_t c1 = 1);

  TObject *Clone(const char *newname = 0);
  void Draw(Option_t *option = "");

  T1T *ToFlatTH1() const;
  void SetBinContentFromFlatTH1(T1T const *h);

  virtual ~TH2Jagged() {}

  ClassDef(TH2Jagged, 1);
};

template class TH2Jagged<TH2D>;
template class TH2Jagged<TH2F>;

using TH2JaggedD = TH2Jagged<TH2D>;
using TH2JaggedF = TH2Jagged<TH2F>;

#endif
