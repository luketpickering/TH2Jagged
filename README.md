# TH2Jagged

A TH2 subclass that allows the binning along one axis to be a function of
the other axis. The solution to this is often a TH2Poly, which for rectangular
binning is overkill.

The fill and find bin should be equivalent to standard TH2s. N.B. you cannot
use the Get[X,Y]Axis() methods, instead GetUniformAxis() and
GetNonUniformAxis(int gbin) are provided. A FindFixBin(x,y) is also provided
for convenience, as such code requiring access to the bin information cannot
use the TH2 interface, but TH2::Fill, TH2::[Get,Set]Bin[Content,Error], and
TH2::Scale should work.

Reading and writing to root files is available, and if libTH2Jagged.so is
loaded via `.L` in a cint/cling session you will even be able to draw them to a
TCanvas (via a transient TH2Poly used for visualization).
