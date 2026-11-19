#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TH1F.h"

//========================================================================
// [Exercise 3.3] Energy-deposit (dE/dx) spectrum
//
// The pulse-height spectrum of scintillator A (in volts) is converted to
// energy with a calibration curve, and its mean and most probable value
// (MPV) are calculated.
//
// Usage : ./dedx group_dedx.txt
// Input : output of step2.mkHists, one line per bin: "lowEdge[V] entries"
// Output: mean and MPV on stdout, dedx.pdf
//
// Tasks (modify the blocks marked "(k)"):
//   (1) Energy calibration: voltage [V] -> energy [MeV].
//   (2) Mean energy of the spectrum.
//
// Things to think about:
//   Q1. How do you obtain the calibration curve from your measurement?
//   Q2. How is the mean of a histogram calculated from its bin contents?
//       Which energy represents each bin?
//
// Check your answer:
//   - The MPV should be close to the expected energy loss of a minimum
//     ionizing particle in the scintillator.
//========================================================================

//======================================================================
//== (1) Energy calibration: voltage [V] -> energy [MeV]
//==     (must be monotonically increasing)
float v2eCalibration (float v) {
  return -0.31f + 10.3f * v;  // <=== modify here
}
//======================================================================

int main(int argc, char* argv[]) { // 2025.10.24 Y.Kamiya, rev. 2026.10.16

  if (argc != 2) {
    std::cout << "usage: dedx group_dedx.txt\n";
    return 1;
  }

  // === Constants (must match nBinDedx and dv in step2.mkHists) ===
  constexpr int nBins = 256;
  constexpr float dv = 0.008f; // [V]

  // === Load data ===
  std::ifstream ifs(argv[1]);
  if (!ifs) {
    std::cerr << "Error: Cannot open " << argv[1] << "\n";
    return 1;
  }

  std::vector<float> voltage(nBins);
  std::vector<float> entries(nBins);
  for (int i = 0; i < nBins; ++i) {
    ifs >> voltage[i] >> entries[i];
  }
  ifs.close();

  // === Histogram in energy ===
  // Each voltage bin [i*dv, (i+1)*dv] maps to the energy bin
  // [v2e(i*dv), v2e((i+1)*dv)], so input and output bins correspond one to one.
  std::vector<double> edges(nBins + 1);
  for (int i = 0; i <= nBins; ++i) {
    edges[i] = v2eCalibration(i * dv); // [MeV]
  }
  auto* h = new TH1F("h", "h;energy [MeV];entries", nBins, edges.data());

  // === Mean and most probable value ===
  double totalEnt = 0.0;
  double mean = 0.0;
  float mpv = 0.0f;
  float mpvEnt = 0.0f;

  for (int i = 0; i < nBins; ++i) {
    h->SetBinContent(i + 1, entries[i]);
    const float energy = v2eCalibration(voltage[i] + 0.5f * dv); // bin center [MeV]
    totalEnt += entries[i];
    //======================================================================
    //== (2) Mean energy: accumulate here, divided by totalEnt below
    mean += 1.0;  // <=== modify here
    //======================================================================

    if (mpvEnt < entries[i]) {
      mpvEnt = entries[i];
      mpv = energy;
    }
  }
  h->SetEntries(totalEnt);

  if (totalEnt <= 0.0) {
    std::cerr << "Error: no entries in " << argv[1] << "\n";
    return 1;
  }
  mean /= totalEnt;

  std::cout << "mean = " << mean << " [MeV]\n";
  std::cout << "mpv  = " << mpv << " [MeV]\n";

  // === Draw ===
  auto* c = new TCanvas("c", "Energy Distribution", 800, 600);
  h->Draw("hist");
  c->Print("dedx.pdf");

  return 0;
}
