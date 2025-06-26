#include "LUTMaker.h"

#include <chrono>

#include "Misc_Functions.h"

float LUTMaker::LUTValues[LUTMaker::SNSTEPS_TRANS][LUTMaker::SNSTEPS_RC]
                         [LUTMaker::SNSTEPS_DRIFT][LUTMaker::SNSTEPS_D]
                         [LUTMaker::SNSTEPS_PHI];
=
// Constructor
LUTMaker::LUTMaker()
{
   v_phi = linspace(0, 90, SNSTEPS_PHI);
   v_impact = linspace(0, PAD_DIAG / 2, SNSTEPS_D);
   v_drift = linspace(0, 1000, SNSTEPS_DRIFT);
   v_RC = linspace(RCmin, RCmax, SNSTEPS_RC);
   v_Dt = {DtwithB / pow(10, 3.5), DtwithoutB / pow(10, 3.5)};

   fp_trackmodel = new TrackModel();

   c1 = new TCanvas("c1", "c1", 1800, 1350);
   ptstyle = SetMyStyle();
   gROOT->SetStyle(ptstyle->GetName());
   gStyle->SetPadLeftMargin(0.04);
   gStyle->SetPadRightMargin(0.11);
   gStyle->SetPadBottomMargin(0.06);
   gStyle->SetPadTopMargin(0.08);
   gPad->UseCurrentStyle();
   gStyle->SetOptStat(0);
}

// Destructor
LUTMaker::~LUTMaker()
{
   delete fp_trackmodel;
   delete ptstyle;
   delete c1;
}

// Display all settings of the LUT maker
void LUTMaker::PrintSettings()
{
   std::cout << "MAKING LUT WITH SETTINGS:" << std::endl;
   std::cout << "Peaking time: " << fp_trackmodel->GetPeakingTime() << std::endl;
   std::cout << "ws:           " << fp_trackmodel->Getws() << std::endl;
   std::cout << "expfactor:    " << fp_trackmodel->Getexpfactor() << std::endl;
   std::cout << "arg:          " << fp_trackmodel->Getarg() << std::endl;
   std::cout << std::endl;

   std::cout << "Number of steps in phi:    " << SNSTEPS_PHI << std::endl;
   std::cout << "Number of steps in impact: " << SNSTEPS_D << std::endl;
   std::cout << "Number of steps in drift:  " << SNSTEPS_DRIFT << std::endl;
   std::cout << "Number of steps in RC:     " << SNSTEPS_RC << std::endl;
   std::cout << "Number of steps in Dt:     " << SNSTEPS_TRANS << std::endl;
   std::cout << "Step size in phi:          " << sSTEP_PHI << std::endl;
   std::cout << "Step size in impact:       " << sSTEP_IMPACT << std::endl;
   std::cout << "Step size in drift:        " << sSTEP_DRIFT << std::endl;
   std::cout << "Step size in RC:           " << sSTEP_RC << std::endl;
   std::cout << "Step size in Dt:           " << stepSizeTrans << std::endl;

   std::cout << "Dt values:                 ";
   for (const auto &val : v_Dt) {
      std::cout << (int)(val * pow(10, 3.5)) << " ";
   }
   std::cout << std::endl;
   std::cout << "RC values:                 ";
   for (const auto &val : v_RC) {
      std::cout << val << " ";
   }
   std::cout << std::endl;
}

// Compute length in the pad for a given (phi, d) pair
double LUTMaker::ComputeLength(const double &phi_rad, const double &d)
{
   std::vector<float> v_x, v_y;

   // Compute intersections with lines defining the pad borders
   float y_xmin = GetY(0, phi_rad, d);
   float y_xmax = GetY(SXWIDTH, phi_rad, d);
   float x_ymin = GetX(0, phi_rad, d);
   float x_ymax = GetX(SYHEIGHT, phi_rad, d);
   if (0 <= y_xmin and y_xmin < SYHEIGHT) {
      v_x.push_back(0);
      v_y.push_back(y_xmin);
   }
   if (0 <= y_xmax and y_xmax < SYHEIGHT) {
      v_x.push_back(SXWIDTH);
      v_y.push_back(y_xmax);
   }
   if (0 <= x_ymin and x_ymin < SXWIDTH) {
      v_x.push_back(x_ymin);
      v_y.push_back(0);
   }
   if (0 <= x_ymax and x_ymax < SXWIDTH) {
      v_x.push_back(x_ymax);
      v_y.push_back(SYHEIGHT);
   }
   if (v_x.size() == 2 and v_y.size() == 2) {
      return std::sqrt(std::pow(v_x[1] - v_x[0], 2) + std::pow(v_y[1] - v_y[0], 2));
   } else
      return 0;
}

// Compute the length in the pad for all (phi, d) pairs
void LUTMaker::ComputeLengthMap()
{
   for (int i = 0; i < SNSTEPS_PHI; i++)
      for (int j = 0; j < SNSTEPS_D; j++)
         arr_length[i][j] = ComputeLength(v_phi[i] / 180 * M_PI, v_impact[j]);
}

// Compute scale factor
void LUTMaker::MakeLUT()
{
   int peakingTime = (int)fp_trackmodel->GetPeakingTime();
   TFile tfileLUT(Form("Output_LUT/dEdx_XP_LUT_Dt%d_%d_RC_%d_%d_nDrift_%d_PT%d.root",
                       (int)(v_Dt[0] * pow(10, 3.5)), (int)(v_Dt[1] * pow(10, 3.5)),
                       RCmin, RCmax, SNSTEPS_DRIFT, peakingTime),
                  "RECREATE");
   TTree ttreeLUT("treeLUTdEdx", "treeLUTdEdx");
   double d, phi;
   float length, RC, transvDiff, drift, scalefactor;
   ttreeLUT.Branch("impact", &d);
   ttreeLUT.Branch("angle", &phi);
   ttreeLUT.Branch("RC", &RC);
   ttreeLUT.Branch("transvDiff", &transvDiff);
   ttreeLUT.Branch("drift", &drift);
   ttreeLUT.Branch("scalefactor", &scalefactor);

   auto start = std::chrono::high_resolution_clock::now();
   double time = 0;

   std::cout << "Starting LUT looping..." << std::endl;
   for (int itransvDiff = 0; itransvDiff < SNSTEPS_TRANS; itransvDiff++) {
      transvDiff = v_Dt[itransvDiff];
      std::cout << " -- transvDiff: " << transvDiff << "(" << transvDiff * pow(10, 3.5)
                << ")" << std::endl;
      auto start_transvDiff = std::chrono::high_resolution_clock::now();
      for (int iRC = 0; iRC < SNSTEPS_RC; iRC++) {
         std::cout << " ---- RC: " << v_RC[iRC] << std::endl;
         auto start_RC = std::chrono::high_resolution_clock::now();
         RC = v_RC[iRC];
         for (int idrift = 0; idrift < SNSTEPS_DRIFT; idrift++) {
            std::cout << " ------ drift: " << v_drift[idrift];
            auto start_drift = std::chrono::high_resolution_clock::now();
            drift = v_drift[idrift];
            for (int iphi = 0; iphi < SNSTEPS_PHI; iphi++) {
               phi = v_phi[iphi];
               for (int id = 0; id < SNSTEPS_D; id++) {
                  d = v_impact[id];
                  length = arr_length[iphi][id];
                  fp_trackmodel->ComputeAmplitudeLoss(length, d, phi, RC, drift,
                                                      transvDiff);
                  scalefactor = 1 / fp_trackmodel->GetAmplitudeLoss();
                  if (std::isnan(scalefactor) || std::isinf(scalefactor) ||
                      scalefactor <= 0)
                     continue;
                  ttreeLUT.Fill();
               }
            }
            auto end_drift = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_drift = end_drift - start_drift;
            time += duration_drift.count();
            std::cout << " => " << duration_drift.count() << " s" << std::endl;
         }
         auto end_RC = std::chrono::high_resolution_clock::now();
         std::chrono::duration<double> duration_RC = end_RC - start_RC;
         time += duration_RC.count();
         std::cout << " ---- RC " << v_RC[iRC] << " => " << duration_RC.count() << " s"
                   << std::endl;
      }
      auto end_transvDiff = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration_transvDiff =
         end_transvDiff - start_transvDiff;
      time += duration_transvDiff.count();
      std::cout << " -- transvDiff " << transvDiff << " => "
                << duration_transvDiff.count() << " s" << std::endl;
   }

   tfileLUT.Write();
   tfileLUT.Close();

   auto end = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> duration = end - start;
   time += duration.count();
   std::cout << "LUT done in: " << duration.count() << " s" << std::endl;
}

// Load LUT
void LUTMaker::LoadLUT(std::string LUTpath)
{
   // Initialidrifte the LUT with zeros
   for (int itransvDiff = 0; itransvDiff < SNSTEPS_TRANS; itransvDiff++)
      for (int iRC = 0; iRC < SNSTEPS_RC; iRC++)
         for (int idrift = 0; idrift < SNSTEPS_DRIFT; idrift++)
            for (int id = 0; id < SNSTEPS_D; id++)
               for (int iphi = 0; iphi < SNSTEPS_PHI; iphi++)
                  LUTValues[itransvDiff][iRC][idrift][id][iphi] = 0;

   // Load the LUT from the file
   TFile tfile(LUTpath.c_str(), "READ");
   TTree &ptree = *(TTree *)tfile.Get("treeLUTdEdx");
   double d, phi;
   float transvDiff, RC, drift, scalefactor;
   ptree.SetBranchAddress("transvDiff", &transvDiff);
   ptree.SetBranchAddress("RC", &RC);
   ptree.SetBranchAddress("drift", &drift);
   ptree.SetBranchAddress("impact", &d);
   ptree.SetBranchAddress("angle", &phi);
   ptree.SetBranchAddress("scalefactor", &scalefactor);
   // TH2F h2_scalefactor("h2_scalefactor", "Scale factor", SNSTEPS_PHI, 0, 90+1e-6,
   // SNSTEPS_D, 0, PAD_DIAG/2+1e-6); std::ofstream debug_log("debug1.log");
   for (int i = 0; i < ptree.GetEntries(); i++) {
      ptree.GetEntry(i);
      int itransvDiff =
         (int)std::round((transvDiff * pow(10, 3.5) - v_Dt[0]) / stepSizeTrans);
      int iRC = (int)std::round((RC - RCmin) / sSTEP_RC);
      int idrift = (int)std::round(drift / sSTEP_DRIFT);
      int iphi = (int)std::round(phi / sSTEP_PHI);
      int id = (int)std::round(d / sSTEP_IMPACT);
      // std::cout << transvDiff << " " << RC << " " << drift << " " << d << " " << phi <<
      // " " << scalefactor << std::endl; std::cout << itransvDiff << " " << iRC << " " <<
      // idrift << " " << id << " " << iphi << std::endl;
      LUTValues[itransvDiff][iRC][idrift][id][iphi] = scalefactor;
      // h2_scalefactor.SetBinContent(h2_scalefactor.FindBin(phi, d), scalefactor);
      // if(scalefactor!=0) debug_log << phi << " " << d << " " << scalefactor <<
      // std::endl;
   }
   // debug_log.close();
   // h2_scalefactor.Draw("COLZ");
   // c1->SaveAs("Output_PDF/scalefactor_map1.pdf");
   tfile.Close();
}

// Get factor from LUT
float LUTMaker::GetFactorFromLUT(const double &transvDiff, const double &RC,
                                 const double &drift, const double &d, const double &phi)
{ // keep double
   int itransvDiff = (int)(transvDiff * pow(10, 3.5) - v_Dt[0]) / stepSizeTrans;
   int iRC = (int)(RC - RCmin) / sSTEP_RC;
   float idrift = drift / sSTEP_DRIFT;
   float idrift_min =
      std::min(std::floor(drift / sSTEP_DRIFT), (double)SNSTEPS_DRIFT - 1);
   float idrift_max = std::max(std::ceil(drift / sSTEP_DRIFT), 0.);
   float iphi = phi / sSTEP_PHI;
   float iphi_min = std::min(std::floor(phi / sSTEP_PHI), (double)SNSTEPS_PHI - 1);
   float iphi_max = std::max(std::ceil(phi / sSTEP_PHI), 0.);
   float id = d / sSTEP_IMPACT;
   float id_min = std::min(std::floor(d / sSTEP_IMPACT), (double)SNSTEPS_D - 1);
   float id_max = std::max(std::ceil(d / sSTEP_IMPACT), 0.);

   // weights
   double w_drift, w_d, w_phi;
   if (idrift_min == idrift_max)
      w_drift = 1;
   else
      w_drift = 1 - (idrift - idrift_min) / (idrift_max - idrift_min);
   if (id_min == id_max)
      w_d = 1;
   else
      w_d = 1 - (id - id_min) / (id_max - id_min);
   if (iphi_min == iphi_max)
      w_phi = 1;
   else
      w_phi = 1 - (iphi - iphi_min) / (iphi_max - iphi_min);

   // std::cout << "transvDiff =    " << transvDiff << " | itransvDiff = " << itransvDiff
   // << std::endl; std::cout << "RC =    " << RC << " | iRC = " << iRC << std::endl;
   // std::cout << "drift = " << drift << " | idrift_min = " << idrift_min << " | idrift =
   // " << idrift
   //           << " | idrift_max = " << idrift_max << " | w_drift = " << w_drift <<
   //           std::endl;
   // std::cout << "d =     " << d << " | id_min = " << id_min << " | id = " << id << " |
   // id_max = " << id_max
   //           << " | w_d = " << w_d << std::endl;
   // std::cout << "phi =   " << phi << " | iphi_min = " << iphi_min << " | iphi = " <<
   // iphi
   //           << " | iphi_max = " << iphi_max << " | w_phi = " << w_phi << std::endl;

   // Interpolation
   float factor = 0;
   factor += w_drift * w_d * w_phi *
             LUTValues[itransvDiff][iRC][(int)idrift_min][(int)id_min][(int)iphi_min];
   factor += w_drift * w_d * (1 - w_phi) *
             LUTValues[itransvDiff][iRC][(int)idrift_min][(int)id_min][(int)iphi_max];
   factor += w_drift * (1 - w_d) * w_phi *
             LUTValues[itransvDiff][iRC][(int)idrift_min][(int)id_max][(int)iphi_min];
   factor += w_drift * (1 - w_d) * (1 - w_phi) *
             LUTValues[itransvDiff][iRC][(int)idrift_min][(int)id_max][(int)iphi_max];
   factor += (1 - w_drift) * w_d * w_phi *
             LUTValues[itransvDiff][iRC][(int)idrift_max][(int)id_min][(int)iphi_min];
   factor += (1 - w_drift) * w_d * (1 - w_phi) *
             LUTValues[itransvDiff][iRC][(int)idrift_max][(int)id_min][(int)iphi_max];
   factor += (1 - w_drift) * (1 - w_d) * w_phi *
             LUTValues[itransvDiff][iRC][(int)idrift_max][(int)id_max][(int)iphi_min];
   factor += (1 - w_drift) * (1 - w_d) * (1 - w_phi) *
             LUTValues[itransvDiff][iRC][(int)idrift_max][(int)id_max][(int)iphi_max];

   // std::cout << std::endl;
   // std::cout << LUTValues[itransvDiff][iRC][(int)idrift_min + 0][(int)id_min +
   // 1][(int)iphi_min + 1] << "  "
   //           << LUTValues[itransvDiff][iRC][(int)idrift_min + 1][(int)id_min +
   //           1][(int)iphi_min + 1] << std::endl;
   // std::cout << LUTValues[itransvDiff][iRC][(int)idrift_min + 0][(int)id_min +
   // 0][(int)iphi_min + 1] << "  "
   //           << LUTValues[itransvDiff][iRC][(int)idrift_min + 1][(int)id_min +
   //           0][(int)iphi_min + 1] << std::endl;
   // std::cout << std::endl;
   // std::cout << LUTValues[itransvDiff][iRC][(int)idrift_min + 0][(int)id_min +
   // 1][(int)iphi_min + 0] << "  "
   //           << LUTValues[itransvDiff][iRC][(int)idrift_min + 1][(int)id_min +
   //           1][(int)iphi_min + 0] << std::endl;
   // std::cout << LUTValues[itransvDiff][iRC][(int)idrift_min + 0][(int)id_min +
   // 0][(int)iphi_min + 0] << "  "
   //           << LUTValues[itransvDiff][iRC][(int)idrift_min + 1][(int)id_min +
   //           0][(int)iphi_min + 0] << std::endl;
   // std::cout << std::endl;

   return factor;
}

// DRAWING FUNCTIONS
// -----------------------------------------------------------------------------------------------

// Draw length map
void LUTMaker::DrawLengthMap()
{
   TH2F h2("h2", "Length map", v_phi.size(), 0, 90, v_impact.size(), 0, PAD_DIAG / 2);
   for (int i = 0; i < (int)v_phi.size(); ++i) {
      for (int j = 0; j < (int)v_impact.size(); ++j) {
         h2.SetBinContent(i + 1, j + 1, arr_length[i][j]);
      }
   }
   h2.Draw("COLZ");
   c1->SaveAs("Output_PDF/length_map.pdf");
}

// Draw LUT
void LUTMaker::DrawLUT(const double &RC, const double &drift, const double &transvDiff)
{
   TH2F h2_scalefactor("h2_scalefactor", "Scale factor", SNSTEPS_PHI, 0, 90 + 1e-6,
                       SNSTEPS_D, 0, PAD_DIAG / 2 + 1e-6);
   std::ofstream debug_log("debug2.log");
   for (int i = 0; i < SNSTEPS_PHI; i++) {
      for (int j = 0; j < SNSTEPS_D; j++) {
         double phi = v_phi[i];  // keep double
         double d = v_impact[j]; // keep double
         float scalefactor = GetFactorFromLUT(transvDiff, RC, drift, d, phi);
         h2_scalefactor.SetBinContent(h2_scalefactor.FindBin(phi, d), scalefactor);
         if (scalefactor != 0)
            debug_log << phi << " " << d << " " << scalefactor << std::endl;
      }
   }
   debug_log.close();
   h2_scalefactor.Draw("COLZ");
   c1->SaveAs("Output_PDF/scalefactor_map2.pdf");
}

// Draw difference between interpolated and exact values on a given range
void LUTMaker::DrawDiffInterpolExact(const double &transvDiff, const double &RC,
                                     const double &drift, const double &phimin,
                                     const double &phimax, const double &dmin,
                                     const double &dmax)
{
   int nbins = 100;
   TH2F h2_diff(
      "h2_diff",
      Form("Relative difference between interpolated and exact values (%%, %d bins)",
           nbins),
      nbins, phimin, phimax, nbins, dmin, dmax);
   for (int i = 0; i < nbins; i++) {
      for (int j = 0; j < nbins; j++) {
         double phi = phimin + i * (phimax - phimin) / nbins;
         double d = dmin + j * (dmax - dmin) / nbins;
         double L = arr_length[(int)std::round(phi / sSTEP_PHI)]
                              [(int)std::round(d / sSTEP_IMPACT)];
         float scalefactor = GetFactorFromLUT(transvDiff, RC, drift, d, phi);
         if (L < 2) {
            h2_diff.SetBinContent(h2_diff.FindBin(phi, d), -10);
            continue;
         }
         fp_trackmodel->ComputeAmplitudeLoss(L, d, phi, v_RC[0], v_drift[0], v_Dt[0]);
         float scalefactor_exact = 1 / fp_trackmodel->GetAmplitudeLoss();
         float diff = (scalefactor - scalefactor_exact) / scalefactor_exact * 100;
         if (scalefactor_exact == 0 || scalefactor == 0)
            continue;
         h2_diff.SetBinContent(h2_diff.FindBin(phi, d), diff);
      }
   }
   h2_diff.SetMaximum(8);
   h2_diff.SetMinimum(-8);
   h2_diff.Draw("colz");
   c1->SaveAs(Form("Output_PDF/diff_interpol_exact_%dphi%d_%dd%d_%dbins.pdf", (int)phimin,
                   (int)phimax, (int)dmin, (int)dmax, nbins));
}