#include "Process.h"
#include "Misc_Functions.h"

#include <fstream>

ClassImp(Process)

   // Constructor
   Process::Process()
{

   std::vector<int> channel2iD = {24, 30, 28, 19, 21, 13, 9,  2,  26, 17, 23, 29, 1,  10, 11, 3,   /*bottom HATPC*/
                                  47, 16, 14, 15, 42, 45, 37, 36, 20, 38, 7,  44, 43, 39, 41, 46}; /*top	HATPC*/

   for (int i = 0; i < 32; i++) {
      vmod_fph1f_XP.push_back(new TH1F(Form("fph1f_XP_%d", i),
                                       Form("Energy loss in ERAM %i;dE/dx (ADC counts/cm);Count", channel2iD[i]), 100,
                                       0, dEdxmax));
      vmod_fph1f_WF.push_back(new TH1F(Form("fph1f_WF_%d", i),
                                       Form("Energy loss in ERAM %i;dE/dx (ADC counts/cm);Count", channel2iD[i]), 100,
                                       0, dEdxmax));
      vmod_fph2f_XPtmean.push_back(
         new TH2F(Form("fph2f_XPtmean_%d", i),
                  Form("Energy loss (XP) vs mean time in ERAM %i;mean time (ns);dE/dx (ADC counts/cm)", channel2iD[i]),
                  510, 0, 510, 100, 0, dEdxmax));
      vmod_fph2f_XPdrift.push_back(new TH2F(
         Form("fph2f_XPdrift_%d", i),
         Form("Energy loss (XP) vs drift time in ERAM %i;drift time (timebins);dE/dx (ADC counts/cm)", channel2iD[i]),
         510, 0, 510, 100, 0, dEdxmax));
   }

   for (int i = 0; i < nmombins; i++) {
      // int mommin = i * mombinwidth - momrange;
      // int mommax = (i + 1) * mombinwidth - momrange;
      int mommin = 0;
      int mommax = momrange;
      vmom_fph1f_WF.push_back(new TH1F(Form("fph1f_mom_WF_%d_%d", mommin, mommax),
                                       Form("Energy loss | %d < p < %d; dE/dx (ADC counts/cm); Count", mommin, mommax),
                                       100, 0, dEdxmax));
      vmom_fph1f_XP.push_back(new TH1F(Form("fph1f_mom_XP_%d_%d", mommin, mommax),
                                       Form("Energy loss | %d < p < %d; dE/dx (ADC counts/cm); Count", mommin, mommax),
                                       100, 0, dEdxmax));
      vmom_fph1f_pullmu.push_back(new TH1F(Form("fph1f_pull_%d_%d", mommin, mommax),
                                           Form("Muon pull | %d < p < %d; Pull; Count", mommin, mommax), 100, -20, 20));
      vmom_fph1f_pullelec.push_back(new TH1F(Form("fph1f_pullelec_%d_%d", mommin, mommax),
                                             Form("Electron pull | %d < p < %d; Pull; Count", mommin, mommax), 100, -20,
                                             20));
   }

   for (int i = 0; i < nddbins; i++) {
      int ddmin = i * ddbinwidth;
      int ddmax = (i + 1) * ddbinwidth;
      vdd_fph1f_WF.push_back(new TH1F(
         Form("ph1f_dd_WF_%d_%d", ddmin, ddmax),
         Form("Energy loss | %d < drift distance < %d; dE/dx (ADC counts/cm); Count", ddmin, ddmax), 100, 0, dEdxmax));
      vdd_fph1f_XP.push_back(new TH1F(
         Form("ph1f_dd_XP_%d_%d", ddmin, ddmax),
         Form("Energy loss | %d < drift distance < %d; dE/dx (ADC counts/cm); Count", ddmin, ddmax), 100, 0, dEdxmax));
   }

   for (int i = 0; i < nphibins; i++) {
      int phimin = i * phibinwidth - phirange;
      int phimax = (i + 1) * phibinwidth - phirange;
      vphi_fph1f_WF.push_back(new TH1F(
         Form("fph1f_phi_WF_%d_%d", phimin, phimax),
         Form("Energy loss | %d < #varphi < %d; dE/dx (ADC counts/cm); Count", phimin, phimax), 100, 0, dEdxmax));
      vphi_fph1f_XP.push_back(new TH1F(
         Form("fph1f_phi_XP_%d_%d", phimin, phimax),
         Form("Energy loss | %d < #varphi < %d; dE/dx (ADC counts/cm); Count", phimin, phimax), 100, 0, dEdxmax));
   }

   for (int i = 0; i < nthetabins; i++) {
      int thetamin = i * thetabinwidth - thetarange;
      int thetamax = (i + 1) * thetabinwidth - thetarange;
      vtheta_fph1f_WF.push_back(new TH1F(
         Form("fph1f_theta_WF_%d_%d", thetamin, thetamax),
         Form("Energy loss | %d < #theta < %d; dE/dx (ADC counts/cm); Count", thetamin, thetamax), 100, 0, dEdxmax));
      vtheta_fph1f_XP.push_back(new TH1F(
         Form("fph1f_theta_XP_%d_%d", thetamin, thetamax),
         Form("Energy loss | %d < #theta < %d; dE/dx (ADC counts/cm); Count", thetamin, thetamax), 100, 0, dEdxmax));
      vtheta_fph1f_pullmu.push_back(new TH1F(Form("fph1f_pullmu_%d_%d", thetamin, thetamax),
                                             Form("Muon pull | %d < #theta < %d; Pull; Count", thetamin, thetamax), 100,
                                             -20, 20));
      vtheta_fph1f_pullelec.push_back(
         new TH1F(Form("fph1f_pullelec_%d_%d", thetamin, thetamax),
                  Form("Electron pull | %d < #theta < %d; Pull; Count", thetamin, thetamax), 100, -20, 20));
   }
}

// Destructor
Process::~Process()
{
   // Don't delete pointers to TFile and TTree, ROOT handles them

   // Delete single pointers
   delete fptf1_WF;
   delete fptf1_XP;

   delete ptge_mom_mean_WF;
   delete ptge_mom_mean_XP;
   delete ptge_mom_reso_WF;
   delete ptge_mom_reso_XP;

   delete ptge_mom_mean_pullmu;
   delete ptge_mom_mean_pullelec;
   delete ptge_mom_std_pullmu;
   delete ptge_mom_std_pullelec;

   delete ptge_dd_mean_WF;
   delete ptge_dd_mean_XP;
   delete ptge_dd_reso_WF;
   delete ptge_dd_reso_XP;

   delete ptge_phi_mean_WF;
   delete ptge_phi_mean_XP;
   delete ptge_phi_reso_WF;
   delete ptge_phi_reso_XP;

   delete ptge_theta_mean_WF;
   delete ptge_theta_mean_XP;
   delete ptge_theta_reso_WF;
   delete ptge_theta_reso_XP;

   delete ptge_theta_mean_pullmu;
   delete ptge_theta_mean_pullelec;
   delete ptge_theta_std_pullmu;
   delete ptge_theta_std_pullelec;

   delete fph1f_WF;
   delete fph1f_XP;
   delete fph2f_WFXP;
   delete fph2f_XPdrift;
   delete fph2f_XPlen;
   delete fph2f_XPphi;
   delete fph2f_XPtheta;
   delete fph2f_XPmom;
   delete fph2f_WFmom;
   delete fph1i_mom;
   delete fph1i_mom_tHAT;
   delete fph1i_mom_bHAT;
   delete fph2f_XZ;
   delete fph1f_phi;
   delete fph1f_theta;
   delete fph2f_phitheta;
   delete fph1f_trklen;
   delete fph1f_chi2;
   delete fph2f_chi2mom;
   delete fph1i_tminBotCath;
   delete fph1i_tminEP0;
   delete fph1i_tminEP1;
   delete fph1i_tmaxBotCath;
   delete fph1i_tmaxEP0;
   delete fph1i_tmaxEP1;
   delete fph1i_tminTopCath;
   delete fph1i_tminEP2;
   delete fph1i_tminEP3;
   delete fph1i_tmaxTopCath;
   delete fph1i_tmaxEP2;
   delete fph1i_tmaxEP3;
   delete fph2f_pullelecmu;
   delete fph2f_chi2ndfphi;
   delete fph2f_momtheta;
   delete fph2f_momR;
   delete fph2f_chi2ndfR;
   delete fph2f_lentheta;
   delete fph2f_lenphi;

   // Clear vector pointers
   for (auto ptr : vmod_fph1f_XP)
      delete ptr;
   for (auto ptr : vmod_fph1f_WF)
      delete ptr;
   for (auto ptr : vmod_fph2f_XPtmean)
      delete ptr;
   for (auto ptr : vmod_fph2f_XPdrift)
      delete ptr;
   for (auto ptr : vmom_fph1f_WF)
      delete ptr;
   for (auto ptr : vmom_fph1f_XP)
      delete ptr;
   for (auto ptr : vmom_fph1f_pullmu)
      delete ptr;
   for (auto ptr : vmom_fph1f_pullelec)
      delete ptr;
   for (auto ptr : vdd_fph1f_WF)
      delete ptr;
   for (auto ptr : vdd_fph1f_XP)
      delete ptr;
   for (auto ptr : vphi_fph1f_WF)
      delete ptr;
   for (auto ptr : vphi_fph1f_XP)
      delete ptr;
   for (auto ptr : vtheta_fph1f_WF)
      delete ptr;
   for (auto ptr : vtheta_fph1f_XP)
      delete ptr;
   for (auto ptr : vtheta_fph1f_pullmu)
      delete ptr;
   for (auto ptr : vtheta_fph1f_pullelec)
      delete ptr;

   // Clear vectors after deleting
   vmod_fph1f_XP.clear();
   vmod_fph1f_WF.clear();
   vmod_fph2f_XPtmean.clear();
   vmod_fph2f_XPdrift.clear();
   vmom_fph1f_WF.clear();
   vmom_fph1f_XP.clear();
   vmom_fph1f_pullmu.clear();
   vmom_fph1f_pullelec.clear();
   vdd_fph1f_WF.clear();
   vdd_fph1f_XP.clear();
   vphi_fph1f_WF.clear();
   vphi_fph1f_XP.clear();
   vtheta_fph1f_WF.clear();
   vtheta_fph1f_XP.clear();
   vtheta_fph1f_pullmu.clear();
   vtheta_fph1f_pullelec.clear();
}

void Process::SetCuts()
{
   // nclmin = 32;
   // fcutslist += ("_" + std::to_string(nclmin) + "ncl");
   // apmcutlow = 2, apmcuthigh = 4;
   // fcutslist += ("_" + std::to_string(apmcutlow) + "apm" + std::to_string(apmcuthigh));
   // momcutlow = 1e2;
   // fcutslist += ("_" + std::to_string(momcutlow) + "mom");
   // tcutmin = 0, tcutmax = 75;
   // fcutslist += ("_" + std::to_string(tcutmin) + "tmin" + std::to_string(tcutmax));
   // nclmin = 50, nclmax = 150;
   // fcutslist += ("_" + std::to_string(nclmin) + "ncl" + std::to_string(nclmax));
   // dxmin = 50, ddEdxmax = 150;
   // fcutslist += ("_" + std::to_string(dxmin) + "dx" + std::to_string(ddEdxmax));
   // pullmumax = 2;
   // fcutslist += ("_pullmu" + std::to_string(pullmumax));
   // phimin = 80;
   // phimax = 90;
   // fcutslist += ("_" + std::to_string(phimin) + "phi" + std::to_string(phimax));
   // thetamin = 0;
   // thetamax = 10;
   // fcutslist += ("_" + std::to_string(thetamin) + "theta" + std::to_string(thetamax));
   // chi2max = 5;
   // fcutslist += ("_chi2ndf" + std::to_string(chi2max));
   // dxmin = 25;
   // fcutslist += ("_" + std::to_string(dxmin) + "dx");
   // momcutlow = 100, momcuthigh = 600;
   // fcutslist += ("_" + std::to_string(momcutlow) + "mom" + std::to_string(momcuthigh));
   fcutslist += "_dir1>0flip";
   // hat = -1;
   // fcutslist += std::string("_") + (hat == -1 ? "bHAT" : "tHAT");

   // Output files
   ffileName += fcutslist;
   fRealpathROOT = foutputROOTfolder + "/" + ffileName + ".root";
}

void Process::Run()
{

   // Input file
   TFile *fpInputFile = TFile::Open(finputFile.c_str());
   TTree *fpInputTree = (TTree *)fpInputFile->Get("outTree");
   int nentries = fpInputTree->GetEntries();

   // Tree branches
   fpInputTree->SetBranchAddress("dEdx_WF", &wf);
   fpInputTree->SetBranchAddress("dEdx_XP", &xp);
   fpInputTree->SetBranchAddress("track_length", &dx);
   fpInputTree->SetBranchAddress("pos", &pos);
   fpInputTree->SetBranchAddress("dir", &dir);
   fpInputTree->SetBranchAddress("chi2", &chi2);
   fpInputTree->SetBranchAddress("curv", &curv);
   fpInputTree->SetBranchAddress("NDF", &NDF);
   fpInputTree->SetBranchAddress("mom", &mom_og);
   fpInputTree->SetBranchAddress("avg_pad_mult", &APM);
   fpInputTree->SetBranchAddress("starttime", &start_time);
   fpInputTree->SetBranchAddress("mean_time", &mean_time);
   fpInputTree->SetBranchAddress("endtime", &end_time);
   fpInputTree->SetBranchAddress("ncl", &ncl);
   fpInputTree->SetBranchAddress("EndPlate_ID", &endplate);
   fpInputTree->SetBranchAddress("eram_channel", &eram_channel);
   fpInputTree->SetBranchAddress("eram_ID", &eram_ID);
   fpInputTree->SetBranchAddress("pull_muon", &pull_muon);
   fpInputTree->SetBranchAddress("pull_ele", &pull_ele);
   fpInputTree->SetBranchAddress("pull_prot", &pull_proton);

   // Debugging variables
   int dirYnegbHAT = 0, dirYposbHAT = 0, dirYnegtHAT = 0, dirYpostHAT = 0;
   int nbotneg = 0, nbotpos = 0, ntopneg = 0, ntoppos = 0;

   for (int i = 0; i < nentries; i++) {
      fpInputTree->GetEntry(i);
      mom = mom_og;
      // FOR COSMICS if the track is going up, the sign must be flipped
      if (fRealpathROOT.find("cosmic") != std::string::npos and fcutslist.find("dir1>0flip") != std::string::npos)
         dir[1] > 0 ? mom = -mom_og : mom = mom_og;
      float dirxy = TMath::Sqrt(dir[0] * dir[0] + dir[1] * dir[1]);
      float diryz = TMath::Sqrt(dir[1] * dir[1] + dir[2] * dir[2]);
      phi = TMath::ATan(dirxy / dir[2]) * 180 / TMath::Pi();
      theta = TMath::ATan(dir[0] / diryz) * 180 / TMath::Pi();
      if (std::isnan(mom))
         continue;

      if (nclmin > ncl or ncl > nclmax)
         continue;
      if (dxmin > dx / 10 or ddEdxmax < dx / 10)
         continue;
      if (fabs(phi) < phimin or phimax < fabs(phi))
         continue;
      if (fabs(theta) < thetamin or thetamax < fabs(theta))
         continue;
      if (dx < dxmin)
         continue;
      if (APM < apmcutlow or APM > apmcuthigh)
         continue;
      if (fabs(mom) < momcutlow or fabs(mom) > momcuthigh)
         continue;
      if (mean_time < tcutmin or mean_time > tcutmax)
         continue;
      if (chi2 / NDF > chi2max)
         continue;
      if (hat != 0 and sign(hat) != sign(pos[1]))
         continue;
      // if(fabs(pull_muon) > fabs(pull_ele)) continue;
      if (fabs(pull_muon) > pullmumax)
         continue;

      fph1f_WF->Fill(wf / 1.019);
      fph1f_XP->Fill(xp);
      fph2f_WFXP->Fill(wf / 1.019, xp);
      fph2f_XPlen->Fill(dx / 10, xp);
      fph2f_XPdrift->Fill(mean_time, xp);
      fph2f_XPphi->Fill(phi, xp);
      fph2f_XPtheta->Fill(theta, xp);

      vmod_fph2f_XPdrift[eram_channel]->Fill(mean_time, xp);
      vmod_fph1f_WF[eram_channel]->Fill(wf / 1.019);
      vmod_fph1f_XP[eram_channel]->Fill(xp);
      vmod_fph2f_XPtmean[eram_channel]->Fill(mean_time, xp);

      fph2f_XZ->Fill(dir[0], dir[2]);
      fph1f_phi->Fill(phi);
      fph1f_theta->Fill(theta);
      fph2f_phitheta->Fill(theta, phi);
      fph1f_trklen->Fill(dx / 10);
      fph1f_chi2->Fill(chi2 / NDF);
      fph2f_chi2mom->Fill(mom, chi2 / NDF);

      fph1i_mom->Fill(mom);
      if (pos[1] < 0)
         fph1i_mom_bHAT->Fill(mom);
      if (pos[1] > 0)
         fph1i_mom_tHAT->Fill(mom);
      fph2f_WFmom->Fill(mom, wf / 1.019);
      fph2f_XPmom->Fill(mom, xp);

      momindex = (int)std::round(fabs(mom) / mombinwidth);
      if (fabs(mom) < momrange) {
         vmom_fph1f_WF[momindex]->Fill(wf / 1.019);
         vmom_fph1f_XP[momindex]->Fill(xp);
         vmom_fph1f_pullmu[momindex]->Fill(pull_muon);
         vmom_fph1f_pullelec[momindex]->Fill(pull_ele);
      }

      ddindex = (int)std::round(mean_time / ddbinwidth);
      if (ddindex < nddbins) {
         vdd_fph1f_WF[ddindex]->Fill(wf / 1.019);
         vdd_fph1f_XP[ddindex]->Fill(xp);
      }

      phiindex = (int)std::round(phi / phibinwidth) + nphibins / 2;
      if (phiindex < nphibins) {
         vphi_fph1f_WF[phiindex]->Fill(wf / 1.019);
         vphi_fph1f_XP[phiindex]->Fill(xp);
      }

      thetaindex = (int)std::round(theta / thetabinwidth) + nthetabins / 2;
      if (thetaindex < nthetabins) {
         vtheta_fph1f_WF[thetaindex]->Fill(wf / 1.019);
         vtheta_fph1f_XP[thetaindex]->Fill(xp);
         vtheta_fph1f_pullmu[thetaindex]->Fill(pull_muon);
         vtheta_fph1f_pullelec[thetaindex]->Fill(pull_ele);
      }

      // Drift time in ERAM modules
      if (eram_channel < 16) {
         if (endplate == 8) {
            fph1i_tminEP0->Fill(start_time);
            fph1i_tmaxEP0->Fill(end_time);
         }
         if (endplate == 4) {
            fph1i_tminEP1->Fill(start_time);
            fph1i_tmaxEP1->Fill(end_time);
         }
         if (endplate == 12) {
            fph1i_tminBotCath->Fill(start_time);
            fph1i_tmaxBotCath->Fill(end_time);
         }
      } else {
         if (endplate == 8) {
            fph1i_tminEP2->Fill(start_time);
            fph1i_tmaxEP2->Fill(end_time);
         }
         if (endplate == 4) {
            fph1i_tminEP3->Fill(start_time);
            fph1i_tmaxEP3->Fill(end_time);
         }
         if (endplate == 12) {
            fph1i_tminTopCath->Fill(start_time);
            fph1i_tmaxTopCath->Fill(end_time);
         }
      }

      // Debug
      if (dir[1] < 0 and pos[1] < 0)
         dirYnegbHAT++;
      if (dir[1] > 0 and pos[1] < 0)
         dirYposbHAT++;
      if (dir[1] < 0 and pos[1] > 0)
         dirYnegtHAT++;
      if (dir[1] > 0 and pos[1] > 0)
         dirYpostHAT++;

      if (pos[1] < 0 and sign(mom) < 0)
         nbotneg++;
      if (pos[1] < 0 and sign(mom) > 0)
         nbotpos++;
      if (pos[1] > 0 and sign(mom) < 0)
         ntopneg++;
      if (pos[1] > 0 and sign(mom) > 0)
         ntoppos++;
      fph2f_pullelecmu->Fill(pull_muon, pull_ele);
      fph2f_chi2ndfphi->Fill(phi, chi2 / NDF);
      fph2f_momtheta->Fill(theta, mom);
      fph2f_momphi->Fill(phi, mom);
      fph2f_momR->Fill(1 / curv, mom);
      fph2f_chi2ndfR->Fill(1 / curv, chi2 / NDF);
      fph2f_lentheta->Fill(theta, dx / 10);
      fph2f_lenphi->Fill(phi, dx / 10);
      fph1f_pullmu->Fill(pull_muon);
      fph1f_pullelec->Fill(pull_ele);
      fph1f_pullproton->Fill(pull_proton);
   }
   // std::cout << "Y direction:" << std::endl;
   // std::cout << "tHAT: negative => " << dirYnegtHAT << " (" << dirYnegtHAT*100.0/(dirYnegtHAT+dirYpostHAT) << "%)
   // positive => " << dirYpostHAT << " (" << dirYpostHAT*100.0/(dirYnegtHAT+dirYpostHAT) << "%)" << std::endl;
   // std::cout << "bHAT: negative => " << dirYnegbHAT << " (" << dirYnegbHAT*100.0/(dirYnegbHAT+dirYposbHAT) << "%)
   // positive => " << dirYposbHAT << " (" << dirYposbHAT*100.0/(dirYnegbHAT+dirYposbHAT) << "%)" << std::endl;
   // std::cout << "Total: negative => " << dirYnegbHAT+dirYnegtHAT << " (" <<
   // (dirYnegbHAT+dirYnegtHAT)*100.0/(dirYnegbHAT+dirYposbHAT+dirYnegtHAT+dirYpostHAT) << "%) positive => " <<
   // dirYposbHAT+dirYpostHAT << " (" <<
   // (dirYposbHAT+dirYpostHAT)*100.0/(dirYnegbHAT+dirYposbHAT+dirYnegtHAT+dirYpostHAT) << "%)" << std::endl;

   // std::cout << "Momentum:" << std::endl;
   // std::cout << "Negative: " << ntopneg << " (top) => " << nbotneg << " (bottom) => top-bottom = " << ntopneg-nbotneg
   // << std::endl; std::cout << "Positive: " << ntoppos << " (top) => " << nbotpos << " (bottom) => top-bottom = " <<
   // ntoppos-nbotpos << std::endl;

   // TGraph filling
   int ivalid = 0;
   // Momentum
   for (int i = 0; i < nmombins; i++) {
      int nentries_here = vmom_fph1f_WF[i]->GetEntries();
      if (nentries_here < 50)
         continue;
      fptf1_WF = Fit1Gauss(vmom_fph1f_WF[i]);
      fptf1_XP = Fit1Gauss(vmom_fph1f_XP[i]);
      fptf1_pullmu = Fit1Gauss(vmom_fph1f_pullmu[i]);
      fptf1_pullelec = Fit1Gauss(vmom_fph1f_pullelec[i]);

      float mean_WF = fptf1_WF->GetParameter(1);
      float mean_XP = fptf1_XP->GetParameter(1);
      float mean_pullmu = fptf1_pullmu->GetParameter(1);
      float mean_pullelec = fptf1_pullelec->GetParameter(1);
      float dmean_WF = fptf1_WF->GetParError(1);
      float dmean_XP = fptf1_XP->GetParError(1);
      float dmean_pullmu = fptf1_pullmu->GetParError(1);
      float dmean_pullelec = fptf1_pullelec->GetParError(1);

      float std_WF = fptf1_WF->GetParameter(2);
      float std_XP = fptf1_XP->GetParameter(2);
      float std_pullmu = fptf1_pullmu->GetParameter(2);
      float std_pullelec = fptf1_pullelec->GetParameter(2);
      float dstd_pullmu = fptf1_pullmu->GetParError(2);
      float dstd_pullelec = fptf1_pullelec->GetParError(2);

      float reso_WF = std_WF / mean_WF * 100;
      float reso_XP = std_XP / mean_XP * 100;
      float dreso_WF = GetResoError(fptf1_WF);
      float dreso_XP = GetResoError(fptf1_XP);

      if (mean_WF == 0 || mean_XP == 0)
         continue;

      ptge_mom_mean_WF->SetPoint(ivalid, i * mombinwidth, mean_WF);
      ptge_mom_mean_WF->SetPointError(ivalid, mombinwidth / 2, dmean_WF);
      ptge_mom_mean_XP->SetPoint(ivalid, i * mombinwidth, mean_XP);
      ptge_mom_mean_XP->SetPointError(ivalid, mombinwidth / 2, dmean_XP);
      ptge_mom_mean_pullmu->SetPoint(ivalid, i * mombinwidth - momrange, mean_pullmu);
      ptge_mom_mean_pullmu->SetPointError(ivalid, mombinwidth / 2, dmean_pullmu);
      ptge_mom_mean_pullelec->SetPoint(ivalid, i * mombinwidth - momrange, mean_pullelec);
      ptge_mom_mean_pullelec->SetPointError(ivalid, mombinwidth / 2, dmean_pullelec);

      ptge_mom_reso_WF->SetPoint(ivalid, (i * mombinwidth)/1000, reso_WF);
      ptge_mom_reso_WF->SetPointError(ivalid, 0, 0);
      ptge_mom_reso_XP->SetPoint(ivalid, (i * mombinwidth)/1000, reso_XP);
      ptge_mom_reso_XP->SetPointError(ivalid, 0, 0);
      ptge_mom_std_pullmu->SetPoint(ivalid, i * mombinwidth - momrange, std_pullmu);
      ptge_mom_std_pullmu->SetPointError(ivalid, mombinwidth / 2, dstd_pullmu);
      ptge_mom_std_pullelec->SetPoint(ivalid, i * mombinwidth - momrange, std_pullelec);
      ptge_mom_std_pullelec->SetPointError(ivalid, mombinwidth / 2, dstd_pullelec);
      ivalid++;
   }

   // Drift distance
   ivalid = 0;
   for (int i = 0; i < nddbins; i++) {
      int nentries_here = vdd_fph1f_WF[i]->GetEntries();
      if (nentries_here < 50)
         continue;
      fptf1_WF = Fit1Gauss(vdd_fph1f_WF[i]);
      fptf1_XP = Fit1Gauss(vdd_fph1f_XP[i]);

      float mean_WF = fptf1_WF->GetParameter(1);
      float mean_XP = fptf1_XP->GetParameter(1);
      float dmean_WF = fptf1_WF->GetParError(1);
      float dmean_XP = fptf1_XP->GetParError(1);

      float std_WF = fptf1_WF->GetParameter(2);
      float std_XP = fptf1_XP->GetParameter(2);

      float reso_WF = std_WF / mean_WF * 100;
      float reso_XP = std_XP / mean_XP * 100;
      float dreso_WF = GetResoError(fptf1_WF);
      float dreso_XP = GetResoError(fptf1_XP);

      if (mean_WF == 0 || mean_XP == 0)
         continue;

      ptge_dd_mean_WF->SetPoint(ivalid, i * ddbinwidth, mean_WF);
      ptge_dd_mean_WF->SetPointError(ivalid, ddbinwidth / 2, dmean_WF);
      ptge_dd_mean_XP->SetPoint(ivalid, i * ddbinwidth, mean_XP);
      ptge_dd_mean_XP->SetPointError(ivalid, ddbinwidth / 2, dmean_XP);

      ptge_dd_reso_WF->SetPoint(ivalid, i * ddbinwidth, reso_WF);
      ptge_dd_reso_WF->SetPointError(ivalid, ddbinwidth / 2, dreso_WF);
      ptge_dd_reso_XP->SetPoint(ivalid, i * ddbinwidth, reso_XP);
      ptge_dd_reso_XP->SetPointError(ivalid, ddbinwidth / 2, dreso_XP);
      ivalid++;
   }

   // Phi angle
   ivalid = 0;
   for (int i = 0; i < nphibins; i++) {
      int nentries_here = vphi_fph1f_WF[i]->GetEntries();
      if (nentries_here < 50)
         continue;
      fptf1_WF = Fit1Gauss(vphi_fph1f_WF[i]);
      fptf1_XP = Fit1Gauss(vphi_fph1f_XP[i]);

      float mean_WF = fptf1_WF->GetParameter(1);
      float mean_XP = fptf1_XP->GetParameter(1);
      float dmean_WF = fptf1_WF->GetParError(1);
      float dmean_XP = fptf1_XP->GetParError(1);

      float std_WF = fptf1_WF->GetParameter(2);
      float std_XP = fptf1_XP->GetParameter(2);

      float reso_WF = std_WF / mean_WF * 100;
      float reso_XP = std_XP / mean_XP * 100;
      float dreso_WF = GetResoError(fptf1_WF);
      float dreso_XP = GetResoError(fptf1_XP);

      if (mean_WF == 0 || mean_XP == 0)
         continue;

      ptge_phi_mean_WF->SetPoint(ivalid, i * phibinwidth - phirange, mean_WF);
      ptge_phi_mean_WF->SetPointError(ivalid, phibinwidth / 2, dmean_WF);
      ptge_phi_mean_XP->SetPoint(ivalid, i * phibinwidth - phirange, mean_XP);
      ptge_phi_mean_XP->SetPointError(ivalid, phibinwidth / 2, dmean_XP);

      ptge_phi_reso_WF->SetPoint(ivalid, i * phibinwidth - phirange, reso_WF);
      ptge_phi_reso_WF->SetPointError(ivalid, phibinwidth / 2, dreso_WF);
      ptge_phi_reso_XP->SetPoint(ivalid, i * phibinwidth - phirange, reso_XP);
      ptge_phi_reso_XP->SetPointError(ivalid, phibinwidth / 2, dreso_XP);
      ivalid++;
   }

   // Theta angle
   ivalid = 0;
   for (int i = 0; i < nthetabins; i++) {
      int nentries_here = vtheta_fph1f_WF[i]->GetEntries();
      if (nentries_here < 50)
         continue;
      fptf1_WF = Fit1Gauss(vtheta_fph1f_WF[i]);
      fptf1_XP = Fit1Gauss(vtheta_fph1f_XP[i]);
      fptf1_pullmu = Fit1Gauss(vtheta_fph1f_pullmu[i]);
      fptf1_pullelec = Fit1Gauss(vtheta_fph1f_pullelec[i]);

      float mean_WF = fptf1_WF->GetParameter(1);
      float mean_XP = fptf1_XP->GetParameter(1);
      float mean_pullmu = fptf1_pullmu->GetParameter(1);
      float mean_pullelec = fptf1_pullelec->GetParameter(1);
      float dmean_WF = fptf1_WF->GetParError(1);
      float dmean_XP = fptf1_XP->GetParError(1);
      float dmean_pullmu = fptf1_pullmu->GetParError(1);
      float dmean_pullelec = fptf1_pullelec->GetParError(1);

      float std_WF = fptf1_WF->GetParameter(2);
      float std_XP = fptf1_XP->GetParameter(2);
      float std_pullmu = fptf1_pullmu->GetParameter(2);
      float std_pullelec = fptf1_pullelec->GetParameter(2);
      float dstd_pullmu = fptf1_pullmu->GetParError(2);
      float dstd_pullelec = fptf1_pullelec->GetParError(2);

      float reso_WF = std_WF / mean_WF * 100;
      float reso_XP = std_XP / mean_XP * 100;
      float dreso_WF = GetResoError(fptf1_WF);
      float dreso_XP = GetResoError(fptf1_XP);

      if (mean_WF == 0 || mean_XP == 0)
         continue;

      ptge_theta_mean_WF->SetPoint(ivalid, i * thetabinwidth - thetarange, mean_WF);
      ptge_theta_mean_WF->SetPointError(ivalid, thetabinwidth / 2, dmean_WF);
      ptge_theta_mean_XP->SetPoint(ivalid, i * thetabinwidth - thetarange, mean_XP);
      ptge_theta_mean_XP->SetPointError(ivalid, thetabinwidth / 2, dmean_XP);
      ptge_theta_mean_pullmu->SetPoint(ivalid, i * thetabinwidth - thetarange, mean_pullmu);
      ptge_theta_mean_pullmu->SetPointError(ivalid, thetabinwidth / 2, dmean_pullmu);
      ptge_theta_mean_pullelec->SetPoint(ivalid, i * thetabinwidth - thetarange, mean_pullelec);
      ptge_theta_mean_pullelec->SetPointError(ivalid, thetabinwidth / 2, dmean_pullelec);

      ptge_theta_reso_WF->SetPoint(ivalid, i * thetabinwidth - thetarange, reso_WF);
      ptge_theta_reso_WF->SetPointError(ivalid, thetabinwidth / 2, dreso_WF);
      ptge_theta_reso_XP->SetPoint(ivalid, i * thetabinwidth - thetarange, reso_XP);
      ptge_theta_reso_XP->SetPointError(ivalid, thetabinwidth / 2, dreso_XP);
      ptge_theta_std_pullmu->SetPoint(ivalid, i * thetabinwidth - thetarange, std_pullmu);
      ptge_theta_std_pullmu->SetPointError(ivalid, thetabinwidth / 2, dstd_pullmu);
      ptge_theta_std_pullelec->SetPoint(ivalid, i * thetabinwidth - thetarange, std_pullelec);
      ptge_theta_std_pullelec->SetPointError(ivalid, thetabinwidth / 2, dstd_pullelec);
      ivalid++;
   }
   fpInputFile->Close();
}
