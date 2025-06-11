#include "Process.h"
#include "Misc_Functions.h"

#include <fstream>

ClassImp(Process)

   // Constructor
   Process::Process()
{

   std::vector<int> channel2iD = {
      24, 30, 28, 19, 21, 13, 9,  2,  26, 17, 23, 29, 1,  10, 11, 3,   /*bottom HATPC*/
      47, 16, 14, 15, 42, 45, 37, 36, 20, 38, 7,  44, 43, 39, 41, 46}; /*top	HATPC*/

   // dEdx per ERAM module
   for (int i = 0; i < 32; i++) {
      vmod_fph1f_XP.push_back(new TH1F(
         Form("fph1f_XP_%d", i),
         Form("Energy loss in ERAM %i;dE/dx (ADC counts/cm);Count", channel2iD[i]), 100,
         0, dEdxmax));
      vmod_fph1f_WF.push_back(new TH1F(
         Form("fph1f_WF_%d", i),
         Form("Energy loss in ERAM %i;dE/dx (ADC counts/cm);Count", channel2iD[i]), 100,
         0, dEdxmax));
      vmod_fph2f_XPdrift.push_back(
         new TH2F(Form("fph2f_XPdrift_%d", i),
                  Form("Energy loss (XP) vs drift time in ERAM %i;drift time "
                       "(timebins);dE/dx (ADC counts/cm)",
                       channel2iD[i]),
                  510, 0, 510, 100, 0, dEdxmax));
   }

   // dEdx vs absolute momentum (mean)
   for (int i = 0; i < nabsmommeanbins; i++) {
      int absmommeanmin = i * absmommeanbinwidth;
      int absmommeanmax = (i + 1) * absmommeanbinwidth;
      vabsmommean_fph1f_WF.push_back(
         new TH1F(Form("fph1f_absmommean_WF_%d_%d", absmommeanmin, absmommeanmax),
                  Form("Energy loss | %d < p < %d; dE/dx (ADC counts/cm); Count",
                       absmommeanmin, absmommeanmax),
                  100, 0, dEdxmax));
      vabsmommean_fph1f_XP.push_back(
         new TH1F(Form("fph1f_absmommean_XP_%d_%d", absmommeanmin, absmommeanmax),
                  Form("Energy loss | %d < p < %d; dE/dx (ADC counts/cm); Count",
                       absmommeanmin, absmommeanmax),
                  100, 0, dEdxmax));
   }

   // dEdx vs absolute momentum (resolution)
   for (int i = 0; i < nabsmomresobins; i++) {
      int absmomresomin = i * absmomresobinwidth;
      int absmomresomax = (i + 1) * absmomresobinwidth;
      vabsmomreso_fph1f_WF.push_back(
         new TH1F(Form("fph1f_absmomreso_WF_%d_%d", absmomresomin, absmomresomax),
                  Form("Energy loss | %d < p < %d; dE/dx (ADC counts/cm); Count",
                       absmomresomin, absmomresomax),
                  100, 0, dEdxmax));
      vabsmomreso_fph1f_XP.push_back(
         new TH1F(Form("fph1f_absmomreso_XP_%d_%d", absmomresomin, absmomresomax),
                  Form("Energy loss | %d < p < %d; dE/dx (ADC counts/cm); Count",
                       absmomresomin, absmomresomax),
                  100, 0, dEdxmax));
   }

   // dEdx and pulls vs momentum
   for (int i = 0; i < nmombins; i++) {
      int mommin = i * mombinwidth - momrange;
      int mommax = (i + 1) * mombinwidth - momrange;
      vmom_fph1f_WF.push_back(new TH1F(
         Form("fph1f_mom_WF_%d_%d", mommin, mommax),
         Form("Energy loss | %d < p < %d; dE/dx (ADC counts/cm); Count", mommin, mommax),
         100, 0, dEdxmax));
      vmom_fph1f_XP.push_back(new TH1F(
         Form("fph1f_mom_XP_%d_%d", mommin, mommax),
         Form("Energy loss | %d < p < %d; dE/dx (ADC counts/cm); Count", mommin, mommax),
         100, 0, dEdxmax));
      vmom_fph1f_pullmu.push_back(new TH1F(
         Form("fph1f_pull_%d_%d", mommin, mommax),
         Form("Muon pull | %d < p < %d; Pull; Count", mommin, mommax), 100, -20, 20));
      vmom_fph1f_pullelec.push_back(new TH1F(
         Form("fph1f_pullelec_%d_%d", mommin, mommax),
         Form("Electron pull | %d < p < %d; Pull; Count", mommin, mommax), 100, -20, 20));
   }

   // dEdx vs X position
   for (int i = 0; i < nxposbins; i++) {
      int Xmin = i * xposbinwidth - xposrange;
      int Xmax = (i + 1) * xposbinwidth - xposrange;
      vX_fph1f_WF.push_back(new TH1F(
         Form("fph1f_X_WF_%d_%d", Xmin, Xmax),
         Form("Energy loss | %d < X < %d; dE/dx (ADC counts/cm); Count", Xmin, Xmax), 100,
         0, dEdxmax));
      vX_fph1f_XP.push_back(new TH1F(
         Form("fph1f_X_XP_%d_%d", Xmin, Xmax),
         Form("Energy loss | %d < X < %d; dE/dx (ADC counts/cm); Count", Xmin, Xmax), 100,
         0, dEdxmax));
   }

   // dEdx vs drift distance (mean)
   for (int i = 0; i < nddmeanbins; i++) {
      int ddmin = i * ddmeanbinwidth;
      int ddmax = (i + 1) * ddmeanbinwidth;
      vddmean_fph1f_WF.push_back(new TH1F(
         Form("ph1f_dd_WF_%d_%d", ddmin, ddmax),
         Form("Energy loss | %d < drift distance < %d; dE/dx (ADC counts/cm); Count",
              ddmin, ddmax),
         100, 0, dEdxmax));
      vddmean_fph1f_XP.push_back(new TH1F(
         Form("ph1f_dd_XP_%d_%d", ddmin, ddmax),
         Form("Energy loss | %d < drift distance < %d; dE/dx (ADC counts/cm); Count",
              ddmin, ddmax),
         100, 0, dEdxmax));
   }

   // dEdx vs drift distance (resolution)
   for (int i = 0; i < nddresobins; i++) {
      int ddresomin = i * ddresobinwidth;
      int ddresomax = (i + 1) * ddresobinwidth;
      vddreso_fph1f_WF.push_back(new TH1F(
         Form("ph1f_ddreso_WF_%d_%d", ddresomin, ddresomax),
         Form("Energy loss | %d < drift distance < %d; dE/dx (ADC counts/cm); Count",
              ddresomin, ddresomax),
         100, 0, dEdxmax));
      vddreso_fph1f_XP.push_back(new TH1F(
         Form("ph1f_ddreso_XP_%d_%d", ddresomin, ddresomax),
         Form("Energy loss | %d < drift distance < %d; dE/dx (ADC counts/cm); Count",
              ddresomin, ddresomax),
         100, 0, dEdxmax));
   }

   // dEdx vs drift time
   for (int i = 0; i < ndtbins; i++) {
      int dtmin = i * dtbinwidth;
      int dtmax = (i + 1) * dtbinwidth;
      vdt_fph1f_WF.push_back(
         new TH1F(Form("ph1f_dt_WF_%d_%d", dtmin, dtmax),
                  Form("Energy loss | %d < drift time < %d; dE/dx (ADC counts/cm); Count",
                       dtmin, dtmax),
                  100, 0, dEdxmax));
      vdt_fph1f_XP.push_back(
         new TH1F(Form("ph1f_dt_XP_%d_%d", dtmin, dtmax),
                  Form("Energy loss | %d < drift time < %d; dE/dx (ADC counts/cm); Count",
                       dtmin, dtmax),
                  100, 0, dEdxmax));
   }

   // dEdx vs track length (mean)
   for (int i = 0; i < ntrklenmeanbins; i++) {
      int trklenmin = i * trklenmeanbinwidth;
      int trklenmax = (i + 1) * trklenmeanbinwidth;
      vtrklenmean_fph1f_WF.push_back(new TH1F(
         Form("ph1f_trklen_WF_%d_%d", trklenmin, trklenmax),
         Form("Energy loss | %d < track length < %d; dE/dx (ADC counts/cm); Count",
              trklenmin, trklenmax),
         100, 0, dEdxmax));
      vtrklenmean_fph1f_XP.push_back(new TH1F(
         Form("ph1f_trklen_XP_%d_%d", trklenmin, trklenmax),
         Form("Energy loss | %d < track length < %d; dE/dx (ADC counts/cm); Count",
              trklenmin, trklenmax),
         100, 0, dEdxmax));
   }

   // dEdx vs track length (resolution)
   for (int i = 0; i < ntrklenresobins; i++) {
      int trklenresomin = i * trklenresobinwidth;
      int trklenresomax = (i + 1) * trklenresobinwidth;
      vtrklenreso_fph1f_WF.push_back(new TH1F(
         Form("ph1f_trklenreso_WF_%d_%d", trklenresomin, trklenresomax),
         Form("Energy loss | %d < track length < %d; dE/dx (ADC counts/cm); Count",
              trklenresomin, trklenresomax),
         100, 0, dEdxmax));
      vtrklenreso_fph1f_XP.push_back(new TH1F(
         Form("ph1f_trklenreso_XP_%d_%d", trklenresomin, trklenresomax),
         Form("Energy loss | %d < track length < %d; dE/dx (ADC counts/cm); Count",
              trklenresomin, trklenresomax),
         100, 0, dEdxmax));
   }

   // dEdx vs absolute phi angle (mean)
   for (int i = 0; i < nabsphimeanbins; i++) {
      int absphimin = i * absphimeanbinwidth;
      int absphimax = (i + 1) * absphimeanbinwidth;
      vabsphimean_fph1f_WF.push_back(
         new TH1F(Form("fph1f_absphi_WF_%d_%d", absphimin, absphimax),
                  Form("Energy loss | %d < #varphi < %d; dE/dx (ADC counts/cm); Count",
                       absphimin, absphimax),
                  100, 0, dEdxmax));
      vabsphimean_fph1f_XP.push_back(
         new TH1F(Form("fph1f_absphi_XP_%d_%d", absphimin, absphimax),
                  Form("Energy loss | %d < #varphi < %d; dE/dx (ADC counts/cm); Count",
                       absphimin, absphimax),
                  100, 0, dEdxmax));
   }

   // dEdx vs absolute phi angle (resolution)
   for (int i = 0; i < nabsphiresobins; i++) {
      int absphiresomin = i * absphiresobinwidth;
      int absphiresomax = (i + 1) * absphiresobinwidth;
      vabsphireso_fph1f_WF.push_back(
         new TH1F(Form("fph1f_absphireso_WF_%d_%d", absphiresomin, absphiresomax),
                  Form("Energy loss | %d < #varphi < %d; dE/dx (ADC counts/cm); Count",
                       absphiresomin, absphiresomax),
                  100, 0, dEdxmax));
      vabsphireso_fph1f_XP.push_back(
         new TH1F(Form("fph1f_absphireso_XP_%d_%d", absphiresomin, absphiresomax),
                  Form("Energy loss | %d < #varphi < %d; dE/dx (ADC counts/cm); Count",
                       absphiresomin, absphiresomax),
                  100, 0, dEdxmax));
   }

   // dEdx vs phi angle
   for (int i = 0; i < nphibins; i++) {
      int phimin = i * phibinwidth - phirange;
      int phimax = (i + 1) * phibinwidth - phirange;
      vphi_fph1f_WF.push_back(
         new TH1F(Form("fph1f_phi_WF_%d_%d", phimin, phimax),
                  Form("Energy loss | %d < #varphi < %d; dE/dx (ADC counts/cm); Count",
                       phimin, phimax),
                  100, 0, dEdxmax));
      vphi_fph1f_XP.push_back(
         new TH1F(Form("fph1f_phi_XP_%d_%d", phimin, phimax),
                  Form("Energy loss | %d < #varphi < %d; dE/dx (ADC counts/cm); Count",
                       phimin, phimax),
                  100, 0, dEdxmax));
   }

   // dEdx vs absolute theta angle (mean)
   for (int i = 0; i < nabsthetameanbins; i++) {
      int absthetamin = i * absthetameanbinwidth;
      int absthetamax = (i + 1) * absthetameanbinwidth;
      vabsthetamean_fph1f_WF.push_back(
         new TH1F(Form("fph1f_abstheta_WF_%d_%d", absthetamin, absthetamax),
                  Form("Energy loss | %d < #theta < %d; dE/dx (ADC counts/cm); Count",
                       absthetamin, absthetamax),
                  100, 0, dEdxmax));
      vabsthetamean_fph1f_XP.push_back(
         new TH1F(Form("fph1f_abstheta_XP_%d_%d", absthetamin, absthetamax),
                  Form("Energy loss | %d < #theta < %d; dE/dx (ADC counts/cm); Count",
                       absthetamin, absthetamax),
                  100, 0, dEdxmax));
   }

   // dEdx vs absolute theta angle (resolution)
   for (int i = 0; i < nabsthetaresobins; i++) {
      int absthetaresomin = i * absthetaresobinwidth;
      int absthetaresomax = (i + 1) * absthetaresobinwidth;
      vabsthetareso_fph1f_WF.push_back(
         new TH1F(Form("fph1f_absthetareso_WF_%d_%d", absthetaresomin, absthetaresomax),
                  Form("Energy loss | %d < #theta < %d; dE/dx (ADC counts/cm); Count",
                       absthetaresomin, absthetaresomax),
                  100, 0, dEdxmax));
      vabsthetareso_fph1f_XP.push_back(
         new TH1F(Form("fph1f_absthetareso_XP_%d_%d", absthetaresomin, absthetaresomax),
                  Form("Energy loss | %d < #theta < %d; dE/dx (ADC counts/cm); Count",
                       absthetaresomin, absthetaresomax),
                  100, 0, dEdxmax));
   }

   // dEdx vs theta angle
   for (int i = 0; i < nthetabins; i++) {
      int thetamin = i * thetabinwidth - thetarange;
      int thetamax = (i + 1) * thetabinwidth - thetarange;
      vtheta_fph1f_WF.push_back(
         new TH1F(Form("fph1f_theta_WF_%d_%d", thetamin, thetamax),
                  Form("Energy loss | %d < #theta < %d; dE/dx (ADC counts/cm); Count",
                       thetamin, thetamax),
                  100, 0, dEdxmax));
      vtheta_fph1f_XP.push_back(
         new TH1F(Form("fph1f_theta_XP_%d_%d", thetamin, thetamax),
                  Form("Energy loss | %d < #theta < %d; dE/dx (ADC counts/cm); Count",
                       thetamin, thetamax),
                  100, 0, dEdxmax));
      vtheta_fph1f_pullmu.push_back(
         new TH1F(Form("fph1f_pullmu_%d_%d", thetamin, thetamax),
                  Form("Muon pull | %d < #theta < %d; Pull; Count", thetamin, thetamax),
                  100, -20, 20));
      vtheta_fph1f_pullelec.push_back(new TH1F(
         Form("fph1f_pullelec_%d_%d", thetamin, thetamax),
         Form("Electron pull | %d < #theta < %d; Pull; Count", thetamin, thetamax), 100,
         -20, 20));
   }
}

// Destructor
Process::~Process()
{
   // Don't delete pointers to TFile and TTree, ROOT handles them

   // Delete single pointers
   delete fptf1_WFmean;
   delete fptf1_XPmean;
   delete fptf1_WFreso;
   delete fptf1_XPreso;

   delete ptge_absmom_mean_WF;
   delete ptge_absmom_mean_XP;
   delete ptge_absmom_reso_WF;
   delete ptge_absmom_reso_XP;

   delete ptge_mom_mean_WF;
   delete ptge_mom_mean_XP;
   delete ptge_mom_reso_WF;
   delete ptge_mom_reso_XP;

   delete ptge_mom_mean_pullmu;
   delete ptge_mom_mean_pullelec;
   delete ptge_mom_std_pullmu;
   delete ptge_mom_std_pullelec;

   delete ptge_X_mean_WF;
   delete ptge_X_mean_XP;
   delete ptge_X_reso_WF;
   delete ptge_X_reso_XP;

   delete ptge_dd_mean_WF;
   delete ptge_dd_mean_XP;
   delete ptge_dd_reso_WF;
   delete ptge_dd_reso_XP;

   delete ptge_dt_mean_WF;
   delete ptge_dt_mean_XP;
   delete ptge_dt_reso_WF;
   delete ptge_dt_reso_XP;

   delete ptge_trklen_mean_WF;
   delete ptge_trklen_mean_XP;
   delete ptge_trklen_reso_WF;
   delete ptge_trklen_reso_XP;

   delete ptge_absphi_mean_WF;
   delete ptge_absphi_mean_XP;
   delete ptge_absphi_reso_WF;
   delete ptge_absphi_reso_XP;

   delete ptge_phi_mean_WF;
   delete ptge_phi_mean_XP;
   delete ptge_phi_reso_WF;
   delete ptge_phi_reso_XP;

   delete ptge_abstheta_mean_WF;
   delete ptge_abstheta_mean_XP;
   delete ptge_abstheta_reso_WF;
   delete ptge_abstheta_reso_XP;

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
   delete fph2f_WFdrift;
   delete fph2f_XPdrift;
   delete fph2f_WFX;
   delete fph2f_XPX;
   delete fph2f_XPlen;
   delete fph2f_XPphi;
   delete fph2f_XPtheta;
   delete fph2f_XPabsmommean;
   delete fph2f_WFabsmommean;
   delete fph2f_XPmom;
   delete fph2f_WFmom;
   delete fph1i_mom;
   delete fph1i_mom_tHAT;
   delete fph1i_mom_bHAT;
   delete fph2f_XZ;
   delete fph2f_YZ;
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
   delete fph1f_pullmu;
   delete fph1f_pullelec;
   delete fph1f_pullproton;
   delete fph2f_momphi;
   delete fph2f_timeX;
   delete fph2f_pullelecmu;
   delete fph2f_chi2ndfphi;
   delete fph2f_momtheta;
   delete fph2f_momR;
   delete fph2f_chi2ndfR;
   delete fph2f_lentheta;
   delete fph2f_lenphi;
   delete fph1f_dir0;
   delete fph1f_dir1;
   delete fph1f_dir2;
   delete fph2f_momlen;
   delete fph2f_momncl;

   // Clear vector pointers
   for (auto ptr : vmod_fph1f_XP)
      delete ptr;
   for (auto ptr : vmod_fph1f_WF)
      delete ptr;
   for (auto ptr : vmod_fph2f_XPdrift)
      delete ptr;
   for (auto ptr : vmom_fph1f_WF)
      delete ptr;
   for (auto ptr : vabsmommean_fph1f_XP)
      delete ptr;
   for (auto ptr : vabsmommean_fph1f_WF)
      delete ptr;
   for (auto ptr : vabsmomreso_fph1f_XP)
      delete ptr;
   for (auto ptr : vabsmomreso_fph1f_WF)
      delete ptr;
   for (auto ptr : vmom_fph1f_XP)
      delete ptr;
   for (auto ptr : vmom_fph1f_pullmu)
      delete ptr;
   for (auto ptr : vmom_fph1f_pullelec)
      delete ptr;
   for (auto ptr : vX_fph1f_WF)
      delete ptr;
   for (auto ptr : vX_fph1f_XP)
      delete ptr;
   for (auto ptr : vddmean_fph1f_WF)
      delete ptr;
   for (auto ptr : vddmean_fph1f_XP)
      delete ptr;
   for (auto ptr : vddreso_fph1f_WF)
      delete ptr;
   for (auto ptr : vddreso_fph1f_XP)
      delete ptr;
   for (auto ptr : vdt_fph1f_WF)
      delete ptr;
   for (auto ptr : vdt_fph1f_XP)
      delete ptr;
   for (auto ptr : vtrklenmean_fph1f_WF)
      delete ptr;
   for (auto ptr : vtrklenmean_fph1f_XP)
      delete ptr;
   for (auto ptr : vtrklenreso_fph1f_WF)
      delete ptr;
   for (auto ptr : vtrklenreso_fph1f_XP)
      delete ptr;
   for (auto ptr : vabsphimean_fph1f_WF)
      delete ptr;
   for (auto ptr : vabsphimean_fph1f_XP)
      delete ptr;
   for (auto ptr : vabsphireso_fph1f_WF)
      delete ptr;
   for (auto ptr : vabsphireso_fph1f_XP)
      delete ptr;
   for (auto ptr : vphi_fph1f_WF)
      delete ptr;
   for (auto ptr : vphi_fph1f_XP)
      delete ptr;
   for (auto ptr : vtheta_fph1f_WF)
      delete ptr;
   for (auto ptr : vtheta_fph1f_XP)
      delete ptr;
   for (auto ptr : vabsthetamean_fph1f_WF)
      delete ptr;
   for (auto ptr : vabsthetamean_fph1f_XP)
      delete ptr;
   for (auto ptr : vabsthetareso_fph1f_WF)
      delete ptr;
   for (auto ptr : vabsthetareso_fph1f_XP)
      delete ptr;
   for (auto ptr : vtheta_fph1f_pullmu)
      delete ptr;
   for (auto ptr : vtheta_fph1f_pullelec)
      delete ptr;

   // Clear vectors after deleting
   vmod_fph1f_XP.clear();
   vmod_fph1f_WF.clear();
   vmod_fph2f_XPdrift.clear();
   vabsmommean_fph1f_WF.clear();
   vabsmommean_fph1f_XP.clear();
   vabsmomreso_fph1f_WF.clear();
   vabsmomreso_fph1f_XP.clear();
   vmom_fph1f_WF.clear();
   vmom_fph1f_XP.clear();
   vmom_fph1f_pullmu.clear();
   vmom_fph1f_pullelec.clear();
   vX_fph1f_WF.clear();
   vX_fph1f_XP.clear();
   vdt_fph1f_WF.clear();
   vdt_fph1f_XP.clear();
   vddmean_fph1f_WF.clear();
   vddmean_fph1f_XP.clear();
   vddreso_fph1f_WF.clear();
   vddreso_fph1f_XP.clear();
   vtrklenmean_fph1f_WF.clear();
   vtrklenmean_fph1f_XP.clear();
   vtrklenreso_fph1f_WF.clear();
   vtrklenreso_fph1f_XP.clear();
   vabsphimean_fph1f_WF.clear();
   vabsphimean_fph1f_XP.clear();
   vabsphireso_fph1f_WF.clear();
   vabsphireso_fph1f_XP.clear();
   vphi_fph1f_WF.clear();
   vphi_fph1f_XP.clear();
   vabsthetamean_fph1f_WF.clear();
   vabsthetamean_fph1f_XP.clear();
   vabsthetareso_fph1f_WF.clear();
   vabsthetareso_fph1f_XP.clear();
   vtheta_fph1f_WF.clear();
   vtheta_fph1f_XP.clear();
   vtheta_fph1f_pullmu.clear();
   vtheta_fph1f_pullelec.clear();
}

void Process::SetCuts()
{
   int nclmincosmics = 0, nclminbeam = 0;
   if (ffileName.find("dog1") != std::string::npos or
       ffileName.find("cosmics") != std::string::npos) {
      nclmin = 60;
      nclmincosmics = 60;
   } else if (ffileName.find("beam") != std::string::npos or
              ffileName.find("sandmu") != std::string::npos) {
      nclmin = 130;
      nclminbeam = 130;
   }
   fcutslist +=
      ("_B" + std::to_string(nclminbeam) + "C" + std::to_string(nclmincosmics) + "ncl");

   // nclmin = 50, nclmax = 150;
   // fcutslist += ("_" + std::to_string(nclmin) + "ncl" + std::to_string(nclmax));

   // tmeanmin = 0, tmeanmax = 75;
   // fcutslist += ("_" + std::to_string(tmeanmin) + "tmin" + std::to_string(tmeanmax));

   // dxmin = 50, dxmax = 150;
   // fcutslist += ("_" + std::to_string(dxmin) + "dx" + std::to_string(dxmax));

   // pullmumax = 2;
   // fcutslist += ("_pullmu" + std::to_string(pullmumax));

   // chi2max = 5;
   // fcutslist += ("_chi2ndf" + std::to_string(chi2max));

   // nclmin = 130;
   // fcutslist += ("_" + std::to_string(nclmin) + "ncl");

   momcutlow = 300, momcuthigh = 500;
   fcutslist += ("_" + std::to_string(momcutlow) + "mom" + std::to_string(momcuthigh));

   // xcutmin = -981, xcutmax = 981;
   // fcutslist += ("_" + std::to_string(xcutmin) + "x" + std::to_string(xcutmax));

   // phimin = 0;
   // phimax = 60;
   // fcutslist += ("_" + std::to_string(phimin) + "phi" + std::to_string(phimax));

   // thetamin = 0;
   // thetamax = 5;
   // fcutslist += ("_" + std::to_string(thetamin) + "theta" + std::to_string(thetamax));

   // momcutlow = 100;
   // fcutslist += ("_" + std::to_string(momcutlow) + "mom");

   // if (ffileName.find("dog1") != std::string::npos or
   //     ffileName.find("cosmics") != std::string::npos){
   //    tstartmin = 100, tstartmax = 1000;
   //    fcutslist += ("_" + std::to_string(tstartmin) + "tstart" +
   //    std::to_string(tstartmax));
   // }

   // fcutslist += "_dir1>0flip";
   // fcutslist += "flipZ";

   // hat = 1;
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
   fpInputTree->SetBranchAddress("T0type", &T0type);

   // Debugging variables
   int dirYnegbHAT = 0, dirYposbHAT = 0, dirYnegtHAT = 0, dirYpostHAT = 0;
   int nbotneg = 0, nbotpos = 0, ntopneg = 0, ntoppos = 0;

   for (int i = 0; i < nentries; i++) {
      fpInputTree->GetEntry(i);

      // Flip momentum sign if necessary ----------------------------------------------
      mom = mom_og;
      // FOR COSMICS if the track is going up, the sign must be flipped
      if (ffileName.find("dog1") != std::string::npos or
          ffileName.find("cosmic") != std::string::npos)
         dir[1] > 0 ? mom = -mom_og : mom = mom_og;
      // For BEAM if the track is going downstream, the sign must be flipped
      else if (ffileName.find("beam") != std::string::npos or
               ffileName.find("sandmu") != std::string::npos) {
         if (dir[2] < 0)
            mom = -mom_og;
      }

      // Define angles ---------------------------------------------------------------
      float dirxy = TMath::Sqrt(dir[0] * dir[0] + dir[1] * dir[1]);
      float diryz = TMath::Sqrt(dir[1] * dir[1] + dir[2] * dir[2]);
      phi = TMath::ATan(dir[1] / dir[2]) * 180 / TMath::Pi();
      theta = TMath::ATan(dir[0] / diryz) * 180 / TMath::Pi();

      // Cuts --------------------------------------------------------------------------
      if (ffileName.find("beam") != std::string::npos or
          ffileName.find("sandmu") != std::string::npos) {
         if (T0type == 0 or T0type > 24)
            continue;
      }
      if (std::isnan(mom))
         continue;
      if (fabs(phi) < phimin or phimax < fabs(phi))
         continue;
      if (fabs(theta) < thetamin or thetamax < fabs(theta))
         continue;
      if (mean_time < tmeanmin or mean_time > tmeanmax)
         continue;
      if (start_time < tstartmin or start_time > tstartmax)
         continue;
      if (pos[0] < xcutmin or pos[0] > xcutmax)
         continue;
      if (chi2 / NDF > chi2max)
         continue;
      if (hat != 0 and sign(hat) != sign(pos[1]))
         continue;
      // if(fabs(pull_muon) > fabs(pull_ele)) continue;
      if (fabs(pull_muon) > pullmumax)
         continue;

      // Momentum related histograms ----------------------------------------------------
      // apply locally length cut
      if (nclmin < ncl and ncl < nclmax and dxmin < dx / 10 and dx / 10 < dxmax) {
         fph1i_mom->Fill(mom);
         if (pos[1] < 0)
            fph1i_mom_bHAT->Fill(mom);
         if (pos[1] > 0)
            fph1i_mom_tHAT->Fill(mom);
         fph2f_WFabsmommean->Fill(fabs(mom), wf);
         fph2f_XPabsmommean->Fill(fabs(mom), xp);
         fph2f_WFmom->Fill(mom, wf);
         fph2f_XPmom->Fill(mom, xp);

         // dEdx vs absolute momentum
         absmomresoindex = (int)std::floor(fabs(mom) / absmomresobinwidth);
         absmommeanindex = (int)std::floor(fabs(mom) / absmommeanbinwidth);
         if (fabs(mom) < absmomrange) {
            vabsmommean_fph1f_WF[absmommeanindex]->Fill(wf);
            vabsmommean_fph1f_XP[absmommeanindex]->Fill(xp);
            vabsmomreso_fph1f_WF[absmomresoindex]->Fill(wf);
            vabsmomreso_fph1f_XP[absmomresoindex]->Fill(xp);
         }

         // dEdx vs momentum
         momindex = (int)std::round(mom / mombinwidth) + nmombins / 2;
         if (fabs(mom) < momrange) {
            vmom_fph1f_WF[momindex]->Fill(wf);
            vmom_fph1f_XP[momindex]->Fill(xp);
            vmom_fph1f_pullmu[momindex]->Fill(pull_muon);
            vmom_fph1f_pullelec[momindex]->Fill(pull_ele);
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
         fph2f_chi2mom->Fill(mom, chi2 / NDF);
      }

      // Cut on momentum range -----------------------------------------------------------
      if (fabs(mom) < momcutlow or fabs(mom) > momcuthigh)
         continue;

      // Track length related histograms ----------------------------------------------
      // dEdx vs track length
      trklenmeanindex = (int)std::floor(dx / trklenmeanbinwidth);
      trklenresoindex = (int)std::floor(dx / trklenresobinwidth);
      if (trklenmeanindex < ntrklenmeanbins) {
         vtrklenmean_fph1f_WF[trklenmeanindex]->Fill(wf);
         vtrklenmean_fph1f_XP[trklenmeanindex]->Fill(xp);
         vtrklenreso_fph1f_WF[trklenresoindex]->Fill(wf);
         vtrklenreso_fph1f_XP[trklenresoindex]->Fill(xp);
      }
      // Other track length histograms
      fph2f_XPlen->Fill(dx / 10, xp);
      fph1f_trklen->Fill(dx / 10);
      fph2f_momlen->Fill(dx / 10, mom);
      fph2f_momncl->Fill(ncl, mom);
      // Debug
      fph2f_lentheta->Fill(theta, dx / 10);
      fph2f_lenphi->Fill(phi, dx / 10);

      // Cut on track length -----------------------------------------------------------
      if (nclmin > ncl or ncl > nclmax)
         continue;
      if (dxmin > dx / 10 or dxmax < dx / 10)
         continue;

      // Rest of the histograms --------------------------------------------------------
      fph1f_WF->Fill(wf);
      fph1f_XP->Fill(xp);
      fph2f_WFXP->Fill(wf, xp);
      fph2f_WFdrift->Fill(mean_time, wf);
      fph2f_XPdrift->Fill(mean_time, xp);
      fph2f_WFX->Fill(pos[0], wf);
      fph2f_XPX->Fill(pos[0], xp);
      fph2f_XPphi->Fill(phi, xp);
      fph2f_XPtheta->Fill(theta, xp);

      vmod_fph2f_XPdrift[eram_channel]->Fill(mean_time, xp);
      vmod_fph1f_WF[eram_channel]->Fill(wf);
      vmod_fph1f_XP[eram_channel]->Fill(xp);

      fph2f_XZ->Fill(dir[0], dir[2]);
      fph2f_YZ->Fill(dir[1], dir[2]);
      fph1f_phi->Fill(phi);
      fph1f_theta->Fill(theta);
      fph2f_phitheta->Fill(theta, phi);
      fph1f_chi2->Fill(chi2 / NDF);

      // dEdx vs X position
      xposindex = (int)std::round(pos[0] / xposbinwidth) + nxposbins / 2;
      if (xposindex < nxposbins) {
         vX_fph1f_WF[xposindex]->Fill(wf);
         vX_fph1f_XP[xposindex]->Fill(xp);
      }

      // dEdx vs drift distance
      ddmeanindex = (int)std::floor((981 - fabs(pos[0])) / ddmeanbinwidth);
      ddresoindex = (int)std::floor((981 - fabs(pos[0])) / ddresobinwidth);
      if (ddmeanindex < nddmeanbins and ddmeanindex >= 0) {
         if (fabs(pos[0]) > 0) {
            vddmean_fph1f_WF[ddmeanindex]->Fill(wf);
            vddmean_fph1f_XP[ddmeanindex]->Fill(xp);
            vddreso_fph1f_WF[ddresoindex]->Fill(wf);
            vddreso_fph1f_XP[ddresoindex]->Fill(xp);
         }
      }

      // dEdx vs drift time
      dtindex = (int)std::floor(mean_time / dtbinwidth);
      if (dtindex < ndtbins) {
         vdt_fph1f_WF[dtindex]->Fill(wf);
         vdt_fph1f_XP[dtindex]->Fill(xp);
      }

      // dEdx vs absolute phi angle
      absphimeanindex = (int)std::floor(fabs(phi) / absphimeanbinwidth);
      absphiresoindex = (int)std::floor(fabs(phi) / absphiresobinwidth);
      if (fabs(phi) < absphirange) {
         vabsphimean_fph1f_WF[absphimeanindex]->Fill(wf);
         vabsphimean_fph1f_XP[absphimeanindex]->Fill(xp);
         vabsphireso_fph1f_WF[absphiresoindex]->Fill(wf);
         vabsphireso_fph1f_XP[absphiresoindex]->Fill(xp);
      }

      // dEdx vs phi angle
      phiindex = (int)std::round(phi / phibinwidth) + nphibins / 2;
      if (phiindex < nphibins) {
         vphi_fph1f_WF[phiindex]->Fill(wf);
         vphi_fph1f_XP[phiindex]->Fill(xp);
      }

      // dEdx vs absolute theta angle
      absthetameanindex = (int)std::floor(fabs(theta) / absthetameanbinwidth);
      absthetaresoindex = (int)std::floor(fabs(theta) / absthetaresobinwidth);
      if (fabs(theta) < absthetarange) {
         vabsthetamean_fph1f_WF[absthetameanindex]->Fill(wf);
         vabsthetamean_fph1f_XP[absthetameanindex]->Fill(xp);
         vabsthetareso_fph1f_WF[absthetaresoindex]->Fill(wf);
         vabsthetareso_fph1f_XP[absthetaresoindex]->Fill(xp);
      }

      // dEdx vs theta angle
      thetaindex = (int)std::round(theta / thetabinwidth) + nthetabins / 2;
      if (thetaindex < nthetabins) {
         vtheta_fph1f_WF[thetaindex]->Fill(wf);
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
      fph2f_chi2ndfR->Fill(1 / curv, chi2 / NDF);
      fph1f_pullmu->Fill(pull_muon);
      fph1f_pullelec->Fill(pull_ele);
      fph1f_pullproton->Fill(pull_proton);
      fph2f_timeX->Fill(pos[0], mean_time);
      fph1f_dir0->Fill(dir[0]);
      fph1f_dir1->Fill(dir[1]);
      fph1f_dir2->Fill(dir[2]);
   }

   // TGraph filling
   // absolute momentum (mean) ----------------------------------------------------------
   int ivalid = 0;
   for (int i = 0; i < nabsmommeanbins; i++) {
      int nentries_here = vabsmommean_fph1f_WF[i]->GetEntries();
      if (nentries_here < 50)
         continue;
      fptf1_WFmean = Fit1Gauss(vabsmommean_fph1f_WF[i]);
      fptf1_XPmean = Fit1Gauss(vabsmommean_fph1f_XP[i]);

      float mean_WF = fptf1_WFmean->GetParameter(1);
      float mean_XP = fptf1_XPmean->GetParameter(1);
      float dmean_WF = fptf1_WFmean->GetParError(1);
      float dmean_XP = fptf1_XPmean->GetParError(1);

      if (mean_WF == 0 || mean_XP == 0)
         continue;

      ptge_absmom_mean_WF->SetPoint(ivalid, (2 * i + 1.) / 2 * absmommeanbinwidth,
                                    mean_WF);
      ptge_absmom_mean_WF->SetPointError(ivalid, absmommeanbinwidth / 2, dmean_WF);
      ptge_absmom_mean_XP->SetPoint(ivalid, (2 * i + 1.) / 2 * absmommeanbinwidth,
                                    mean_XP);
      ptge_absmom_mean_XP->SetPointError(ivalid, absmommeanbinwidth / 2, dmean_XP);
      ivalid++;
   }

   // absolute momentum (resolution) ----------------------------------------------------
   ivalid = 0;
   for (int i = 0; i < nabsmomresobins; i++) {
      int nentries_here = vabsmomreso_fph1f_WF[i]->GetEntries();
      if (nentries_here < 50)
         continue;
      fptf1_WFreso = Fit1Gauss(vabsmomreso_fph1f_WF[i]);
      fptf1_XPreso = Fit1Gauss(vabsmomreso_fph1f_XP[i]);

      float mean_WF = fptf1_WFreso->GetParameter(1);
      float mean_XP = fptf1_XPreso->GetParameter(1);
      float dmean_WF = fptf1_WFreso->GetParError(1);
      float dmean_XP = fptf1_XPreso->GetParError(1);

      float std_WF = fptf1_WFreso->GetParameter(2);
      float std_XP = fptf1_XPreso->GetParameter(2);

      float reso_WF = std_WF / mean_WF * 100;
      float reso_XP = std_XP / mean_XP * 100;
      float dreso_WF = GetResoError(fptf1_WFreso);
      float dreso_XP = GetResoError(fptf1_XPreso);

      if (mean_WF == 0 || mean_XP == 0)
         continue;

      ptge_absmom_reso_WF->SetPoint(ivalid, (2 * i + 1.) / 2 * absmomresobinwidth,
                                    reso_WF);
      ptge_absmom_reso_WF->SetPointError(ivalid, absmomresobinwidth / 2, dreso_WF);
      ptge_absmom_reso_XP->SetPoint(ivalid, (2 * i + 1.) / 2 * absmomresobinwidth,
                                    reso_XP);
      ptge_absmom_reso_XP->SetPointError(ivalid, absmomresobinwidth / 2, dreso_XP);
      ivalid++;
   }

   // Momentum
   ivalid = 0;
   for (int i = 0; i < nmombins; i++) {
      int nentries_here = vmom_fph1f_WF[i]->GetEntries();
      if (nentries_here < 50)
         continue;
      fptf1_WFmean = Fit1Gauss(vmom_fph1f_WF[i]);
      fptf1_XPmean = Fit1Gauss(vmom_fph1f_XP[i]);
      fptf1_pullmu = Fit1Gauss(vmom_fph1f_pullmu[i]);
      fptf1_pullelec = Fit1Gauss(vmom_fph1f_pullelec[i]);

      float mean_WF = fptf1_WFmean->GetParameter(1);
      float mean_XP = fptf1_XPmean->GetParameter(1);
      float mean_pullmu = fptf1_pullmu->GetParameter(1);
      float mean_pullelec = fptf1_pullelec->GetParameter(1);
      float dmean_WF = fptf1_WFmean->GetParError(1);
      float dmean_XP = fptf1_XPmean->GetParError(1);
      float dmean_pullmu = fptf1_pullmu->GetParError(1);
      float dmean_pullelec = fptf1_pullelec->GetParError(1);

      float std_WF = fptf1_WFmean->GetParameter(2);
      float std_XP = fptf1_XPmean->GetParameter(2);
      float std_pullmu = fptf1_pullmu->GetParameter(2);
      float std_pullelec = fptf1_pullelec->GetParameter(2);
      float dstd_pullmu = fptf1_pullmu->GetParError(2);
      float dstd_pullelec = fptf1_pullelec->GetParError(2);

      float reso_WF = std_WF / mean_WF * 100;
      float reso_XP = std_XP / mean_XP * 100;
      float dreso_WF = GetResoError(fptf1_WFmean);
      float dreso_XP = GetResoError(fptf1_XPmean);

      if (mean_WF == 0 || mean_XP == 0)
         continue;

      ptge_mom_mean_WF->SetPoint(ivalid, i * mombinwidth - momrange, mean_WF);
      ptge_mom_mean_WF->SetPointError(ivalid, mombinwidth / 2, dmean_WF);
      ptge_mom_mean_XP->SetPoint(ivalid, i * mombinwidth - momrange, mean_XP);
      ptge_mom_mean_XP->SetPointError(ivalid, mombinwidth / 2, dmean_XP);
      ptge_mom_mean_pullmu->SetPoint(ivalid, i * mombinwidth - momrange, mean_pullmu);
      ptge_mom_mean_pullmu->SetPointError(ivalid, mombinwidth / 2, dmean_pullmu);
      ptge_mom_mean_pullelec->SetPoint(ivalid, i * mombinwidth - momrange, mean_pullelec);
      ptge_mom_mean_pullelec->SetPointError(ivalid, mombinwidth / 2, dmean_pullelec);

      ptge_mom_reso_WF->SetPoint(ivalid, i * mombinwidth - momrange, reso_WF);
      ptge_mom_reso_WF->SetPointError(ivalid, mombinwidth / 2, dreso_WF);
      ptge_mom_reso_XP->SetPoint(ivalid, i * mombinwidth - momrange, reso_XP);
      ptge_mom_reso_XP->SetPointError(ivalid, mombinwidth / 2, dreso_XP);
      ptge_mom_std_pullmu->SetPoint(ivalid, i * mombinwidth - momrange, std_pullmu);
      ptge_mom_std_pullmu->SetPointError(ivalid, mombinwidth / 2, dstd_pullmu);
      ptge_mom_std_pullelec->SetPoint(ivalid, i * mombinwidth - momrange, std_pullelec);
      ptge_mom_std_pullelec->SetPointError(ivalid, mombinwidth / 2, dstd_pullelec);

      ivalid++;
   }

   // X postion
   ivalid = 0;
   for (int i = 0; i < nxposbins; i++) {
      int nentries_here = vX_fph1f_WF[i]->GetEntries();
      if (nentries_here < 50)
         continue;
      fptf1_WFmean = Fit1Gauss(vX_fph1f_WF[i]);
      fptf1_XPmean = Fit1Gauss(vX_fph1f_XP[i]);

      float mean_WF = fptf1_WFmean->GetParameter(1);
      float mean_XP = fptf1_XPmean->GetParameter(1);
      float dmean_WF = fptf1_WFmean->GetParError(1);
      float dmean_XP = fptf1_XPmean->GetParError(1);

      float std_WF = fptf1_WFmean->GetParameter(2);
      float std_XP = fptf1_XPmean->GetParameter(2);

      float reso_WF = std_WF / mean_WF * 100;
      float reso_XP = std_XP / mean_XP * 100;
      float dreso_WF = GetResoError(fptf1_WFmean);
      float dreso_XP = GetResoError(fptf1_XPmean);

      if (mean_WF == 0 || mean_XP == 0)
         continue;

      ptge_X_mean_WF->SetPoint(ivalid, i * xposbinwidth - xposrange, mean_WF);
      ptge_X_mean_WF->SetPointError(ivalid, xposbinwidth / 2, dmean_WF);
      ptge_X_mean_XP->SetPoint(ivalid, i * xposbinwidth - xposrange, mean_XP);
      ptge_X_mean_XP->SetPointError(ivalid, xposbinwidth / 2, dmean_XP);

      ptge_X_reso_WF->SetPoint(ivalid, i * xposbinwidth - xposrange, reso_WF);
      ptge_X_reso_WF->SetPointError(ivalid, xposbinwidth / 2, dreso_WF);
      ptge_X_reso_XP->SetPoint(ivalid, i * xposbinwidth - xposrange, reso_XP);
      ptge_X_reso_XP->SetPointError(ivalid, xposbinwidth / 2, dreso_XP);
      ivalid++;
   }

   // Drift distance (mean)
   ivalid = 0;
   for (int i = 0; i < nddmeanbins; i++) {
      int nentries_here = vddmean_fph1f_WF[i]->GetEntries();
      if (nentries_here < 50)
         continue;
      fptf1_WFmean = Fit1Gauss(vddmean_fph1f_WF[i]);
      fptf1_XPmean = Fit1Gauss(vddmean_fph1f_XP[i]);

      float mean_WF = fptf1_WFmean->GetParameter(1);
      float mean_XP = fptf1_XPmean->GetParameter(1);
      float dmean_WF = fptf1_WFmean->GetParError(1);
      float dmean_XP = fptf1_XPmean->GetParError(1);

      if (mean_WF == 0 || mean_XP == 0)
         continue;

      ptge_dd_mean_WF->SetPoint(ivalid, (2 * i + 1.) / 2 * ddmeanbinwidth, mean_WF);
      ptge_dd_mean_WF->SetPointError(ivalid, ddmeanbinwidth / 2, dmean_WF);
      ptge_dd_mean_XP->SetPoint(ivalid, (2 * i + 1.) / 2 * ddmeanbinwidth, mean_XP);
      ptge_dd_mean_XP->SetPointError(ivalid, ddmeanbinwidth / 2, dmean_XP);
      ivalid++;
   }

   // Drift distance (resolution)
   ivalid = 0;
   for (int i = 0; i < nddresobins; i++) {
      int nentries_here = vddreso_fph1f_WF[i]->GetEntries();
      if (nentries_here < 50)
         continue;
      fptf1_WFreso = Fit1Gauss(vddreso_fph1f_WF[i]);
      fptf1_XPreso = Fit1Gauss(vddreso_fph1f_XP[i]);

      float mean_WF = fptf1_WFreso->GetParameter(1);
      float mean_XP = fptf1_XPreso->GetParameter(1);
      float dmean_WF = fptf1_WFreso->GetParError(1);
      float dmean_XP = fptf1_XPreso->GetParError(1);

      float std_WF = fptf1_WFreso->GetParameter(2);
      float std_XP = fptf1_XPreso->GetParameter(2);

      float reso_WF = std_WF / mean_WF * 100;
      float reso_XP = std_XP / mean_XP * 100;
      float dreso_WF = GetResoError(fptf1_WFreso);
      float dreso_XP = GetResoError(fptf1_XPreso);

      if (mean_WF == 0 || mean_XP == 0)
         continue;

      ptge_dd_reso_WF->SetPoint(ivalid, (2 * i + 1.) / 2 * ddresobinwidth, reso_WF);
      ptge_dd_reso_WF->SetPointError(ivalid, ddresobinwidth / 2, dreso_WF);
      ptge_dd_reso_XP->SetPoint(ivalid, (2 * i + 1.) / 2 * ddresobinwidth, reso_XP);
      ptge_dd_reso_XP->SetPointError(ivalid, ddresobinwidth / 2, dreso_XP);
      ivalid++;
   }

   // Drift time
   ivalid = 0;
   for (int i = 0; i < ndtbins; i++) {
      int nentries_here = vdt_fph1f_WF[i]->GetEntries();
      if (nentries_here < 50)
         continue;
      fptf1_WFmean = Fit1Gauss(vdt_fph1f_WF[i]);
      fptf1_XPmean = Fit1Gauss(vdt_fph1f_XP[i]);

      float mean_WF = fptf1_WFmean->GetParameter(1);
      float mean_XP = fptf1_XPmean->GetParameter(1);
      float dmean_WF = fptf1_WFmean->GetParError(1);
      float dmean_XP = fptf1_XPmean->GetParError(1);

      float std_WF = fptf1_WFmean->GetParameter(2);
      float std_XP = fptf1_XPmean->GetParameter(2);

      float reso_WF = std_WF / mean_WF * 100;
      float reso_XP = std_XP / mean_XP * 100;
      float dreso_WF = GetResoError(fptf1_WFmean);
      float dreso_XP = GetResoError(fptf1_XPmean);

      if (mean_WF == 0 || mean_XP == 0)
         continue;

      ptge_dt_mean_WF->SetPoint(ivalid, (2 * i + 1.) / 2 * dtbinwidth, mean_WF);
      ptge_dt_mean_WF->SetPointError(ivalid, dtbinwidth / 2, dmean_WF);
      ptge_dt_mean_XP->SetPoint(ivalid, (2 * i + 1.) / 2 * dtbinwidth, mean_XP);
      ptge_dt_mean_XP->SetPointError(ivalid, dtbinwidth / 2, dmean_XP);

      ptge_dt_reso_WF->SetPoint(ivalid, (2 * i + 1.) / 2 * dtbinwidth, reso_WF);
      ptge_dt_reso_WF->SetPointError(ivalid, dtbinwidth / 2, dreso_WF);
      ptge_dt_reso_XP->SetPoint(ivalid, (2 * i + 1.) / 2 * dtbinwidth, reso_XP);
      ptge_dt_reso_XP->SetPointError(ivalid, dtbinwidth / 2, dreso_XP);
      ivalid++;
   }

   // track length (mean)
   ivalid = 0;
   for (int i = 0; i < ntrklenmeanbins; i++) {
      int nentries_here = vtrklenmean_fph1f_WF[i]->GetEntries();
      if (nentries_here < 50)
         continue;
      fptf1_WFmean = Fit1Gauss(vtrklenmean_fph1f_WF[i]);
      fptf1_XPmean = Fit1Gauss(vtrklenmean_fph1f_XP[i]);

      float mean_WF = fptf1_WFmean->GetParameter(1);
      float mean_XP = fptf1_XPmean->GetParameter(1);
      float dmean_WF = fptf1_WFmean->GetParError(1);
      float dmean_XP = fptf1_XPmean->GetParError(1);

      if (mean_WF == 0 || mean_XP == 0)
         continue;

      ptge_trklen_mean_WF->SetPoint(ivalid, (2 * i + 1.) / 2 * trklenmeanbinwidth,
                                    mean_WF);
      ptge_trklen_mean_WF->SetPointError(ivalid, trklenmeanbinwidth / 2, dmean_WF);
      ptge_trklen_mean_XP->SetPoint(ivalid, (2 * i + 1.) / 2 * trklenmeanbinwidth,
                                    mean_XP);
      ptge_trklen_mean_XP->SetPointError(ivalid, trklenmeanbinwidth / 2, dmean_XP);
      ivalid++;
   }

   // track length (resolution)
   ivalid = 0;
   for (int i = 0; i < ntrklenresobins; i++) {
      int nentries_here = vtrklenreso_fph1f_WF[i]->GetEntries();
      if (nentries_here < 50)
         continue;
      fptf1_WFreso = Fit1Gauss(vtrklenreso_fph1f_WF[i]);
      fptf1_XPreso = Fit1Gauss(vtrklenreso_fph1f_XP[i]);

      float mean_WF = fptf1_WFreso->GetParameter(1);
      float mean_XP = fptf1_XPreso->GetParameter(1);
      float dmean_WF = fptf1_WFreso->GetParError(1);
      float dmean_XP = fptf1_XPreso->GetParError(1);

      float std_WF = fptf1_WFreso->GetParameter(2);
      float std_XP = fptf1_XPreso->GetParameter(2);

      float reso_WF = std_WF / mean_WF * 100;
      float reso_XP = std_XP / mean_XP * 100;
      float dreso_WF = GetResoError(fptf1_WFreso);
      float dreso_XP = GetResoError(fptf1_XPreso);

      if (mean_WF == 0 || mean_XP == 0)
         continue;

      ptge_trklen_reso_WF->SetPoint(ivalid, (2 * i + 1.) / 2 * trklenresobinwidth,
                                    reso_WF);
      ptge_trklen_reso_WF->SetPointError(ivalid, trklenresobinwidth / 2, dreso_WF);
      ptge_trklen_reso_XP->SetPoint(ivalid, (2 * i + 1.) / 2 * trklenresobinwidth,
                                    reso_XP);
      ptge_trklen_reso_XP->SetPointError(ivalid, trklenresobinwidth / 2, dreso_XP);
      ivalid++;
   }

   // absolute phi angle (mean)
   ivalid = 0;
   for (int i = 0; i < nabsphimeanbins; i++) {
      int nentries_here = vabsphimean_fph1f_WF[i]->GetEntries();
      if (nentries_here < 50)
         continue;
      fptf1_WFmean = Fit1Gauss(vabsphimean_fph1f_WF[i]);
      fptf1_XPmean = Fit1Gauss(vabsphimean_fph1f_XP[i]);

      float mean_WF = fptf1_WFmean->GetParameter(1);
      float mean_XP = fptf1_XPmean->GetParameter(1);
      float dmean_WF = fptf1_WFmean->GetParError(1);
      float dmean_XP = fptf1_XPmean->GetParError(1);

      if (mean_WF == 0 || mean_XP == 0)
         continue;
      ptge_absphi_mean_WF->SetPoint(ivalid, (2 * i + 1.) / 2 * absphimeanbinwidth,
                                    mean_WF);
      ptge_absphi_mean_WF->SetPointError(ivalid, absphimeanbinwidth / 2, dmean_WF);
      ptge_absphi_mean_XP->SetPoint(ivalid, (2 * i + 1.) / 2 * absphimeanbinwidth,
                                    mean_XP);
      ptge_absphi_mean_XP->SetPointError(ivalid, absphimeanbinwidth / 2, dmean_XP);
      ivalid++;
   }

   // absolute phi angle (resolution)
   ivalid = 0;
   for (int i = 0; i < nabsphiresobins; i++) {
      int nentries_here = vabsphireso_fph1f_WF[i]->GetEntries();
      if (nentries_here < 50)
         continue;
      fptf1_WFreso = Fit1Gauss(vabsphireso_fph1f_WF[i]);
      fptf1_XPreso = Fit1Gauss(vabsphireso_fph1f_XP[i]);

      float mean_WF = fptf1_WFreso->GetParameter(1);
      float mean_XP = fptf1_XPreso->GetParameter(1);
      float dmean_WF = fptf1_WFreso->GetParError(1);
      float dmean_XP = fptf1_XPreso->GetParError(1);

      float std_WF = fptf1_WFreso->GetParameter(2);
      float std_XP = fptf1_XPreso->GetParameter(2);

      float reso_XP = std_XP / mean_XP * 100;
      float reso_WF = std_WF / mean_WF * 100;
      float dreso_WF = GetResoError(fptf1_WFreso);
      float dreso_XP = GetResoError(fptf1_XPreso);

      if (mean_WF == 0 || mean_XP == 0)
         continue;

      ptge_absphi_reso_WF->SetPoint(ivalid, (2 * i + 1.) / 2 * absphiresobinwidth,
                                    reso_WF);
      ptge_absphi_reso_WF->SetPointError(ivalid, absphiresobinwidth / 2, dreso_WF);
      ptge_absphi_reso_XP->SetPoint(ivalid, (2 * i + 1.) / 2 * absphiresobinwidth,
                                    reso_XP);
      ptge_absphi_reso_XP->SetPointError(ivalid, absphiresobinwidth / 2, dreso_XP);
      ivalid++;
   }

   // Phi angle
   ivalid = 0;
   for (int i = 0; i < nphibins; i++) {
      int nentries_here = vphi_fph1f_WF[i]->GetEntries();
      if (nentries_here < 50)
         continue;
      fptf1_WFmean = Fit1Gauss(vphi_fph1f_WF[i]);
      fptf1_XPmean = Fit1Gauss(vphi_fph1f_XP[i]);

      float mean_WF = fptf1_WFmean->GetParameter(1);
      float mean_XP = fptf1_XPmean->GetParameter(1);
      float dmean_WF = fptf1_WFmean->GetParError(1);
      float dmean_XP = fptf1_XPmean->GetParError(1);

      float std_WF = fptf1_WFmean->GetParameter(2);
      float std_XP = fptf1_XPmean->GetParameter(2);

      float reso_WF = std_WF / mean_WF * 100;
      float reso_XP = std_XP / mean_XP * 100;
      float dreso_WF = GetResoError(fptf1_WFmean);
      float dreso_XP = GetResoError(fptf1_XPmean);

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

   // absolute theta angle (mean)
   ivalid = 0;
   for (int i = 0; i < nabsthetameanbins; i++) {
      int nentries_here = vabsthetamean_fph1f_WF[i]->GetEntries();
      if (nentries_here < 50)
         continue;
      fptf1_WFmean = Fit1Gauss(vabsthetamean_fph1f_WF[i]);
      fptf1_XPmean = Fit1Gauss(vabsthetamean_fph1f_XP[i]);

      float mean_WF = fptf1_WFmean->GetParameter(1);
      float mean_XP = fptf1_XPmean->GetParameter(1);
      float dmean_WF = fptf1_WFmean->GetParError(1);
      float dmean_XP = fptf1_XPmean->GetParError(1);

      if (mean_WF == 0 || mean_XP == 0)
         continue;

      ptge_abstheta_mean_WF->SetPoint(ivalid, (2 * i + 1.) / 2 * absthetameanbinwidth,
                                      mean_WF);
      ptge_abstheta_mean_WF->SetPointError(ivalid, absthetameanbinwidth / 2, dmean_WF);
      ptge_abstheta_mean_XP->SetPoint(ivalid, (2 * i + 1.) / 2 * absthetameanbinwidth,
                                      mean_XP);
      ptge_abstheta_mean_XP->SetPointError(ivalid, absthetameanbinwidth / 2, dmean_XP);
      ivalid++;
   }

   // absolute theta angle (resolution)
   ivalid = 0;
   for (int i = 0; i < nabsthetaresobins; i++) {
      int nentries_here = vabsthetareso_fph1f_WF[i]->GetEntries();
      if (nentries_here < 50)
         continue;
      fptf1_WFreso = Fit1Gauss(vabsthetareso_fph1f_WF[i]);
      fptf1_XPreso = Fit1Gauss(vabsthetareso_fph1f_XP[i]);

      float mean_WF = fptf1_WFreso->GetParameter(1);
      float mean_XP = fptf1_XPreso->GetParameter(1);
      float dmean_WF = fptf1_WFreso->GetParError(1);
      float dmean_XP = fptf1_XPreso->GetParError(1);

      float std_WF = fptf1_WFreso->GetParameter(2);
      float std_XP = fptf1_XPreso->GetParameter(2);

      float reso_XP = std_XP / mean_XP * 100;
      float reso_WF = std_WF / mean_WF * 100;
      float dreso_WF = GetResoError(fptf1_WFreso);
      float dreso_XP = GetResoError(fptf1_XPreso);

      if (mean_WF == 0 || mean_XP == 0)
         continue;

      ptge_abstheta_reso_WF->SetPoint(ivalid, (2 * i + 1.) / 2 * absthetaresobinwidth,
                                      reso_WF);
      ptge_abstheta_reso_WF->SetPointError(ivalid, absthetaresobinwidth / 2, dreso_WF);
      ptge_abstheta_reso_XP->SetPoint(ivalid, (2 * i + 1.) / 2 * absthetaresobinwidth,
                                      reso_XP);
      ptge_abstheta_reso_XP->SetPointError(ivalid, absthetaresobinwidth / 2, dreso_XP);
      ivalid++;
   }

   // Theta angle
   ivalid = 0;
   for (int i = 0; i < nthetabins; i++) {
      int nentries_here = vtheta_fph1f_WF[i]->GetEntries();
      if (nentries_here < 50)
         continue;
      fptf1_WFmean = Fit1Gauss(vtheta_fph1f_WF[i]);
      fptf1_XPmean = Fit1Gauss(vtheta_fph1f_XP[i]);
      fptf1_pullmu = Fit1Gauss(vtheta_fph1f_pullmu[i]);
      fptf1_pullelec = Fit1Gauss(vtheta_fph1f_pullelec[i]);

      float mean_WF = fptf1_WFmean->GetParameter(1);
      float mean_XP = fptf1_XPmean->GetParameter(1);
      float mean_pullmu = fptf1_pullmu->GetParameter(1);
      float mean_pullelec = fptf1_pullelec->GetParameter(1);
      float dmean_WF = fptf1_WFmean->GetParError(1);
      float dmean_XP = fptf1_XPmean->GetParError(1);
      float dmean_pullmu = fptf1_pullmu->GetParError(1);
      float dmean_pullelec = fptf1_pullelec->GetParError(1);

      float std_WF = fptf1_WFmean->GetParameter(2);
      float std_XP = fptf1_XPmean->GetParameter(2);
      float std_pullmu = fptf1_pullmu->GetParameter(2);
      float std_pullelec = fptf1_pullelec->GetParameter(2);
      float dstd_pullmu = fptf1_pullmu->GetParError(2);
      float dstd_pullelec = fptf1_pullelec->GetParError(2);

      float reso_WF = std_WF / mean_WF * 100;
      float reso_XP = std_XP / mean_XP * 100;
      float dreso_WF = GetResoError(fptf1_WFmean);
      float dreso_XP = GetResoError(fptf1_XPmean);

      if (mean_WF == 0 || mean_XP == 0)
         continue;

      ptge_theta_mean_WF->SetPoint(ivalid, i * thetabinwidth - thetarange, mean_WF);
      ptge_theta_mean_WF->SetPointError(ivalid, thetabinwidth / 2, dmean_WF);
      ptge_theta_mean_XP->SetPoint(ivalid, i * thetabinwidth - thetarange, mean_XP);
      ptge_theta_mean_XP->SetPointError(ivalid, thetabinwidth / 2, dmean_XP);
      ptge_theta_mean_pullmu->SetPoint(ivalid, i * thetabinwidth - thetarange,
                                       mean_pullmu);
      ptge_theta_mean_pullmu->SetPointError(ivalid, thetabinwidth / 2, dmean_pullmu);
      ptge_theta_mean_pullelec->SetPoint(ivalid, i * thetabinwidth - thetarange,
                                         mean_pullelec);
      ptge_theta_mean_pullelec->SetPointError(ivalid, thetabinwidth / 2, dmean_pullelec);

      ptge_theta_reso_WF->SetPoint(ivalid, i * thetabinwidth - thetarange, reso_WF);
      ptge_theta_reso_WF->SetPointError(ivalid, thetabinwidth / 2, dreso_WF);
      ptge_theta_reso_XP->SetPoint(ivalid, i * thetabinwidth - thetarange, reso_XP);
      ptge_theta_reso_XP->SetPointError(ivalid, thetabinwidth / 2, dreso_XP);
      ptge_theta_std_pullmu->SetPoint(ivalid, i * thetabinwidth - thetarange, std_pullmu);
      ptge_theta_std_pullmu->SetPointError(ivalid, thetabinwidth / 2, dstd_pullmu);
      ptge_theta_std_pullelec->SetPoint(ivalid, i * thetabinwidth - thetarange,
                                        std_pullelec);
      ptge_theta_std_pullelec->SetPointError(ivalid, thetabinwidth / 2, dstd_pullelec);
      ivalid++;
   }
   fpInputFile->Close();
}
