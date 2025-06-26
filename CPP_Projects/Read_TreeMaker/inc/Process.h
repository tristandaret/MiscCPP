#ifndef PROCESS_H
#define PROCESS_H

#include "Util.h"
#include <TObject.h>

class Process : public TObject {
public:
   ClassDef(Process, 1);
   // Constructor
   Process();
   // Destructor
   ~Process();
   // Member functions
   void SetCuts();
   void Run();

   void SetType(const std::string &type) { ftype = type; }
   void SetRun(const std::string &run) { frun = run; }
   void SetTag(const std::string &tag) { ftag = tag; }
   void SetComment(const std::string &comment) { fcomment = comment; }
   void SetLegend(const std::string &legend) { flegend = legend; }
   void SetFileName(const std::string &fileName) { ffileName = fileName; }
   void SetInputFile(const std::string &inputfile) { finputFile = inputfile; }
   void SetOutputROOTfolder(const std::string &outputROOTfolder)
   {
      foutputROOTfolder = outputROOTfolder;
   }

   std::string GetROOTOutputPath() { return fRealpathROOT; }
   float GetMean() { return fph1f_WF->GetMean(); }
   std::string GetCutsList() { return fcutslist; }

private:
   // Data file settings
   std::string ftype;
   std::string frun;
   std::string ftag;
   std::string fcomment;
   std::string flegend;

   std::string finputFile;
   std::string foutputROOTfolder;

   std::string ffileName;
   std::string fcutslist;
   std::string fRealpathROOT;

   // Analysis settings
   int dEdxmax = 1600;
   int nbinsreso = 25;
   int nbinsmean = 30;
   int nbinsresoangle = 30;
   int nbinsmeanangle = 30;

   // Absolute momentum
   float absmomrange = 2000;
   // Mean plots
   int absmommeanindex = 25;
   float nabsmommeanbins = 100;
   float absmommeanbinwidth = absmomrange / nabsmommeanbins;
   // Reso plots
   int absmomresoindex = 0;
   float nabsmomresobins = nbinsreso;
   float absmomresobinwidth = absmomrange / nabsmomresobins;

   // Signed momentum
   float nmombins = 201;
   float momrange = 2000;
   float mombinwidth = 2 * momrange / (nmombins - 1);
   int momindex = 0;

   // X position
   float nxposbins = nbinsmean;
   float xposrange = 2400;
   float xposbinwidth = 2 * xposrange / (nxposbins - 1);
   int xposindex = 0;

   // drift distance
   float ddrange = 1000;
   // Mean plots
   int ddmeanindex = 0;
   float nddmeanbins = nbinsmean;
   float ddmeanbinwidth = ddrange / nddmeanbins;
   // Reso plots
   int ddresoindex = 0;
   float nddresobins = nbinsreso;
   float ddresobinwidth = ddrange / nddresobins;

   // drift time
   float ndtbins = nbinsmean;
   float dtrange = 510;
   float dtbinwidth = dtrange / (ndtbins - 1);
   int dtindex = 0;

   // Track length
   float trklenrange = 1800;
   // Mean plots
   int trklenmeanindex = 0;
   float ntrklenmeanbins = nbinsmean;
   float trklenmeanbinwidth = trklenrange / ntrklenmeanbins;
   // Reso plots
   int trklenresoindex = 0;
   float ntrklenresobins = nbinsreso;
   float trklenresobinwidth = trklenrange / ntrklenresobins;

   // absolute phi angle
   float absphirange = 90;
   // Mean plots
   int absphimeanindex = 0;
   float nabsphimeanbins = nbinsmeanangle;
   float absphimeanbinwidth = absphirange / nabsphimeanbins;
   // Reso plots
   int absphiresoindex = 0;
   float nabsphiresobins = nbinsresoangle;
   float absphiresobinwidth = absphirange / nabsphiresobins;

   // phi angle
   float nphibins = nbinsmeanangle;
   float phirange = 90;
   float phibinwidth = 2 * phirange / (nphibins - 1);
   int phiindex = 0;

   // absolute theta angle
   float absthetarange = 90;
   // Mean plots
   int absthetameanindex = 0;
   float nabsthetameanbins = nbinsmeanangle;
   float absthetameanbinwidth = absthetarange / nabsthetameanbins;
   // Reso plots
   int absthetaresoindex = 0;
   float nabsthetaresobins = nbinsresoangle;
   float absthetaresobinwidth = absthetarange / nabsthetaresobins;

   // theta angle
   float nthetabins = nbinsmeanangle;
   float thetarange = 90;
   float thetabinwidth = 2 * thetarange / (nthetabins - 1);
   int thetaindex = 0;

   // Histograms
   TF1 *fptf1_WFmean = nullptr;
   TF1 *fptf1_XPmean = nullptr;
   TF1 *fptf1_WFreso = nullptr;
   TF1 *fptf1_XPreso = nullptr;
   TF1 *fptf1_pullmu = nullptr;
   TF1 *fptf1_pullelec = nullptr;
   // Vectors for dE/dx vs ERAM module
   std::vector<TH1F *> vmod_fph1f_XP;
   std::vector<TH1F *> vmod_fph1f_WF;
   std::vector<TH2F *> vmod_fph2f_XPdrift;

   // Vectors for mean dE/dx vs absolute momentum bin
   std::vector<TH1F *> vabsmommean_fph1f_WF;
   std::vector<TH1F *> vabsmommean_fph1f_XP;
   TGraphErrors *ptge_absmom_mean_WF = new TGraphErrors();
   TGraphErrors *ptge_absmom_mean_XP = new TGraphErrors();
   std::vector<TH1F *> vabsmomreso_fph1f_WF;
   std::vector<TH1F *> vabsmomreso_fph1f_XP;
   TGraphErrors *ptge_absmom_reso_WF = new TGraphErrors();
   TGraphErrors *ptge_absmom_reso_XP = new TGraphErrors();

   // Vectors for dE/dx vs signed momentum bin
   std::vector<TH1F *> vmom_fph1f_WF;
   std::vector<TH1F *> vmom_fph1f_XP;
   TGraphErrors *ptge_mom_mean_WF = new TGraphErrors();
   TGraphErrors *ptge_mom_mean_XP = new TGraphErrors();
   TGraphErrors *ptge_mom_reso_WF = new TGraphErrors();
   TGraphErrors *ptge_mom_reso_XP = new TGraphErrors();

   // Vectors for pulls vs momentum bin
   std::vector<TH1F *> vmom_fph1f_pullmu;
   std::vector<TH1F *> vmom_fph1f_pullelec;
   TGraphErrors *ptge_mom_mean_pullmu = new TGraphErrors();
   TGraphErrors *ptge_mom_mean_pullelec = new TGraphErrors();
   TGraphErrors *ptge_mom_std_pullmu = new TGraphErrors();
   TGraphErrors *ptge_mom_std_pullelec = new TGraphErrors();

   // Vectors for dE/dx vs X position
   std::vector<TH1F *> vX_fph1f_WF;
   std::vector<TH1F *> vX_fph1f_XP;
   TGraphErrors *ptge_X_mean_WF = new TGraphErrors();
   TGraphErrors *ptge_X_mean_XP = new TGraphErrors();
   TGraphErrors *ptge_X_reso_WF = new TGraphErrors();
   TGraphErrors *ptge_X_reso_XP = new TGraphErrors();

   // Vectors for dE/dx vs drift distance
   std::vector<TH1F *> vddmean_fph1f_WF;
   std::vector<TH1F *> vddmean_fph1f_XP;
   TGraphErrors *ptge_dd_mean_WF = new TGraphErrors();
   TGraphErrors *ptge_dd_mean_XP = new TGraphErrors();
   std::vector<TH1F *> vddreso_fph1f_WF;
   std::vector<TH1F *> vddreso_fph1f_XP;
   TGraphErrors *ptge_dd_reso_WF = new TGraphErrors();
   TGraphErrors *ptge_dd_reso_XP = new TGraphErrors();

   // Vectors for dE/dx vs drift time bin
   std::vector<TH1F *> vdt_fph1f_WF;
   std::vector<TH1F *> vdt_fph1f_XP;
   TGraphErrors *ptge_dt_mean_WF = new TGraphErrors();
   TGraphErrors *ptge_dt_mean_XP = new TGraphErrors();
   TGraphErrors *ptge_dt_reso_WF = new TGraphErrors();
   TGraphErrors *ptge_dt_reso_XP = new TGraphErrors();

   // Vectors for dE/dx vs Track length
   std::vector<TH1F *> vtrklenmean_fph1f_WF;
   std::vector<TH1F *> vtrklenmean_fph1f_XP;
   TGraphErrors *ptge_trklen_mean_WF = new TGraphErrors();
   TGraphErrors *ptge_trklen_mean_XP = new TGraphErrors();
   std::vector<TH1F *> vtrklenreso_fph1f_WF;
   std::vector<TH1F *> vtrklenreso_fph1f_XP;
   TGraphErrors *ptge_trklen_reso_WF = new TGraphErrors();
   TGraphErrors *ptge_trklen_reso_XP = new TGraphErrors();

   // Vectors for dE/dx vs absolute phi angle bin
   std::vector<TH1F *> vabsphimean_fph1f_WF;
   std::vector<TH1F *> vabsphimean_fph1f_XP;
   TGraphErrors *ptge_absphi_mean_WF = new TGraphErrors();
   TGraphErrors *ptge_absphi_mean_XP = new TGraphErrors();
   std::vector<TH1F *> vabsphireso_fph1f_WF;
   std::vector<TH1F *> vabsphireso_fph1f_XP;
   TGraphErrors *ptge_absphi_reso_WF = new TGraphErrors();
   TGraphErrors *ptge_absphi_reso_XP = new TGraphErrors();

   // Vectors for dE/dx vs phi angle bin
   std::vector<TH1F *> vphi_fph1f_WF;
   std::vector<TH1F *> vphi_fph1f_XP;
   TGraphErrors *ptge_phi_mean_WF = new TGraphErrors();
   TGraphErrors *ptge_phi_mean_XP = new TGraphErrors();
   TGraphErrors *ptge_phi_reso_WF = new TGraphErrors();
   TGraphErrors *ptge_phi_reso_XP = new TGraphErrors();

   // Vectors for dE/dx vs absolute theta angle bin
   std::vector<TH1F *> vabsthetamean_fph1f_WF;
   std::vector<TH1F *> vabsthetamean_fph1f_XP;
   TGraphErrors *ptge_abstheta_mean_WF = new TGraphErrors();
   TGraphErrors *ptge_abstheta_mean_XP = new TGraphErrors();
   std::vector<TH1F *> vabsthetareso_fph1f_WF;
   std::vector<TH1F *> vabsthetareso_fph1f_XP;
   TGraphErrors *ptge_abstheta_reso_WF = new TGraphErrors();
   TGraphErrors *ptge_abstheta_reso_XP = new TGraphErrors();

   // Vectors for dE/dx vs theta angle bin
   std::vector<TH1F *> vtheta_fph1f_WF;
   std::vector<TH1F *> vtheta_fph1f_XP;
   TGraphErrors *ptge_theta_mean_WF = new TGraphErrors();
   TGraphErrors *ptge_theta_mean_XP = new TGraphErrors();
   TGraphErrors *ptge_theta_reso_WF = new TGraphErrors();
   TGraphErrors *ptge_theta_reso_XP = new TGraphErrors();

   // Vectors for pulls vs theta angle bin
   std::vector<TH1F *> vtheta_fph1f_pullmu;
   std::vector<TH1F *> vtheta_fph1f_pullelec;
   TGraphErrors *ptge_theta_mean_pullmu = new TGraphErrors();
   TGraphErrors *ptge_theta_mean_pullelec = new TGraphErrors();
   TGraphErrors *ptge_theta_std_pullmu = new TGraphErrors();
   TGraphErrors *ptge_theta_std_pullelec = new TGraphErrors();

   // Base
   int nbinsdEdx = 500;
   TH1F *fph1f_WF = new TH1F("fph1f_WF", ";dE/dx [ADC counts/cm];Count", 100, 0, dEdxmax);
   TH1F *fph1f_XP = new TH1F("fph1f_XP", ";dE/dx [ADC counts/cm];Count", 100, 0, dEdxmax);
   // 2D with dE/dx
   TH2F *fph2f_WFXP = new TH2F(
      "fph1f_WFXP", ";dE/dx with WF [ADC counts/cm];dE/dx with XP [ADC counts/cm]", 100,
      0, 1000, 100, 0, 1000);
   TH2F *fph2f_WFdrift =
      new TH2F("fph2f_WFdrift", ";drift time (timebins);dE/dx with WF [ADC counts/cm]",
               510, 0, 510, nbinsdEdx, 0, 1000);
   TH2F *fph2f_XPdrift =
      new TH2F("fph2f_XPdrift", ";drift time (timebins);dE/dx with XP [ADC counts/cm]",
               510, 0, 510, nbinsdEdx, 0, 1000);
   TH2F *fph2f_WFX =
      new TH2F("fph2f_WFX", ";Track X position [mm];dE/dx with WF [ADC counts/cm]", 100,
               -1200, 1200, nbinsdEdx, 0, 1000);
   TH2F *fph2f_XPX =
      new TH2F("fph2f_XPX", ";Track X position [mm];dE/dx with XP [ADC counts/cm]", 100,
               -1200, 1200, nbinsdEdx, 0, 1000);
   TH2F *fph2f_XPlen =
      new TH2F("fph1f_lenXP", ";Track length [cm];dE/dx with XP [ADC counts/cm]", 171, 0,
               170, nbinsdEdx, 0, 1000);
   TH2F *fph2f_XPphi = new TH2F("fph2f_XPphi", ";#phi; dE/dx with XP [ADC counts/cm]",
                                5 * nphibins, -90, 90, nbinsdEdx, 0, 1000);
   TH2F *fph2f_XPtheta = new TH2F(
      "fph2f_XPtheta", ";Track angle #theta [#circ]; dE/dx with XP [ADC counts/cm]",
      5 * nthetabins, -90, 90, nbinsdEdx, 0, 1000);
   TH2F *fph2f_XPabsmommean =
      new TH2F("fph2f_XPabsmommean", ";Momentum [MeV/c];dE/dx with XP [ADC counts/cm]",
               5 * nabsmommeanbins / 3, 0, absmomrange, nbinsdEdx / 3, 0, 1000);
   TH2F *fph2f_WFabsmommean =
      new TH2F("fph2f_WFabsmommean", ";Momentum [MeV/c];dE/dx with WF [ADC counts/cm]",
               5 * nabsmommeanbins / 3, 0, absmomrange, nbinsdEdx / 3, 0, 1000);
   TH2F *fph2f_XPmom =
      new TH2F("fph2f_XPmom", ";Signed momentum (MeV);dE/dx with XP [ADC counts/cm]",
               5 * nmombins, -momrange, momrange, nbinsdEdx, 0, 1000);
   TH2F *fph2f_WFmom =
      new TH2F("fph2f_WFmom", ";Signed momentum (MeV);dE/dx with WF [ADC counts/cm]",
               5 * nmombins, -momrange, momrange, nbinsdEdx, 0, 1000);

   // Position
   TH1F *fph1f_Xpos =
      new TH1F("fph1f_Xpos", ";Track X position [mm];Count", 200, -1200, 1200);

   // Momentum
   TH1I *fph1i_mom = new TH1I("fph1i_mom", ";Signed momentum (MeV);Count", 3 * nmombins,
                              -momrange, momrange);
   TH1I *fph1i_mom_tHAT = new TH1I("fph1i_mom_thAT", ";Signed momentum (MeV);Count",
                                   3 * nmombins, -momrange, momrange);
   TH1I *fph1i_mom_bHAT = new TH1I("fph1i_mom_bhAT", ";Signed momentum (MeV);Count",
                                   3 * nmombins, -momrange, momrange);
   // Geometry
   TH2F *fph2f_XZ = new TH2F("fph1f_XZ", ";X;Z", 100, -1, 1, 100, -1, 1);
   TH2F *fph2f_YZ = new TH2F("fph1f_YZ", ";Y;Z", 100, -1, 1, 100, -1, 1);
   TH1F *fph1f_phi =
      new TH1F("fph1f_phi", ";Track angle #varphi [#circ] angle;Count", 100, -90, 90);
   TH1F *fph1f_theta =
      new TH1F("fph1f_theta", ";Track angle #theta [#circ];Count", 100, -90, 90);
   TH2F *fph2f_phitheta = new TH2F(
      "fph2f_phitheta", ";Track angle #theta [#circ];Track angle #varphi [#circ]", 100,
      -90, 90, 100, -90, 90);
   TH1F *fph1f_trklen = new TH1F("fph1f_trklen", ";Track length [cm];Count", 171, 0, 170);
   TH1F *fph1f_chi2 = new TH1F("fph1f_chi2", ";#chi^{2};Count", 1000, 0, 50);
   TH2F *fph2f_chi2mom = new TH2F("fph2f_chi2mom", ";Signed momentum (MeV);#chi^{2}",
                                  nmombins, -momrange, momrange, 1000, 0, 50);
   // Time in bHAT
   TH1I *fph1i_tminBotCath =
      new TH1I("fph1i_tminBotCath", ";time bin;Count", 510, 0, 510);
   TH1I *fph1i_tminEP0 =
      new TH1I("fph1i_tminEP0", "Start time in EP0;time bin;Count", 510, 0, 510);
   TH1I *fph1i_tminEP1 =
      new TH1I("fph1i_tminEP1", "Start time in EP1;time bin;Count", 510, 0, 510);
   TH1I *fph1i_tmaxBotCath =
      new TH1I("fph1i_tmaxBotCath", ";time bin;Count", 510, 0, 510);
   TH1I *fph1i_tmaxEP0 =
      new TH1I("fph1i_tmaxEP0", "End time in EP0;time bin;Count", 510, 0, 510);
   TH1I *fph1i_tmaxEP1 =
      new TH1I("fph1i_tmaxEP1", "End time in EP1;time bin;Count", 510, 0, 510);
   // Time in tHAT
   TH1I *fph1i_tminTopCath =
      new TH1I("fph1i_tminTopCath", "Start time in tHATPC;time bin;Count", 510, 0, 510);
   TH1I *fph1i_tminEP2 =
      new TH1I("fph1i_tminEP2", "Start time in EP2;time bin;Count", 510, 0, 510);
   TH1I *fph1i_tminEP3 =
      new TH1I("fph1i_tminEP3", "Start time in EP3;time bin;Count", 510, 0, 510);
   TH1I *fph1i_tmaxTopCath =
      new TH1I("fph1i_tmaxTopCath", "End time in tHATPC;time bin;Count", 510, 0, 510);
   TH1I *fph1i_tmaxEP2 =
      new TH1I("fph1i_tmaxEP2", "End time in EP2;time bin;Count", 510, 0, 510);
   TH1I *fph1i_tmaxEP3 =
      new TH1I("fph1i_tmaxEP3", "End time in EP3;time bin;Count", 510, 0, 510);
   // Debug
   TH1F *fph1f_pullmu = new TH1F("fph1f_pullmu", ";Pull #mu;Count", 100, -25, 25);
   TH1F *fph1f_pullelec = new TH1F("fph1f_pullelec", ";Pull e;Count", 100, -25, 25);
   TH1F *fph1f_pullproton = new TH1F("fph1f_pullproton", ";Pull p;Count", 100, -25, 25);
   TH2F *fph2f_pullelecmu =
      new TH2F("ph2f_pullelecmu", ";Pull #mu;Pull e", 100, -25, 25, 100, -25, 25);
   TH2F *fph2f_chi2ndfphi =
      new TH2F("ph2f_chi2ndfphi", ";Track angle #varphi [#circ];#chi^{2}/NDF", 100, -90,
               90, 100, 0, 50);
   TH2F *fph2f_momtheta =
      new TH2F("ph2f_momtheta", ";Track angle #theta [#circ];Signed momentum (MeV)",
               5 * nthetabins, -90, 90, 5 * nmombins, -momrange, momrange);
   TH2F *fph2f_momphi =
      new TH2F("ph2f_momphi", ";Track angle #varphi [#circ];Signed momentum (MeV)",
               5 * nphibins, -90, 90, 5 * nmombins, -momrange, momrange);
   TH2F *fph2f_momR = new TH2F("ph2f_momR", ";R;Signed momentum (MeV)", 300, -5e4, 5e4,
                               3 * nmombins, -momrange, momrange);
   TH2F *fph2f_chi2ndfR =
      new TH2F("ph2f_chi2ndfR", ";R;#chi^{2}/NDF", 300, -5e4, 5e4, 300, 0, 50);
   TH2F *fph2f_lentheta =
      new TH2F("ph2f_lentheta", ";Track angle #theta [#circ];Track length [cm]",
               5 * nthetabins, -90, 90, 171, 0, 170);
   TH2F *fph2f_lenphi =
      new TH2F("ph2f_lenphi", ";Track angle #varphi [#circ];Track length [cm]",
               5 * nphibins, -90, 90, 171, 0, 170);
   TH2F *fph2f_timeX = new TH2F("ph2f_timeX", ";X;time", 100, -1200, 1200, 100, 0, 510);
   TH1F *fph1f_dir0 = new TH1F("ph1f_dir0", ";X direction;Counts", 100, 0, 1);
   TH1F *fph1f_dir1 = new TH1F("ph1f_dir1", ";Y direction;Counts", 100, 0, 1);
   TH1F *fph1f_dir2 = new TH1F("ph1f_dir2", ";Z direction;Counts", 100, 0, 1);
   TH2F *fph2f_momlen =
      new TH2F("ph2f_momlen", ";Track length [cm];Signed momentum (MeV)", 171, 0, 170,
               3 * nmombins, 0, momrange);
   TH2F *fph2f_momncl = new TH2F("ph2f_momncl", ";N_{clusters};Signed momentum (MeV)",
                                 200, 0, 200, 3 * nmombins, 0, momrange);

   // Tree variables
   Double_t wf;
   Double_t xp;
   Double_t xpSystFitRelat;
   Double_t xpSystRCRelat;
   Double_t xpSystLUTRelat;
   Double_t dx;
   Double_t pos[3];
   Double_t dir[3];
   float phi;
   float theta;
   Double_t curv;
   Double_t chi2;
   Double_t NDF;
   Double_t mom_og;
   Double_t mom;
   Double_t start_time;
   Double_t mean_time;
   Double_t end_time;
   Double_t pull_muon;
   Double_t pull_ele;
   Double_t pull_proton;
   Int_t ncl;
   Int_t endplate;
   Int_t eram_channel;
   Int_t eram_ID;
   Int_t T0type;

   // Cuts
   int nclmin{0};
   int nclmax{1000};
   int dxmin{0};
   int dxmax{1000000};
   int ncroscut{0};
   int momcutlow{0};
   int momcuthigh{1000000000};
   int xcutmin{-1200};
   int xcutmax{1200};
   int tstartmin{0};
   int tstartmax{1000};
   int tmeanmin{0};
   int tmeanmax{1000};
   int chi2max{1000000000};
   int hat{0};
   int pullmumax{1000000};
   int phimin{-90};
   int phimax{90};
   int thetamin{-90};
   int thetamax{90};

   // Friend declarations to allow inner classes to access private members
   friend class Draw;
};

#endif