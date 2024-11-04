#ifndef PROCESS_H
#define PROCESS_H


#include "Util.h"

class Process
{
	public:
		Process();
		~Process();
		// Member functions
		void Datafile(const std::string &comment, const std::string &datafilepath, const std::string &drawfolderpath, const std::string &filename);

	private:
		
		// Data file settings
		std::string fdrawfolderpath;
		std::string ffilename;
		std::string fdrawfile;
		std::string fcomment;
		TFile *fpFile;
		TTree *fpTree;

		// Analysis settings
		int xmax =								1300;

		int nmombins =							101;
		int momrange =							3000;
		int mombinwidth =						2*momrange/(nmombins-1);
		int momindex =							0;

		int nddbins =							101;
		int ddrange =							510;
		int ddbinwidth =						ddrange/(nddbins-1);
		int ddindex =							0;

		int nphibins =							91;
		int phirange =							90;
		int phibinwidth =						2*phirange/(nphibins-1);
		int phiindex =							0;

		int nthetabins =						91;
		int thetarange =						90;
		int thetabinwidth =						2*thetarange/(nthetabins-1);
		int thetaindex =						0;


		// Histograms
		TF1  *fptf1_WF;
		TF1  *fptf1_XP;
		// Vectors for dE/dx per ERAM module
		std::vector<TH1F*>	vmod_fph1f_XP;
		std::vector<TH1F*>	vmod_fph1f_WF;
		std::vector<TH2F*>	vmod_fph2f_XPtmean;
		std::vector<TH2F*>	vmod_fph2f_XPdrift;

		// Vectors for dE/dx per momentum bin
		std::vector<TH1F*> vmom_fph1f_WF;
		std::vector<TH1F*> vmom_fph1f_XP;
		TGraphErrors *ptge_mom_mean_WF =		new TGraphErrors();
		TGraphErrors *ptge_mom_mean_XP =		new TGraphErrors();
		TGraphErrors *ptge_mom_reso_WF =		new TGraphErrors();
		TGraphErrors *ptge_mom_reso_XP =		new TGraphErrors();

		// Vectors for dE/dx per drift distance bin
		std::vector<TH1F*> vdd_fph1f_WF;
		std::vector<TH1F*> vdd_fph1f_XP;
		TGraphErrors *ptge_dd_mean_WF =			new TGraphErrors();
		TGraphErrors *ptge_dd_mean_XP =			new TGraphErrors();
		TGraphErrors *ptge_dd_reso_WF =			new TGraphErrors();
		TGraphErrors *ptge_dd_reso_XP =			new TGraphErrors();

		// Vectors for dE/dx per phi angle bin
		std::vector<TH1F*> vphi_fph1f_WF;
		std::vector<TH1F*> vphi_fph1f_XP;
		TGraphErrors *ptge_phi_mean_WF =		new TGraphErrors();
		TGraphErrors *ptge_phi_mean_XP =		new TGraphErrors();
		TGraphErrors *ptge_phi_reso_WF =		new TGraphErrors();
		TGraphErrors *ptge_phi_reso_XP =		new TGraphErrors();

		// Vectors for dE/dx per theta angle bin
		std::vector<TH1F*> vtheta_fph1f_WF;
		std::vector<TH1F*> vtheta_fph1f_XP;
		TGraphErrors *ptge_theta_mean_WF =		new TGraphErrors();
		TGraphErrors *ptge_theta_mean_XP =		new TGraphErrors();
		TGraphErrors *ptge_theta_reso_WF =		new TGraphErrors();
		TGraphErrors *ptge_theta_reso_XP =		new TGraphErrors();


		// Base
		int nbinsdEdx =							500;
		TH1F *fph1f_WF =						new TH1F("fph1f_WF",		";dE/dx(ADC counts/cm);Count", 100, 0, xmax);
		TH1F *fph1f_XP =						new TH1F("fph1f_XP",		";dE/dx(ADC counts/cm);Count", 100, 0, xmax);
		// 2D with dE/dx		
		TH2F *fph2f_WFXP =						new TH2F("fph1f_WFXP",		";dE/dx with WF (ADC counts/cm);dE/dx with XP (ADC counts/cm)", 100, 0, 1000, 100, 0, 1000);
		TH2F *fph2f_XPdrift =					new TH2F("fph2f_XPdrift",	";drift time (timebins);dE/dx with WF (ADC counts/cm)", 510, 0, 510, nbinsdEdx, 0, 1000);
		TH2F *fph2f_XPlen =						new TH2F("fph1f_lenXP",		";track length (cm);dE/dx with XP (ADC counts/cm)", 171, 0, 170, nbinsdEdx, 0, 1000);
		TH2F *fph2f_XPphi = 					new TH2F("fph2f_XPphi",		";#phi; dE/dx with XP (ADC counts/cm)", 5*nphibins, -90, 90, nbinsdEdx, 0, 1000);
		TH2F *fph2f_XPtheta =					new TH2F("fph2f_XPtheta",	";#theta; dE/dx with XP (ADC counts/cm)", 5*nthetabins, -90, 90, nbinsdEdx, 0, 1000);
		TH2F *fph2f_XPmom =						new TH2F("fph2f_XPmom",		";momentum (MeV);dE/dx with XP (ADC counts/cm)", 5*nmombins, -momrange, momrange, nbinsdEdx, 0, 1000);
		TH2F *fph2f_WFmom =						new TH2F("fph2f_WFmom",		";momentum (MeV);dE/dx with WF (ADC counts/cm)", 5*nmombins, -momrange, momrange, nbinsdEdx, 0, 1000);
		// Momentum		
		TH1I *fph1i_mom =						new TH1I("fph1i_mom",		";momentum (MeV);Count", 3*nmombins, -momrange, momrange);
		TH1I *fph1i_mom_tHAT =					new TH1I("fph1i_mom_thAT",	";momentum (MeV);Count", 3*nmombins, -momrange, momrange);
		TH1I *fph1i_mom_bHAT =					new TH1I("fph1i_mom_bhAT",	";momentum (MeV);Count", 3*nmombins, -momrange, momrange);
		// Geometry		
		TH2F *fph2f_XZ =						new TH2F("fph1f_XZ",		";X;Z", 100, -1, 1, 100, -1, 1);
		TH1F *fph1f_phi =						new TH1F("fph1f_phi",		";#varphi angle;Count", 100, -90, 90);
		TH1F *fph1f_theta =						new TH1F("fph1f_theta",		";#theta angle;Count", 100, -90, 90);
		TH2F *fph2f_phitheta =					new TH2F("fph2f_phitheta",	";#theta;#varphi", 100, -90, 90, 100, -90, 90);
		TH1F *fph1f_trklen = 					new TH1F("fph1f_trklen",	";track length (cm);Count", 171, 0, 170);
		TH1F *fph1f_chi2 =						new TH1F("fph1f_chi2",		";#chi^{2};Count", 1000, 0, 50);
		TH2F *fph2f_chi2mom = 					new TH2F("fph2f_chi2mom",	";momentum (MeV);#chi^{2}", nmombins, -momrange, momrange, 1000, 0, 50);
		// Time in bHAT		
		TH1I *fph1i_tminBotCath =				new TH1I("fph1i_tminBotCath",";time bin;Count", 510, 0, 510);
		TH1I *fph1i_tminEP0 =					new TH1I("fph1i_tminEP0",	"Start time in EP0;time bin;Count", 510, 0, 510);
		TH1I *fph1i_tminEP1 =					new TH1I("fph1i_tminEP1",	"Start time in EP1;time bin;Count", 510, 0, 510);
		TH1I *fph1i_tmaxBotCath =				new TH1I("fph1i_tmaxBotCath",";time bin;Count", 510, 0, 510);
		TH1I *fph1i_tmaxEP0 =					new TH1I("fph1i_tmaxEP0",	"End time in EP0;time bin;Count", 510, 0, 510);
		TH1I *fph1i_tmaxEP1 =					new TH1I("fph1i_tmaxEP1",	"End time in EP1;time bin;Count", 510, 0, 510);
		// Time in tHAT		
		TH1I *fph1i_tminTopCath =				new TH1I("fph1i_tminTopCath","Start time in tHATPC;time bin;Count", 510, 0, 510);
		TH1I *fph1i_tminEP2 =					new TH1I("fph1i_tminEP2",	"Start time in EP2;time bin;Count", 510, 0, 510);
		TH1I *fph1i_tminEP3 =					new TH1I("fph1i_tminEP3",	"Start time in EP3;time bin;Count", 510, 0, 510);
		TH1I *fph1i_tmaxTopCath =				new TH1I("fph1i_tmaxTopCath","End time in tHATPC;time bin;Count", 510, 0, 510);
		TH1I *fph1i_tmaxEP2 =					new TH1I("fph1i_tmaxEP2",	"End time in EP2;time bin;Count", 510, 0, 510);
		TH1I *fph1i_tmaxEP3 =					new TH1I("fph1i_tmaxEP3",	"End time in EP3;time bin;Count", 510, 0, 510);
		// Debug
		TH2F *fph2f_pullemu = 					new TH2F("ph2f_pullemu",	";Pull #mu;Pull e", 100, -10, 10, 100, -10, 10);
		TH2F *fph2f_chi2ndfphi = 				new TH2F("ph2f_chi2ndfphi",	";#varphi;#chi^{2}/NDF", 100, -90, 90, 100, 0, 50);
		TH2F *fph2f_momtheta = 					new TH2F("ph2f_momtheta",	";#theta;momentum (MeV)", 5*nthetabins, -90, 90, 5*nmombins, -momrange, momrange);
		TH2F *fph2f_momR = 						new TH2F("ph2f_momR",		";R;momentum (MeV)", 300, -5e4, 5e4, 3*nmombins, -momrange, momrange);
		TH2F *fph2f_chi2ndfR = 					new TH2F("ph2f_chi2ndfR",	";R;#chi^{2}/NDF", 300, -5e4, 5e4, 300, 0, 50);
		TH2F *fph2f_lentheta = 					new TH2F("ph2f_lentheta",	";#theta;track length (cm)", 5*nthetabins, -90, 90, 171, 0, 170);
		TH2F *fph2f_lenphi = 					new TH2F("ph2f_lenphi",		";#varphi;track length (cm)", 5*nphibins, -90, 90, 171, 0, 170);

		// Tree variables
		Double_t wf;
		Double_t xp;
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
		Double_t APM;
		Double_t start_time;
		Double_t mean_time;
		Double_t end_time;
		Double_t pull_muon;
		Double_t pull_ele;
		Int_t	ncl;
		Int_t	endplate;
		Int_t	eram_channel;
		Int_t	eram_ID;

		// Friend declarations to allow inner classes to access private members
		friend class Draw;
};

#endif