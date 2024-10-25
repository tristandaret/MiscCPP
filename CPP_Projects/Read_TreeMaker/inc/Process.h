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
		TGraphErrors *ptge_mom_std_WF =			new TGraphErrors();
		TGraphErrors *ptge_mom_std_XP =			new TGraphErrors();
		TGraphErrors *ptge_mom_reso_WF =		new TGraphErrors();
		TGraphErrors *ptge_mom_reso_XP =		new TGraphErrors();

		// Base
		TH1F *fph1f_WF =						new TH1F("fph1f_WF",		";dE/dx(ADC counts/cm);Count", 100, 0, xmax);
		TH1F *fph1f_XP =						new TH1F("fph1f_XP",		";dE/dx(ADC counts/cm);Count", 100, 0, xmax);
		// 2D with dE/dx		
		TH2F *fph2f_WFXP =						new TH2F("fph1f_WFXP",		";dE/dx with WF (ADC counts/cm);dE/dx with XP (ADC counts/cm)", 100, 0, 1000, 100, 0, 1000);
		TH2F *fph2f_XPdrift =					new TH2F("fph2f_XPdrift",	";drift time (timebins);dE/dx with WF (ADC counts/cm)", 510, 0, 510, 100, 0, 1000);
		TH2F *fph2f_lenXP =						new TH2F("fph1f_lenXP",		";track length (cm);dE/dx with XP (ADC counts/cm)", 171, 0, 170, 100, 0, 1000);
		// Momentum		
		TH2F *fph2f_WFmom =						new TH2F("fph2f_WFmom",		";momentum (MeV);dE/dx with WF (ADC counts/cm)", 5*nmombins, -momrange, momrange, 5*100, 0, 1000);
		TH2F *fph2f_XPmom =						new TH2F("fph2f_XPmom",		";momentum (MeV);dE/dx with XP (ADC counts/cm)", 5*nmombins, -momrange, momrange, 5*100, 0, 1000);
		TH1I *fph1i_mom =						new TH1I("fph1i_mom",		";momentum (MeV);Count", 3*nmombins, -momrange, momrange);
		TH1I *fph1i_mom_tHAT =					new TH1I("fph1i_mom_thAT",	";momentum (MeV);Count", 3*nmombins, -momrange, momrange);
		TH1I *fph1i_mom_bHAT =					new TH1I("fph1i_mom_bhAT",	";momentum (MeV);Count", 3*nmombins, -momrange, momrange);
		// Geometry		
		TH2F *fph2f_XZ =						new TH2F("fph1f_XZ",		";X;Z", 100, -1, 1, 100, -1, 1);
		TH1F *fph1f_dirY =						new TH1F("fph1f_dirY",		";Y direction;Count", 100, -1, 1);
		TH1F *fph1f_trklen = 					new TH1F("fph1f_trklen",	";track length (cm);Count", 171, 0, 170);
		TH1F *fph1f_chi2 =						new TH1F("fph1f_chi2",		";#chi^{2};Count", 1000, 0, 1e4);
		TH2F *fph1f_chi2mom = 					new TH2F("fph1f_chi2mom",	";momentum (MeV);#chi^{2}", nmombins, -momrange, momrange, 1000, 0, 1e4);
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
		TH2F *fph2f_dirYposY = 					new TH2F("fph2f_dirYposY",	";Y direction;Y position", 100, -1, 1, 100, -1e3, 1e3);
		TH2F *fph2f_dirYmom = 					new TH2F("fph2f_dirYmom",	";Y direction;Momentum", 100, -1, 1, nmombins, -momrange, momrange);
		TH2F *fph2f_momposY = 					new TH2F("fph2f_momposY",	";Momentum;Y position", nmombins, -momrange, momrange, 100, -1e3, 1e3);
		TH2F *fph2f_pullemu = 					new TH2F("ph2f_pullemu",	";Pull #mu;Pull e", 100, -10, 10, 100, -10, 10);

		// Tree variables
		Double_t wf;
		Double_t xp;
		Double_t dx;
		Double_t pos[3];
		Double_t dir[3];
		Double_t chi2;
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