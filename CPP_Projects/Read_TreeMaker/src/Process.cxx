#include "Process.h"
#include "Misc_Functions.h"
#include "TFrame.h"
#include "SetStyle.h"

#include <fstream>

// Constructor
Process::Process(){

	std::vector<int> channel2iD =	{	24, 30, 28, 19, 21, 13,	9,	2, 26, 17, 23, 29,	1, 10, 11,	3,	/*bottom HATPC*/ 
										47,	16, 14, 15, 42, 45, 37, 36, 20, 38,	7, 44, 43, 39, 41, 46}; /*top	HATPC*/

	for(int i = 0; i < 32; i++){
		vmod_fph1f_XP.						push_back(new TH1F(Form("fph1f_XP_%d", i), Form("Energy loss in ERAM %i;dE/dx (ADC counts/cm);Count", channel2iD[i]), 100, 0, xmax));
		vmod_fph1f_WF.						push_back(new TH1F(Form("fph1f_WF_%d", i), Form("Energy loss in ERAM %i;dE/dx (ADC counts/cm);Count", channel2iD[i]), 100, 0, xmax));
		vmod_fph2f_XPtmean.					push_back(new TH2F(Form("fph2f_XPtmean_%d", i), Form("Energy loss (XP) vs mean time in ERAM %i;mean time (ns);dE/dx (ADC counts/cm)", channel2iD[i]), 510, 0, 510, 100, 0, xmax));
		vmod_fph2f_XPdrift.					push_back(new TH2F(Form("fph2f_XPdrift_%d", i), Form("Energy loss (XP) vs drift time in ERAM %i;drift time (timebins);dE/dx (ADC counts/cm)", channel2iD[i]), 510, 0, 510, 100, 0, xmax));
	}

	for(int i=0; i<nmombins;i++){
		int mommin =						i*mombinwidth - momrange;
		int mommax =						(i+1)*mombinwidth - momrange;
		vmom_fph1f_WF.						push_back(new TH1F(Form("vmom_fph1f_WF_%d_%d", mommin, mommax), Form("Energy loss | %d < p < %d; dE/dx (ADC counts/cm); Count", mommin, mommax), 100, 0, xmax));
		vmom_fph1f_XP.						push_back(new TH1F(Form("vmom_fph1f_XP_%d_%d", mommin, mommax), Form("Energy loss | %d < p < %d; dE/dx (ADC counts/cm); Count", mommin, mommax), 100, 0, xmax));

		arr_momdd_fph1f_XP.					push_back(std::vector<TH1F*>());
		for(int i=0; i<nddbins;i++){
			int ddmin =							i*ddbinwidth;
			int ddmax =							(i+1)*ddbinwidth;
			arr_momdd_fph1f_XP.back().			push_back(new TH1F(Form("arr_momdd_fph1f_XP_%d_%d_%d_%d", mommin, mommax, ddmin, ddmax), Form("Energy loss | %d < p < %d and %d < dd < %d; dE/dx (ADC counts/cm); Count", mommin, mommax, ddmin, ddmax), 100, 0, xmax));
		}
	}
}

// Destructor
Process::~Process(){
}

void Process::Datafile(const std::string &comment, const std::string &datafilepath, const std::string &drawfolderpath, const std::string &filename){

	fdrawfolderpath =						drawfolderpath;
	ffilename =								filename;
	fdrawfile =								drawfolderpath + "/" + filename;
	fcomment =								comment;
	fpFile =								TFile::Open(datafilepath.c_str());
	fpTree =								(TTree*)fpFile->Get("outTree");
	int nentries =							fpTree->GetEntries();
	int dirYnegbHAT=0, dirYposbHAT=0, dirYnegtHAT=0, dirYpostHAT=0;
	int nbotneg=0, nbotpos=0, ntopneg=0, ntoppos=0;

	fpTree->SetBranchAddress("dEdx_WF", &wf);
	fpTree->SetBranchAddress("dEdx_XP", &xp);
	fpTree->SetBranchAddress("track_length", &dx);
	fpTree->SetBranchAddress("pos", &pos);
	fpTree->SetBranchAddress("dir", &dir);
	fpTree->SetBranchAddress("chi2", &chi2);
	fpTree->SetBranchAddress("NDF", &NDF);
	fpTree->SetBranchAddress("mom", &mom_og);
	fpTree->SetBranchAddress("avg_pad_mult", &APM);
	fpTree->SetBranchAddress("starttime", &start_time);
	fpTree->SetBranchAddress("mean_time", &mean_time);
	fpTree->SetBranchAddress("endtime",	&end_time);
	fpTree->SetBranchAddress("ncl", &ncl);
	fpTree->SetBranchAddress("EndPlate_ID", &endplate);
	fpTree->SetBranchAddress("eram_channel", &eram_channel);
	fpTree->SetBranchAddress("eram_ID", &eram_ID);
	fpTree->SetBranchAddress("pull_muon", &pull_muon);
	fpTree->SetBranchAddress("pull_ele", &pull_ele);

	// Cuts and output file 
	int nclmin=0, nclmax=200, dxmin=0, dxmax=2e3, ncroscut=0;
	int apmcutlow=0, apmcuthigh=20, momcutlow=0, momcuthigh=1e6, tcutmin=0, tcutmax=1e3;
	int chi2max = 1e3, hat=0, pullmumax=100;
	float dirminY=-1, dirmaxY=1;
	// nclmin = 32;							fdrawfile += ("_" + std::to_string(nclmin)		+ "ncl");
	// apmcutlow = 2, apmcuthigh = 4;		fdrawfile += ("_" + std::to_string(apmcutlow)	+ "apm" + std::to_string(apmcuthigh));
	// momcutlow = 250, momcuthigh=600;		fdrawfile += ("_" + std::to_string(momcutlow)	+ "mom" + std::to_string(momcuthigh));
	// momcutlow = 1e2;						fdrawfile += ("_" + std::to_string(momcutlow)	+ "mom");
	// tcutmin = 0, tcutmax = 75;			fdrawfile += ("_" + std::to_string(tcutmin)		+ "tmin" + std::to_string(tcutmax));
	// nclmin=50, nclmax = 150;				fdrawfile += ("_" + std::to_string(nclmin)		+ "ncl" + std::to_string(nclmax));
	// dxmin = 50, dxmax = 150;			fdrawfile += ("_" + std::to_string(dxmin) 		+ "dx" + std::to_string(dxmax));
	dxmin = 25;							fdrawfile += ("_" + std::to_string(dxmin) + "dx");
	// dirminY = 0.9; dirmaxY = 1;			fdrawfile += ("_" + std::to_string(dirminY)		+ "dirY" + std::to_string(dirmaxY));
	chi2max = 5;						fdrawfile += ("_chi2ndf" + std::to_string(chi2max));
	pullmumax = 2;						fdrawfile += ("_pullmu" + std::to_string(pullmumax));
	// fdrawfile += 						"_dir1>0flip";
	// hat = -1;							fdrawfile += std::string("_") + (hat == -1 ? "bHAT" : "tHAT");
	fdrawfile += 						".pdf";

	for(int i = 0; i < nentries; i++){
		fpTree->GetEntry(i);
		// FOR COSMICS if the track is in tHAT and going up, the sign must be flipped
		// if(fdrawfile.find("cosmic") != std::string::npos)
			// (pos[1] > 0 and dir[1] > 0) ? mom = -mom_og : mom = mom_og;
		// dir[1] > 0 ? mom = -mom_og : mom = mom_og;
		// pos[1] > 0 ? mom = -mom_og : mom = mom_og;
		mom = mom_og;
		if(fabs(mom) < 1 || std::isnan(mom)) continue;

		if(nclmin > ncl or ncl > nclmax) continue;
		if(dxmin > dx/10 or dxmax < dx/10) continue;
		if(fabs(dir[1]) < dirminY or dirmaxY < fabs(dir[1])) continue;
		if(dx < dxmin) continue;
		if(APM < apmcutlow or APM > apmcuthigh) continue;
		if(fabs(mom) < momcutlow or fabs(mom) > momcuthigh) continue;
		if(mean_time < tcutmin or mean_time > tcutmax) continue;
		if(chi2/NDF > chi2max) continue;
		if(hat != 0 and sign(hat) != sign(pos[1])) continue;
		// if(fabs(pull_muon) > fabs(pull_ele)) continue;
		if(fabs(pull_muon) > pullmumax) continue;

		fph1f_WF->							Fill(wf/1.019);
		fph1f_XP->							Fill(xp);
		fph2f_WFXP->						Fill(wf/1.019, xp);
		fph2f_XPlen->						Fill(dx/10, xp);
		fph2f_XPdrift->						Fill(mean_time, xp);
		fph2f_XPphi->						Fill(TMath::ASin(dir[1])*180/TMath::Pi(), xp);
		vmod_fph2f_XPdrift[eram_channel]->	Fill(mean_time, xp);
		vmod_fph1f_WF[eram_channel]->		Fill(wf/1.019);
		vmod_fph1f_XP[eram_channel]->		Fill(xp);
		vmod_fph2f_XPtmean[eram_channel]->	Fill(mean_time, xp);

		fph2f_XZ->							Fill(dir[0], dir[2]);
		fph1f_dirY->						Fill(dir[1]);
		fph1f_trklen->						Fill(dx/10);
		fph1f_chi2->						Fill(chi2/NDF);
		fph1f_chi2mom->						Fill(mom, chi2/NDF);

		momindex =							(int)std::round(mom/mombinwidth)+nmombins/2;
		if(fabs(mom) < momrange){
			vmom_fph1f_WF[momindex]->			Fill(wf/1.019);
			vmom_fph1f_XP[momindex]->			Fill(xp);
		}
		fph1i_mom->							Fill(mom);
		if(pos[1] < 0) fph1i_mom_bHAT->		Fill(mom);
		if(pos[1] > 0) fph1i_mom_tHAT->		Fill(mom);
		fph2f_WFmom->						Fill(mom, wf/1.019);
		fph2f_XPmom->						Fill(mom, xp);

		ddindex = (int)std::round(mean_time/ddbinwidth);
		if(mean_time > 0 and mean_time < 300 and fabs(mom) < momrange){
			arr_momdd_fph1f_XP[momindex][ddindex]->Fill(xp);
		}

		if(eram_channel < 16){
			if(endplate == 8){
				fph1i_tminEP0->Fill(start_time);
				fph1i_tmaxEP0->Fill(end_time);
			}
			if(endplate == 4) {
				fph1i_tminEP1->Fill(start_time);
				fph1i_tmaxEP1->Fill(end_time);
			}
			if(endplate == 12){
				fph1i_tminBotCath->Fill(start_time);
				fph1i_tmaxBotCath->Fill(end_time);
			}
		}
		else{
			if(endplate == 8){
				fph1i_tminEP2->Fill(start_time);
				fph1i_tmaxEP2->Fill(end_time);
			}
			if(endplate == 4) {
				fph1i_tminEP3->Fill(start_time);
				fph1i_tmaxEP3->Fill(end_time);
			}
			if(endplate == 12){
				fph1i_tminTopCath->Fill(start_time);
				fph1i_tmaxTopCath->Fill(end_time);
			}
		}

		// Debug
		if(dir[1] < 0 and pos[1] < 0) dirYnegbHAT++;
		if(dir[1] > 0 and pos[1] < 0) dirYposbHAT++;
		if(dir[1] < 0 and pos[1] > 0) dirYnegtHAT++;
		if(dir[1] > 0 and pos[1] > 0) dirYpostHAT++;

		if(pos[1] < 0 and sign(mom) < 0) nbotneg++;
		if(pos[1] < 0 and sign(mom) > 0) nbotpos++;
		if(pos[1] > 0 and sign(mom) < 0) ntopneg++;
		if(pos[1] > 0 and sign(mom) > 0) ntoppos++;
		fph2f_dirYmom->						Fill(dir[1], mom);
		fph2f_momposY->						Fill(mom, pos[1]);
		fph2f_dirYposY->					Fill(dir[1], pos[1]);
		fph2f_pullemu->						Fill(pull_muon, pull_ele);
	}
	std::cout << "Y direction:" << std::endl;
	std::cout << "tHAT: negative => " << dirYnegtHAT << " (" << dirYnegtHAT*100.0/(dirYnegtHAT+dirYpostHAT) << "%) positive => " << dirYpostHAT << " (" << dirYpostHAT*100.0/(dirYnegtHAT+dirYpostHAT) << "%)" << std::endl;
	std::cout << "bHAT: negative => " << dirYnegbHAT << " (" << dirYnegbHAT*100.0/(dirYnegbHAT+dirYposbHAT) << "%) positive => " << dirYposbHAT << " (" << dirYposbHAT*100.0/(dirYnegbHAT+dirYposbHAT) << "%)" << std::endl;
	std::cout << "Total: negative => " << dirYnegbHAT+dirYnegtHAT << " (" << (dirYnegbHAT+dirYnegtHAT)*100.0/(dirYnegbHAT+dirYposbHAT+dirYnegtHAT+dirYpostHAT) << "%) positive => " << dirYposbHAT+dirYpostHAT << " (" << (dirYposbHAT+dirYpostHAT)*100.0/(dirYnegbHAT+dirYposbHAT+dirYnegtHAT+dirYpostHAT) << "%)" << std::endl;

	std::cout << "Momentum:" << std::endl;
	std::cout << "Negative: " << ntopneg << " (top) => " << nbotneg << " (bottom) => top-bottom = " << ntopneg-nbotneg << std::endl;
	std::cout << "Positive: " << ntoppos << " (top) => " << nbotpos << " (bottom) => top-bottom = " << ntoppos-nbotpos << std::endl;

	// TGrafph filling
	// Momentum
	int ivalid = 0;
	for(int i=0; i<nmombins; i++){
		int nentries_here =			vmom_fph1f_WF[i]->GetEntries();
		if(nentries_here < 50) continue;
		fptf1_WF =					Fit1Gauss(vmom_fph1f_WF[i]);
		fptf1_XP =					Fit1Gauss(vmom_fph1f_XP[i]);

		float mean_WF =				fptf1_WF->GetParameter(1);
		float mean_XP =				fptf1_XP->GetParameter(1);
		float dmean_WF =			fptf1_WF->GetParError(1);
		float dmean_XP =			fptf1_XP->GetParError(1);

		float std_WF =				fptf1_WF->GetParameter(2);
		float std_XP =				fptf1_XP->GetParameter(2);
		float dstd_WF =				fptf1_WF->GetParError(2);
		float dstd_XP =				fptf1_XP->GetParError(2);

		float reso_WF =				std_WF/mean_WF*100;
		float reso_XP =				std_XP/mean_XP*100;
		float dreso_WF =			GetResoError(fptf1_WF);
		float dreso_XP =			GetResoError(fptf1_XP);

		if(mean_WF == 0 || mean_XP == 0) continue;

		ptge_mom_mean_WF->			SetPoint(ivalid, i*mombinwidth-momrange, mean_WF);
		ptge_mom_mean_WF->			SetPointError(ivalid, mombinwidth/2, dmean_WF);
		ptge_mom_mean_XP->			SetPoint(ivalid, i*mombinwidth-momrange, mean_XP);
		ptge_mom_mean_XP->			SetPointError(ivalid, mombinwidth/2, dmean_XP);

		ptge_mom_std_WF->			SetPoint(ivalid, i*mombinwidth-momrange, std_WF);
		ptge_mom_std_WF->			SetPointError(ivalid, mombinwidth/2, dstd_WF);
		ptge_mom_std_XP->			SetPoint(ivalid, i*mombinwidth-momrange, std_XP);
		ptge_mom_std_XP->			SetPointError(ivalid, mombinwidth/2, dstd_XP);

		ptge_mom_reso_WF->			SetPoint(ivalid, i*mombinwidth-momrange, reso_WF);
		ptge_mom_reso_WF->			SetPointError(ivalid, mombinwidth/2, dreso_WF);
		ptge_mom_reso_XP->			SetPoint(ivalid, i*mombinwidth-momrange, reso_XP);
		ptge_mom_reso_XP->			SetPointError(ivalid, mombinwidth/2, dreso_XP);
		ivalid++;
	}

	// Drift distance
	for(int i=0; i<nmombins; i++){
		for(int j=0; j<nddbins; j++){
			int nentries_here =			arr_momdd_fph1f_XP[i][j]->GetEntries();
			if(nentries_here < 50) continue;
			fptf1_XP =					Fit1Gauss(arr_momdd_fph1f_XP[i][j]);

			float mean_XP =				fptf1_XP->GetParameter(1);
			float dmean_XP =			fptf1_XP->GetParError(1);

			float std_XP =				fptf1_XP->GetParameter(2);
			float dstd_XP =				fptf1_XP->GetParError(2);

			float reso_XP =				std_XP/mean_XP*100;
			float dreso_XP =			GetResoError(fptf1_XP);

			fph2f_momdd_reso_XP->		Fill(j*ddbinwidth, i*mombinwidth-momrange, reso_XP);
		}
	}

	// 	// Bethe-Bloch curve -------------------------------------------------------------------------------------------------------------------
	// float MeV =						5.9e-3/(224*1493/183);
	// const char *formula =			"[1]/pow(x/sqrt(x*x+[0]*[0]),[4])	*( [2] - pow(x/sqrt(x*x+[0]*[0]),[4]) - log([3]+ pow(x/[0], [5])) )/[6]";
	// TF1 *BBcurve =					new TF1(Form("dEdx_cosmics"), formula, 1, 50e3, "");
	// BBcurve->						SetParameters(105.658, 0.186543e-3, 5.382656, 0.004234, 2.028548, -0.994807, MeV);
}
