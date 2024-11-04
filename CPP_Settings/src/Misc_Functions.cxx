#include "Misc_Functions.h"



// General functions //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Find if a value is in a vector
bool is_in(std::vector<double> v, double val){
	if(std::find(v.begin(), v.end(), val) != v.end()) return true;
	else return false;
}

// Equivalent of numpy linspace (npoints uniformly spaced between start and end)
std::vector<double> linspace(double start, double end, int numPoints) {
	std::vector<double> result(numPoints);
	double step = (end - start) / (numPoints - 1);
	
	for (int i = 0; i < numPoints; ++i) {
		result[i] = start + i * step;
	}
	return result;
}



// General Math //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Mean of vector
float mean(const std::vector<float>& values) {
	float sum = 0.0;
	for (const float& value : values) {
		sum += value;
	}

	return sum / values.size();
}



// Heaviside
float Heaviside(const float& t, const float& t_0){
	if(t_0 >	t)	return 0.0 ;
	if(t_0 == t)	return 0.5 ;
	if(t_0 <	t)	return 1.0 ;
	else			return -999.0 ;
}



TF1* Fit1Gauss(TH1* h1F) {return Fit1Gauss(h1F, 2);}
TF1* Fit1Gauss(TH1* h1F, const float &range){	
	TF1* gausn	= new TF1("gausn", "gausn") ;
	float mean	= h1F->GetBinCenter(h1F->GetMaximumBin()) ;
	float std	 = h1F->GetRMS() ;
	float max	 = mean + 1.8*std ;
	float min	 = mean - 1.8*std ;
	if (h1F->Fit("gausn","RQ","0", min, max) != 0) return gausn;

	TF1* g1		= h1F->GetFunction("gausn") ;
	mean			= g1->GetParameter(1) ;
	std			= g1->GetParameter(2) ;
	max			= mean + 1.7*std ;
	min			= mean - 1.7*std ;
	h1F->Fit("gausn","RQ","0", min, max) ;

	TF1* g2		= h1F->GetFunction("gausn") ;
	mean			= g2->GetParameter(1) ;
	std			= g2->GetParameter(2) ;
	max			= mean + 1.6*std ;
	min			= mean - 1.6*std ;
	h1F->Fit("gausn","RQ","0", min, max) ;

	TF1* g3		= h1F->GetFunction("gausn") ;
	mean			= g3->GetParameter(1) ;
	std			= g3->GetParameter(2) ;
	max			= mean + range*std ;
	min			= mean - range*std ;
	h1F->Fit("gausn","RQ","0", min, max) ;

	TF1* tf1 = 			h1F->GetFunction("gausn") ;
	delete gausn ;

	return tf1 ;
} 


// ROOT /////////////////////////////////////////////////////////////////////////////////////////////////////////

// Draw TH1
void DrawTH1(const std::string& OutDir, TH1* h1){
	TCanvas* pCanTH1				= new TCanvas("pCanTH1", "pCanTH1", 1800, 1200);
	pCanTH1->cd();
	h1->SetLineWidth(4);
	h1->Draw("hist");
	pCanTH1->SaveAs(OutDir.c_str());
	delete pCanTH1;
}




double GetResoError(TF1* tf1){
	return GetResoError(tf1, 1, 2) ;
}
double GetResoError(TF1* tf1, const int& mu, const int& sigma){
	double mean	 = tf1->GetParameter(mu) ;
	double std		= tf1->GetParameter(sigma) ;
	double dmean	= tf1->GetParError(mu) ;
	double dstd	 = tf1->GetParError(sigma) ;
	return (std*dmean/(mean*mean) + dstd/mean) * 100 ;
}



void PrintResolution(TH1* th1, TCanvas* pCanvas){
	return PrintResolution(th1, pCanvas, 0.05, 0.7, kBlack, " ");
}
void PrintResolution(TH1* th1, TCanvas* pCanvas, float NDCx, float NDCy, Color_t color, const std::string& title){
	return PrintResolution(th1, pCanvas, NDCx, NDCy, 0.3, color, title);
}
void PrintResolution(TH1* th1, TCanvas* pCanvas, float NDCx, float NDCy, const float &size, Color_t color, const std::string& title){
	TF1* tf1 =					th1->GetFunction("gausn");
	if (!tf1){
		tf1 =	 				Fit1Gauss(th1);
		if (!tf1) {
			std::cerr << "Error: No fit found for " << th1->GetName() << std::endl;
			return;
		}
	}

	double	xMax =			pCanvas->GetUxmax();
	double	yMax =			pCanvas->GetUymax();
	TPaveText *pPaveText =	new TPaveText(NDCx, NDCy, NDCx+size, NDCy+size, "NDC");
	pPaveText->				SetFillStyle(0);
	pPaveText->				SetTextAlign(12);
	pPaveText->				SetLineColor(color);
	pPaveText->				SetTextColor(color);
	pPaveText->				SetShadowColor(0);
	pPaveText->				SetLineWidth(1.5);

	float mu	=				tf1->GetParameter(1);
	float dmu	=			tf1->GetParError(1);
	float sigma	=			tf1->GetParameter(2);
	float dsigma =			tf1->GetParError(2);
	float reso	= 			tf1->GetParameter(2)/tf1->GetParameter(1) * 100;
	float dreso	= 			GetResoError(tf1);

	pPaveText->				AddText(Form("%s (%d entries)", title.c_str(), (int)th1->GetEntries()));
	// pPaveText->				AddText(Form("%s", title.c_str()));
	pPaveText->				AddText(Form("#frac{#sigma}{#mu}	= %.2f #pm %.2f %%", reso, dreso));
	pPaveText->				AddText(Form("#mu	= %.1f #pm %.1f", mu, dmu));
	pPaveText->				AddText(Form("#sigma	= %.1f #pm %.1f", sigma, dsigma));
	pPaveText->				GetLine(0)->SetTextFont(22);
	pPaveText->				DrawClone();
	delete pPaveText;
}