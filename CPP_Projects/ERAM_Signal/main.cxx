#include "PointModel.h"
#include "SetStyle.h"

#include <TROOT.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>

int main(int argc, char const *argv[])
{
	// Set Style
	float Q0 = 						10;
	int lw = 						4;
	TCanvas* pCanvas =				new TCanvas("pCanvas", "pCanvas", 1600, 1200);
	TStyle* ptstyle =				SetMyStyle();
	gROOT->							SetStyle(ptstyle->GetName());
	gStyle->						SetPadRightMargin(0.05);
	gStyle->						SetPadTopMargin(0.05);
	gPad->							UseCurrentStyle();

	int nbins = 					5e3;
	int peaking = 					std::stod(argv[1]);
	PointModel *p_pointmodel = 		new PointModel(peaking, 11.28/2, 10.19/2, 11.28/2, 10.19/2, 120, Q0);

	// TH1F *ph1f_currentmano =		new TH1F("ph1f_currentmano", "Current", nbins, 0, nbins);
	// TH1F *ph1f_dETFdtmano = 		new TH1F("ph1f_dETFdtmano", "Derivative of electronics", nbins, 0, nbins);
	// TH1F *ph1f_convorediv =			new TH1F("ph1f_convorediv", "Derivative of Q #otimes ETF", nbins, 0, nbins);
	// for(int i = 0; i < nbins-1; i++){
	// 	ph1f_currentmano->			Fill(i, 1e3*(p_pointmodel->ptf1_Charge->Eval(i+1)-p_pointmodel->ptf1_Charge->Eval(i)));
	// 	ph1f_dETFdtmano->			Fill(i, p_pointmodel->ptf1_ETF->Eval(i+1)-p_pointmodel->ptf1_ETF->Eval(i));
	// 	ph1f_convorediv->			Fill(i, p_pointmodel->ptf1_Signal3->Eval(i+1)-p_pointmodel->ptf1_Signal3->Eval(i));
	// }

	// TH1F *ph1f_convomanoIETF  = new TH1F("ph1f_convomanoIETF", ";time (ns);ADC counts", 5e3, 0, 5e3);
	// TH1F *ph1f_convomanoQdETFdt = new TH1F("ph1f_convomanoQdETFdt", ";time (ns);ADC counts", 5e3, 0, 5e3);
    // for(int i = 0; i < nbins; i++){
	// 	Double_t ix = i;
	// 	ph1f_convomanoIETF->Fill(i, Q0/nbins*p_pointmodel->ETF(&ix, 0));
    //     for(int j = 0; j < nbins; j++){
    //         Double_t jx = j;
    //         Double_t kx = i-j;
    //         ph1f_convomanoIETF->Fill(i, p_pointmodel->ETF(&jx, 0)*p_pointmodel->Current(&kx, 0));
	// 		ph1f_convomanoQdETFdt->Fill(i, p_pointmodel->Charge(&jx, 0)*p_pointmodel->dETFdt(&kx, 0));
    //     }
    // }

	p_pointmodel->ptf1_Charge->SetNpx(nbins);
	p_pointmodel->ptf1_Charge->SetLineWidth(lw);
	p_pointmodel->ptf1_Charge->SetLineColor(kBlue);
	p_pointmodel->ptf1_CurrentScaled->SetNpx(nbins);
	p_pointmodel->ptf1_CurrentScaled->SetLineWidth(lw);
	p_pointmodel->ptf1_CurrentScaled->SetLineColor(kRed);
	p_pointmodel->ptf1_ETF->SetNpx(nbins);
	p_pointmodel->ptf1_ETF->SetLineWidth(lw);
	p_pointmodel->ptf1_ETF->SetLineColor(kRed);
	p_pointmodel->ptf1_dETFdt->SetNpx(nbins);
	p_pointmodel->ptf1_dETFdt->SetLineWidth(lw);
	p_pointmodel->ptf1_dETFdt->SetLineColor(kBlue);
	gStyle->SetOptStat(0);

	p_pointmodel->ptf1_Charge->Draw();
	pCanvas->SaveAs(Form("PointSignal_%d.pdf(", peaking));

	pCanvas->Clear();
	p_pointmodel->ptf1_CurrentScaled->Draw();
	// ph1f_currentmano->SetLineStyle(2);
	// ph1f_currentmano->SetLineWidth(lw);
	// ph1f_currentmano->Draw("hist same");
	pCanvas->SaveAs(Form("PointSignal_%d.pdf", peaking));

	pCanvas->Clear();
	p_pointmodel->ptf1_ETF->Draw();
	pCanvas->SaveAs(Form("PointSignal_%d.pdf", peaking));

	pCanvas->Clear();
	p_pointmodel->ptf1_dETFdt->Draw();
	// ph1f_dETFdtmano->SetLineStyle(2);
	// ph1f_dETFdtmano->SetLineWidth(lw);
	// ph1f_dETFdtmano->Draw("hist same");
	pCanvas->SaveAs(Form("PointSignal_%d.pdf", peaking));

	pCanvas->Clear();
	p_pointmodel->ptf1_Signal->SetNpx(nbins);
	p_pointmodel->ptf1_Signal->SetLineWidth(lw);
	p_pointmodel->ptf1_Signal->SetLineColor(kRed);
	p_pointmodel->ptf1_Signal2->SetNpx(nbins);
	p_pointmodel->ptf1_Signal2->SetLineWidth(lw);
	p_pointmodel->ptf1_Signal2->SetLineStyle(2);
	p_pointmodel->ptf1_Signal2->SetLineColor(kBlue);
	p_pointmodel->ptf1_SignalnoDirac->SetNpx(nbins);
	p_pointmodel->ptf1_SignalnoDirac->SetLineWidth(lw);
	p_pointmodel->ptf1_SignalnoDirac->SetLineColor(kGray+2);
	p_pointmodel->ptf1_Signal->SetMaximum(1.1*std::max(p_pointmodel->ptf1_Signal->GetMaximum(), p_pointmodel->ptf1_Signal2->GetMaximum()));
	p_pointmodel->ptf1_Signal->SetMinimum(1.1*std::min(p_pointmodel->ptf1_Signal->GetMinimum(), p_pointmodel->ptf1_Signal2->GetMinimum()));
	p_pointmodel->ptf1_Signal->DrawClone();
	p_pointmodel->ptf1_Signal2->Draw("same");
	p_pointmodel->ptf1_SignalnoDirac->Draw("same");
	TLegend *plegend = new TLegend(0.5, 0.6, 0.9, 0.9);
	p_pointmodel->ptf1_Signal->SetLineStyle(2);
	plegend->AddEntry(p_pointmodel->ptf1_Signal, "I^{H} #otimes T^{H}", "l");
	plegend->AddEntry(p_pointmodel->ptf1_SignalnoDirac, "I #otimes T^{H}", "l");
	plegend->AddEntry(p_pointmodel->ptf1_Signal2, "Q^{H} #otimes #partial_{t}T^{H}", "l");
	plegend->Draw();
	pCanvas->SaveAs(Form("PointSignal_%d.pdf", peaking));

	pCanvas->Clear();
	p_pointmodel->ptf1_Signal->SetNpx(nbins);
	p_pointmodel->ptf1_Signal->SetLineWidth(lw);
	p_pointmodel->ptf1_Signal->SetLineStyle(1);
	p_pointmodel->ptf1_Signal->SetLineColor(kRed);
	p_pointmodel->ptf1_SignalnoDirac->SetNpx(nbins);
	p_pointmodel->ptf1_SignalnoDirac->SetLineWidth(lw);
	p_pointmodel->ptf1_SignalnoDirac->SetLineStyle(2);
	p_pointmodel->ptf1_SignalnoDirac->SetLineColor(kGray+2);
	p_pointmodel->ptf1_ETFscaled->SetNpx(nbins);
	p_pointmodel->ptf1_ETFscaled->SetLineWidth(lw);
	p_pointmodel->ptf1_ETFscaled->SetLineStyle(2);
	p_pointmodel->ptf1_ETFscaled->SetLineColor(kBlue);
	p_pointmodel->ptf1_Signal->SetMaximum(1.1*std::max(p_pointmodel->ptf1_Signal->GetMaximum(), p_pointmodel->ptf1_Signal2->GetMaximum()));
	p_pointmodel->ptf1_Signal->SetMinimum(1.1*std::min(p_pointmodel->ptf1_Signal->GetMinimum(), p_pointmodel->ptf1_Signal2->GetMinimum()));
	p_pointmodel->ptf1_Signal->Draw();
	p_pointmodel->ptf1_SignalnoDirac->Draw("same");
	p_pointmodel->ptf1_ETFscaled->Draw("same");
	delete plegend;
	plegend = new TLegend(0.45, 0.6, 0.9, 0.9);
	plegend->SetTextSize(0.075);
	plegend->AddEntry(p_pointmodel->ptf1_Signal, "Complete signal", "l");
	plegend->AddEntry(p_pointmodel->ptf1_ETFscaled, "Pulse signal", "l");
	plegend->AddEntry(p_pointmodel->ptf1_SignalnoDirac, "Escape signal", "l");
	plegend->Draw();
	pCanvas->SaveAs(Form("PointSignal_%d.pdf)", peaking));
	// pCanvas->SaveAs(Form("PointSignal_%d.pdf", peaking));

	// pCanvas->Clear();
	// ph1f_convomanoIETF->SetLineColor(kRed);
	// ph1f_convomanoIETF->SetLineWidth(lw+3);
	// ph1f_convomanoIETF->SetLineStyle(2);
	// ph1f_convomanoQdETFdt->SetLineColor(kBlue);
	// ph1f_convomanoQdETFdt->SetLineWidth(lw);
	// ph1f_convomanoQdETFdt->SetLineStyle(2);
	// ph1f_convorediv->SetLineColor(kGreen);
	// ph1f_convorediv->SetLineWidth(lw);
	// ph1f_convomanoIETF->SetMinimum(1.1*std::min(ph1f_convomanoIETF->GetMinimum(), ph1f_convomanoQdETFdt->GetMinimum()));
	// ph1f_convomanoIETF->SetMaximum(1.1*std::max(ph1f_convomanoIETF->GetMaximum(), ph1f_convomanoQdETFdt->GetMaximum()));
	// ph1f_convomanoIETF->Draw("hist");
	// ph1f_convorediv->Draw("hist same");
	// ph1f_convomanoQdETFdt->Draw("hist same");
	// plegend->GetListOfPrimitives()->RemoveAt(1);
	// plegend->AddEntry(ph1f_convorediv, "d_{t}(Q^{H} #otimes T^{H})", "l");
	// plegend->Draw();
	// pCanvas->SaveAs(Form("PointSignal_%d.pdf)", peaking));
}