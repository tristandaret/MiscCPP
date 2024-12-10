#include "TCanvas.h"
#include "TPaveText.h"

int main() {

    std::string options = "NDC ARC";
    TPaveText pavetext2(0.15, 0.3, 0.9, 0.7, options.c_str());
	pavetext2.SetCornerRadius(0.1);
    pavetext2.SetTextColor(kGray);
    pavetext2.SetFillColor(kBlue);
    pavetext2.SetLineColor(kRed);
    pavetext2.SetBorderSize(1);
    pavetext2.SetLineWidth(10);
    pavetext2.AddText(("options " + options).c_str());

    TCanvas canvas("canvas", "MWE canvas", 1600, 1200);
    pavetext2.Draw();

    canvas.SaveAs("pavetext_bug_arc.png");
}