#include "PointModel.h"


// Constructor
PointModel::PointModel(float peakingTime, float fxpad, float fypad, float xcharge, float ycharge, float fRC, float fQdrop){
    SetPeakingTime(peakingTime);
    SetPadParam(fxpad, fypad, fRC);
    SetChargeParam(xcharge, ycharge, fQdrop);

    // Defining the border of the pad considered
    fxleft =                     fxpad - fSXWIDTH/2;
    fxright =                    fxpad + fSXWIDTH/2;
    fylow =                      fypad - fSYWIDTH/2;
    fyhigh =                     fypad + fSYWIDTH/2;

    ptf1_ETF =			new TF1(";time (ns);ADC counts",       this, &PointModel::ETF,     0, ftmax, 0);
    ptf1_ETFscaled =	new TF1(";time (ns);ADC counts",       this, &PointModel::ETFscaled,     0, ftmax, 0);
    ptf1_dETFdt =		new TF1(";time (ns); ADC counts/ns",    this, &PointModel::dETFdt,  0, ftmax, 0);
    ptf1_Charge =		new TF1(";time (ns); Charge (fC)",     this, &PointModel::Charge,  0, ftmax, 0);
    ptf1_Current =		new TF1(";time (ns);Current (#muA)",   this, &PointModel::Current, 0, ftmax, 0);
    ptf1_CurrentScaled =new TF1(";time (ns);Current (nA)",     this, &PointModel::CurrentScaled, 0, ftmax, 0);
    ptf1_CurrentnoDirac =	new TF1(";time (ns);Current (#muA)",   this, &PointModel::CurrentnoDirac, 0, ftmax, 0);

    ptf1_Convolution =	new TF1Convolution(ptf1_Current, ptf1_ETF, 0, ftmax, true);
    ptf1_Convolution->	SetRange(-ftmax, ftmax);
    ptf1_Convolution->	SetNofPointsFFT(fnconvolpoints);
	ptf1_Signal =		new TF1(";time (ns);ADC counts", *ptf1_Convolution, 0, ftmax, ptf1_Convolution->GetNpar());

    ptf1_ConvolutionnoDirac =	new TF1Convolution(ptf1_CurrentnoDirac, ptf1_ETF, 0, ftmax, true);
    ptf1_ConvolutionnoDirac->	SetRange(-ftmax, ftmax);
    ptf1_ConvolutionnoDirac->	SetNofPointsFFT(fnconvolpoints);
	ptf1_SignalnoDirac =		new TF1("", *ptf1_ConvolutionnoDirac, 0, ftmax, ptf1_Convolution->GetNpar());

    ptf1_Convolution2 =	new TF1Convolution(ptf1_Charge, ptf1_dETFdt, 0, ftmax, true);
    ptf1_Convolution2->	SetRange(-ftmax, ftmax);
    ptf1_Signal2 =		new TF1("", *ptf1_Convolution2, 0, ftmax, ptf1_Convolution2->GetNpar());

    ptf1_Convolution3 =	new TF1Convolution(ptf1_Charge, ptf1_ETF, 0, ftmax, true);
    ptf1_Convolution3->	SetRange(-ftmax, ftmax);
    ptf1_Signal3 =		new TF1("", *ptf1_Convolution3, 0, ftmax, ptf1_Convolution3->GetNpar());
}


// Destructor
PointModel::~PointModel(){
    delete ptf1_ETF;
    delete ptf1_ETFscaled;
    delete ptf1_dETFdt;
    delete ptf1_Charge;
    delete ptf1_Current;
    delete ptf1_CurrentScaled;
    delete ptf1_CurrentnoDirac;
    delete ptf1_Convolution;
    delete ptf1_ConvolutionnoDirac;
    delete ptf1_Convolution2;
    delete ptf1_Convolution3;
    delete ptf1_Signal;
}



// Settings
void PointModel::SetPeakingTime(float peakingTime){
    fpeakingTime =              peakingTime;
    ws =                        2/fpeakingTime;
    expfactor =                 ws/(2*Q);
    arg =                       ws/2*std::sqrt(4-1/std::pow(Q,2));
}




// Eletronic transfer function
Double_t PointModel::ETF(Double_t *x, Double_t *par){
    Double_t t = x[0];

	if(t <= 0) return 0;
	else return 4096/120*(std::exp(-ws*t) + std::exp(-expfactor*t) * (sinfactor*std::sin(arg*t) - std::cos(arg*t)))/fSETFMAX;
}
Double_t PointModel::ETFscaled(Double_t *x, Double_t *par){
    Double_t t = x[0];

	if(t <= 0) return 0;
	else return 10*4096/120*(std::exp(-ws*t) + std::exp(-expfactor*t) * (sinfactor*std::sin(arg*t) - std::cos(arg*t)))/fSETFMAX;
}

Double_t PointModel::dETFdt(Double_t *x, Double_t *par){
    Double_t t = x[0];

    if(t <= 0) return 0;
    else return 4096/120*(-ws*std::exp(-ws*t) + std::exp(-expfactor*t) * (expfactor*(std::cos(arg*t) - sinfactor*std::sin(arg*t)) + arg*(std::sin(arg*t) + sinfactor*std::cos(arg*t))))/fSETFMAX;
}



// Charge function
Double_t PointModel::Charge(Double_t *x, Double_t *par){
    Double_t t = x[0];
    
    if(t <= 0) return 0;

    float sigma =           std::sqrt(2*t/fRC);
    auto erfRC =            [sigma](float edge, float charge){return std::erf((edge - charge) / (std::sqrt(2)*sigma));};
    
    return fQdrop/4 * (erfRC(fyhigh,  fycharge) - erfRC(fylow,  fycharge)) * (erfRC(fxright, fxcharge) - erfRC(fxleft, fxcharge));
}



// Current function
Double_t PointModel::Current(Double_t *x, Double_t *par){
    Double_t t = x[0];

    if(t<0) return 0;
    if(t<=1/fnconvolpoints) return fQdrop/(2*ftmax)*fnconvolpoints; // normalized to mimic a dirac delta function

    float sigma =   std::sqrt(2*t / fRC);
    float norm =    fQdrop / (4*std::sqrt(2*TMath::Pi())*sigma*t);
    auto erfRC =    [sigma](float edge, float charge){return std::erf((edge - charge) / (std::sqrt(2)*sigma));};
    auto expRC =    [sigma, t](float edge, float charge){return std::exp(-std::pow(edge - charge, 2) / (2 * std::pow(sigma, 2))) * (edge - charge);};

    return  norm * ( (erfRC(fyhigh,  fycharge) - erfRC(fylow,  fycharge)) * (expRC(fxleft, fxcharge) - expRC(fxright, fxcharge)) 
                    +(erfRC(fxright, fxcharge) - erfRC(fxleft, fxcharge)) * (expRC(fylow,  fycharge) - expRC(fyhigh,  fycharge)));
}

Double_t PointModel::CurrentScaled(Double_t *x, Double_t *par){
    Double_t t = x[0];

    if(t<0) return 0;
    if(t<=1/fnconvolpoints) return fQdrop/(2*ftmax)*fnconvolpoints; // normalized to mimic a dirac delta function

    float sigma =   std::sqrt(2*t / fRC);
    float norm =    fQdrop / (4*std::sqrt(2*TMath::Pi())*sigma*t);
    auto erfRC =    [sigma](float edge, float charge){return std::erf((edge - charge) / (std::sqrt(2)*sigma));};
    auto expRC =    [sigma, t](float edge, float charge){return std::exp(-std::pow(edge - charge, 2) / (2 * std::pow(sigma, 2))) * (edge - charge);};

    return  1e3*norm * ( (erfRC(fyhigh,  fycharge) - erfRC(fylow,  fycharge)) * (expRC(fxleft, fxcharge) - expRC(fxright, fxcharge)) 
                    +(erfRC(fxright, fxcharge) - erfRC(fxleft, fxcharge)) * (expRC(fylow,  fycharge) - expRC(fyhigh,  fycharge)));
}

Double_t PointModel::CurrentnoDirac(Double_t *x, Double_t *par){
    Double_t t = x[0];

    if(t<=0) return 0;

    float sigma =   std::sqrt(2*t / fRC);
    float norm =    fQdrop / (4*std::sqrt(2*TMath::Pi())*sigma*t);
    auto erfRC =    [sigma](float edge, float charge){return std::erf((edge - charge) / (std::sqrt(2)*sigma));};
    auto expRC =    [sigma, t](float edge, float charge){return std::exp(-std::pow(edge - charge, 2) / (2 * std::pow(sigma, 2))) * (edge - charge);};

    return  norm * ( (erfRC(fyhigh,  fycharge) - erfRC(fylow,  fycharge)) * (expRC(fxleft, fxcharge) - expRC(fxright, fxcharge)) 
                    +(erfRC(fxright, fxcharge) - erfRC(fxleft, fxcharge)) * (expRC(fylow,  fycharge) - expRC(fyhigh,  fycharge)));
}