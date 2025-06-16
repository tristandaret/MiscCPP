#ifndef MISC_T_H
#define MISC_T_H

#include "Util.h"
#include <type_traits> // for sign template

/* GENERAL FUNCTIONS*/
// Find if a value is in a vector
bool is_in(std::vector<double> v, double val);

// Equivalent of numpy linspace (npoints uniformly spaced between start and end)
std::vector<double> linspace(double start, double end, int numPoints);

// Sign template
template <typename T>
int sign(T x)
{
   return (x >= 0) ? 1 : (x < 0 ? -1 : 0);
}

// Add vectors
template <typename T>
std::vector<T>& operator+=(std::vector<T>& a, const std::vector<T>& b) {
    if (a.size() != b.size()) {
        throw std::runtime_error("Vector sizes don't match in operator+=");
    }
    for (size_t i = 0; i < a.size(); ++i) {
        a[i] += b[i];
    }
    return a;
}

// Substract vectors
template <typename T>
std::vector<T>& operator-=(std::vector<T>& a, const std::vector<T>& b) {
    if (a.size() != b.size()) {
        throw std::runtime_error("Vector sizes don't match in operator-=");
    }
    for (size_t i = 0; i < a.size(); ++i) {
        a[i] -= b[i];
    }
    return a;
}

/* GENERAL MATH*/
// Heaviside
float Heaviside(const float &t, const float &t_0);

// Mean of vector
float mean(const std::vector<float> &values);

// Fit a gaussian but within +/- 2sigma to fit only the peak correctly
TF1 *Fit1Gauss(TH1 *h1F);
TF1 *Fit1Gauss(TH1 *h1F, const float &range);

/* ROOT */
double GetSeparation(const float &mean1, const float &std1, const float &mean2, const float &std2);
double GetSeparation(const TF1 *tf1_1, const TF1 *tf1_2);
double GetSeparationError(const float &mean1, const float &std1, const float &dmean1, const float &dstd1,
                          const float &mean2, const float &std2, const float &dmean2, const float &dstd2);
double GetSeparationError(const TF1 *tf1_1, const TF1 *tf1_2);

// Formula to get resolution error using a TF1
double GetResoError(TF1 *tf1);
double GetResoError(TF1 *tf1, const int &mu, const int &sigma);

// Write resolution value for a 1-gaussian fit
void PrintResolution(TH1 *th1, TCanvas *pCanvas);
void PrintResolution(TH1 *th1, TCanvas *pCanvas, float NDCx, float NDCy, Color_t color, const std::string &title);
void PrintResolution(TH1 *th1, TCanvas *pCanvas, float NDCx, float NDCy, const float &xwidth, const float &ywidth,
                     const std::string &anchor, Color_t color, const std::string &title);

// Draw TH1
void DrawTH1(const std::string &OutDir, TH1 *h1);

#endif // MISC_T_H