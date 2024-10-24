#include "ReadLUT_vTTree.h"
#include "Util.h"

#include <algorithm>
#include <iostream>

#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>

int main(){
  LUT* p_lut = new LUT();
  float phi = 2;
  float d   = 4.978;
  float RC  = 150;
  float z   = 50;
  std::cout << "ratio(" << phi << ", " << d << ", " << RC << ", " << z << ") = " << p_lut->ratio(phi, d, RC, z) << std::endl;

  return 0;
}