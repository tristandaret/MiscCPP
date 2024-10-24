#include "ReadLUT_vTTree.h"
#include "Util.h"
#include <iomanip>
#include <iostream>
#include <sstream>

// Default constructor
LUT::LUT()
{
  m_file = Form("~/Documents/Code/Python/LUT/LUT_Dt350_PT412_nphi%i_nd%i_nRC%i_nZ%i.root", nphi, nd, nRC, nz);

  Load();
}

// Destructor
LUT::~LUT()
{
  p_file->Close();
  delete p_file;
  for(int i=0;i<nphi;i++) for(int j=0;j<nd;j++) for(int k=0;k<nRC;k++) for(int l=0;l<nz;l++) m_value[i][j][k][l] = 0;
}


// Load the ERAM maps
void LUT::Load()
{  
  p_file                      = TFile::Open(m_file.c_str(),"READ");
  p_tree                      = (TTree*) p_file->Get("outTree");
  p_tree->                      SetBranchAddress("weight",        &m_weight);
  p_tree->                      SetBranchAddress("angle",         &m_phi);
  p_tree->                      SetBranchAddress("impact_param",  &m_d);
  p_tree->                      SetBranchAddress("RC",            &m_RC);
  p_tree->                      SetBranchAddress("drift_dist",    &m_z);

  
  // Initializing all values to 0
  for(int i=0;i<nphi;i++) for(int j=0;j<nd;j++) for(int k=0;k<nRC;k++) for(int l=0;l<nz;l++) m_value[i][j][k][l] = 0;

  int nentries = p_tree->GetEntries();
  std::cout << "Number of entries: " << nentries << std::endl;
  int id, iphi, iRC, iz;
  for (int i=0; i<nentries;i++) {
    p_tree->                      GetEntry(i);
    iphi                        = (int)std::round((m_phi-1e-6)/step_phi);
    id                          = (int)std::round(m_d/step_d);
    iRC                         = (int)std::round((m_RC-50)/step_RC);
    iz                          = (int)std::round(m_z/step_z);
    m_value[iphi][id][iRC][iz]  = m_weight;
  }
}


float LUT::ratio(const float& phi, const float& d, const float& RC, const float& z)
{
  // Determine the indices of the 8 points surrounding the point of interest
  float iphi     = (phi-1e-6)/step_phi;
  float iphi_min = std::min(std::floor((phi-1e-6)/step_phi), (double)nphi-1);
  float iphi_max = std::max(std::ceil((phi-1e-6)/step_phi), 0.0);
  float id       = d/step_d;
  float id_min   = std::min(std::floor(d/step_d), (float)nd-1);
  float id_max   = std::max(std::ceil(d/step_d), (float)0);
  float iRC      = (RC-50)/step_RC;
  float iRC_min  = std::min(std::floor((RC-50)/step_RC), (float)(nRC-1));
  float iRC_max  = std::max(std::ceil((RC-50)/step_RC), (float)0);
  float iz       = z/step_z;
  float iz_min   = std::min(std::floor(z/step_z), (float)(nz-1));
  float iz_max   = std::max(std::ceil(z/step_z), (float)0);

  // Check if the point of interest is outside the range of the LUT
  if(iphi_min > nphi-1) iphi_min  = nphi-1;
  if(iphi_max < 0)      iphi_max  = 0;
  if(id_min > nd-1)     id_min    = nd-1;
  if(id_max < 0)        id_max    = 0;
  if(iRC_min > nRC-1)   iRC_min   = nRC-1;
  if(iRC_max < 0)       iRC_max   = 0;
  if(iz_min > nz-1)     iz_min    = nz-1;
  if(iz_max < 0)        iz_max    = 0;

  // Determine the weights for the 8 points
  float w_phi, w_d, w_RC, w_z;
  if(iphi_min == iphi_max)  w_phi = 1;
  else                      w_phi = 1 - (iphi - iphi_min)/(iphi_max - iphi_min);
  if(id_min == id_max)      w_d   = 1;
  else                      w_d   = 1 - (id - id_min)/(id_max - id_min);
  if(iRC_min == iRC_max)    w_RC  = 1;
  else                      w_RC  = 1 - (iRC - iRC_min)/(iRC_max - iRC_min);
  if(iz_min == iz_max)      w_z   = 1;
  else                      w_z   = 1 - (iz - iz_min)/(iz_max - iz_min);

  float interpolated_value = 0;
for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) for (int k = 0; k < 2; ++k) for (int l = 0; l < 2; ++l) 
  interpolated_value += m_value[(int)iphi_min + i][(int)id_min + j][(int)iRC_min + k][(int)iz_min + l] *
  (i == 0 ? w_phi : (1 - w_phi)) *
  (j == 0 ? w_d : (1 - w_d)) *
  (k == 0 ? w_RC : (1 - w_RC)) *
  (l == 0 ? w_z : (1 - w_z));

  return interpolated_value;
}