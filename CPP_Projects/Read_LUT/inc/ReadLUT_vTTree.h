#ifndef lut_H
#define lut_H

#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

#include "Util.h"

class LUT {
public:
   /* Constructor */
   LUT();
   virtual ~LUT();

   float ratio(const float &phi, const float &d, const float &RC, const float &z);

private:
   std::string m_file;
   float m_value[200][100][41][21];

   float pad_diag = sqrt(pow(11.28, 2) + pow(10.19, 2));

   int nphi = 100;
   int nd = 100;
   int nRC = 41;
   int nz = 21;

   float step_phi = (90 - 2e-6) / (nphi - 1);
   float step_d = pad_diag / 2 / (nd - 1);
   float step_RC = 200. / (nRC - 1);
   float step_z = 1000. / (nz - 1);

   std::vector<double> v_iphi = linspace(1e-6, 90 - 1e-6, nphi);
   std::vector<double> v_id = linspace(0., pad_diag / 2, nd);
   std::vector<double> v_iRC = linspace(50., 250., nRC);
   std::vector<double> v_iz = linspace(0., 1000., nz);

   void Load();

   //------------------------------Data Members-----------------------//

   TFile *p_file;
   TTree *p_tree;

   // Leaves
   float m_weight;
   float m_d;
   float m_phi;
   float m_RC;
   float m_z;
};

#endif