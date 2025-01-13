#include <cmath>

double Qpad(double xpad, double ypad, double xcharge, double ycharge, double RC, double t, double Qdrop)
{
   /* 	Distances are in mm
      (xpad, ypad) is the center of the pad considered
      (xcharge, ycharge) is the position of the charge
      Qdrop and Qpad are in whatever you want
      RC is in ns/mm2 */

   // Avoid division by zero from sigma
   if (t == 0)
      t = 1e-6;

   // Defining the border of the pad considered
   constexpr double xwidth = 11.28;
   constexpr double ywidth = 10.19;
   double xleft = xpad - xwidth / 2;
   double xright = xpad + xwidth / 2;
   double ylow = ypad - ywidth / 2;
   double yhigh = ypad + ywidth / 2;

   // Computation
   double sigma = std::sqrt(2 * t / RC); // Does not include transverse diffusion
   double erfx =
      std::erf((xright - xcharge) / (sigma * std::sqrt(2))) - std::erf((xleft - xcharge) / (sigma * std::sqrt(2)));
   double erfy =
      std::erf((yhigh - ycharge) / (sigma * std::sqrt(2))) - std::erf((ylow - ycharge) / (sigma * std::sqrt(2)));

   return Qdrop / 4 * erfx * erfy;
}