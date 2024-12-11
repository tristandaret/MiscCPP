#ifndef GEOMFUNCTIONS_H
#define GEOMFUNCTIONS_H

#include <cmath>

inline float GetX(const float &y, const float &phi_rad, const float &d) { return (y - (d-sin(phi_rad)*11.28/2+cos(phi_rad)*10.19/2)/cos(phi_rad))/tan(phi_rad);}
inline float GetY(const float &x, const float &phi_rad, const float &d) { return tan(phi_rad)*x + (d-sin(phi_rad)*11.28/2+cos(phi_rad)*10.19/2)/cos(phi_rad);}

float ComputeLength(const float &phi, const float &d);

#endif // GEOMFUNCTIONS_H