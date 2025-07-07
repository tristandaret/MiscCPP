#ifndef DRAWERS_H
#define DRAWERS_H

#include <TCanvas.h>
#include <TH1.h>
#include <TLegend.h>

#include "TrackModel2.h"
#include "TrackModel.h"
#include "PointModel.h"

void DrawTrackModel2(double phi=0, double impact=0, double drift=0, double RC=120);

void DrawPointModel();

#endif // DRAWERS_H