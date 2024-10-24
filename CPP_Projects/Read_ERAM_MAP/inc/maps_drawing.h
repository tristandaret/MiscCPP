#ifndef MAPS_DRAWING_H
#define MAPS_DRAWING_H

#include "LUTs.h"
#include "TCanvas.h"
#include "TLegend.h"

class Drawing{
    public:
        Drawing();
        virtual ~Drawing();
        void maps_pdf();
        void maps_csv();
        void Distributions();

    private:
        Reconstruction::ERAMMaps *pERAMMaps;
        std::string OutputFile;
        std::string OutputFile_Beg;
        std::string OutputFile_End;
        TCanvas* pTCanvas;
        TLegend* pTLegend;

};

#endif