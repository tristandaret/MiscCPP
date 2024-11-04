#ifndef DRAW_H
#define DRAW_H

#include "Util.h"
#include "Process.h"

class Draw
{
    public:
        Draw();
        ~Draw();

        // Member functions
        void Run(const Process &process);
        void CompareRuns(const std::vector<Process*> &v_processes, const std::string &compfolderpath, const std::string &runsettagname);

    private:

        // Settings
        std::vector<int> v_fvalues;
        TStyle *fpStyle;
        TCanvas *fpCanvas;
        TLegend *fpLegend;
        std::vector<int> colors = {kRed+2, kBlue+2, kGreen+2, kOrange+2, kMagenta+2, kCyan+2, kYellow+2, kViolet+2};
        std::vector<int> markers = {47, 34, 20, 21, 22, 23, 33, 43, 29};

        int resomin = 0;
        int resomax = 25;
        int meanmin = 200;
        int meanmax = 800;
        int stdmin = 0;
        int stdmax = 100;
};

#endif