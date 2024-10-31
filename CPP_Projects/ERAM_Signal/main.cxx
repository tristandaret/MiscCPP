#include "PointModel.h"
#include "TrackModel.h"
#include <iomanip>

int main(int argc, char const *argv[])
{
	// DrawPointModel();
	TrackModel *p_trackmodel = new TrackModel();
	for(int i=0; i<1e4; i++) p_trackmodel->GetRealCharge(800, 11.28, i*7/1e4, 0, 120, 300);

	std::cout << std::fixed << std::setprecision(5);

	std::cout << "Time: " << p_trackmodel->time << std::endl;
	std::cout << "time per event: " << p_trackmodel->time/1e4 << std::endl;
}