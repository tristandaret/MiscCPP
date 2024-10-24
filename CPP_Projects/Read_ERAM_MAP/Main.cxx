#include "maps_drawing.h"

int main()
{
	Drawing *pDrawing = new Drawing();

	// pDrawing->maps_pdf();
	// pDrawing->maps_csv();
	pDrawing->Distributions();

	delete pDrawing;
	return 0;
}