#include <iostream>
#include "Process.h"
#include "Draw.h"

int main() {
	std::string filename =              "Process.root";
	std::string type =                  "cosmics";
	std::string run =                   "TreeMaker_dog1_00001022";
	std::string tag =                   "_0000-45";
	std::string comment =               "_RCHighLow";
	std::string fileName =           	run + tag + comment;
	std::string inputFolderPath =    	"~/Documents/Code/CPP/CPP_Projects/Read_TreeMaker/ROOT_files/" + type + "/" + run + "/" + comment;
	std::string inputFilePath =     	inputFolderPath + "/" + fileName + ".root";
	std::string outputROOTFolder =   	"~/Documents/Code/CPP/CPP_Projects/Read_TreeMaker/Output_ROOT/" + type + "/" + run + "/" + comment;
	std::string outputPDFFolder =		"~/Documents/Code/CPP/CPP_Projects/Read_TreeMaker/Output_PDF/"  + type + "/" + run + "/" + comment;

	Process process;
	process.             SetType(type);
	process.             SetRun(run);
	process.             SetTag(tag);
	process.             SetComment(comment);
	process.             SetFileName(fileName);
	process.             SetInputFile(inputFilePath);
	process.             SetOutputROOTfolder(outputROOTFolder);
	process.             Run();

	// Write
	TFile filewrite(filename.c_str(), "RECREATE");
	process.Write("Process");
	filewrite.Close();

	// Read
	TFile fileread(filename.c_str(), "READ");
	Process &objread = dynamic_cast<Process&>(*fileread.Get("Process"));
	std::cout << "ROOT path output: " << objread.GetROOTOutputPath() << std::endl;
	std::cout << "Mean of hist fph1f_WF: " << objread.GetMean() << std::endl;
	fileread.Close();

	std::cout << "Draw class section:" << std::endl;
	Draw draw;
	draw.				SetOutputPDFfolder(outputPDFFolder);
	draw.				SetOutputPDFrealpath(objread.GetROOTOutputPath());
	draw.				CompareRuns({process.GetROOTOutputPath()});
	return 0;
}
