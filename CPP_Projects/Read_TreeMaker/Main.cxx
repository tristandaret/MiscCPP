#include "Draw.h"
#include "Process.h"
#include "Includes.h"
#include "Util.h"

int main(){

	int do_run = 1;
	int do_draw = 1;
	int do_draw_compcuts = 0;
   // Output structure: type -> run -> comment -> fileName
	std::string type =                  "cosmics";
	// std::string run =                   "TreeMaker_dog1_00001022";
	std::string run =                   "hattree_1148";
	// std::string tag =                   "_0000-45";
	std::string tag =                   "_full";
	// std::vector<std::string>            v_comments{"_RCHighLow", "_master_T04_GFixed"};
	// std::vector<std::string>            v_comments{"_RCHighLow"};
	std::vector<std::string>            v_comments{""};
   	std::vector<std::string>            v_filepaths;


	Draw draw;
	draw.                               SetOutputComparisonFolder("Output_PDF/" + type + "/" + run);

	for(std::string comment : v_comments)
	{
		std::string fileName =           run + tag + comment;

		// data
		std::string inputFolderPath =    "/local/home/td263283/Documents/Code/CPP/CPP_Projects/Read_TreeMaker/ROOT_files/" + type + "/" + run + "/" + comment;
		MakeMyDir(inputFolderPath);
		std::string inputFilePath =      inputFolderPath + "/" + fileName + ".root";
		if(!(fopen(inputFilePath.c_str(), "r")))
		{
			std::cout << "Not found: " << inputFilePath << std::endl;
			int result =                  system(("scp cca9.in2p3.fr:~/public/Output_root/" + run + "/" + comment + "/" + fileName + ".root " + inputFolderPath).c_str());
			if(result == 0)               std::cout << "Downloaded: " << fileName << ".root" << std::endl;
			else{                         std::cerr << "Error: scp command failed with exit code " << result << std::endl; continue;}
		}

		// draw
		std::string outputComparisonFolder = "Output_PDF/"  + type + "/" + run;
		std::string outputPDFFolder =    "Output_PDF/"  + type + "/" + run + "/" + comment;
		std::string outputROOTFolder =   "Output_ROOT/" + type + "/" + run + "/" + comment;
		MakeMyDir(outputPDFFolder);
		MakeMyDir(outputROOTFolder);

		Process process;
		process.             SetType(type);
		process.             SetRun(run);
		process.             SetTag(tag);
		process.             SetComment(comment);
		process.             SetFileName(fileName);
		process.             SetInputFile(inputFilePath);
		process.             SetOutputROOTfolder(outputROOTFolder);
		process.             SetCuts();
		if(do_run){
			process.		Run();
			TFile tfile(process.GetROOTOutputPath().c_str(), "RECREATE");
			process.		Write();
			tfile.			Close();
		}

		v_filepaths.		push_back(process.GetROOTOutputPath());
		draw.				SetOutputPDFFolder(outputPDFFolder);
		if(do_draw) draw.	Run(v_filepaths.back());
	}
	if(v_comments.size() > 1) draw.	CompareRuns(v_filepaths, "comments");

	if(!do_draw_compcuts) return 0;
	v_filepaths.clear();
	v_filepaths.push_back("Output_ROOT/cosmics/TreeMaker_dog1_00001022/_RCHighLow/TreeMaker_dog1_00001022_0000-45_RCHighLow_dir1>0flip.root");
	v_filepaths.push_back("Output_ROOT/cosmics/TreeMaker_dog1_00001022/_RCHighLow/TreeMaker_dog1_00001022_0000-45_RCHighLow_chi2ndf5_25dx_dir1>0flip.root");
	v_filepaths.push_back("Output_ROOT/cosmics/TreeMaker_dog1_00001022/_RCHighLow/TreeMaker_dog1_00001022_0000-45_RCHighLow_chi2ndf5_25dx_200mom1000_dir1>0flip.root");
	v_filepaths.push_back("Output_ROOT/cosmics/TreeMaker_dog1_00001022/_RCHighLow/TreeMaker_dog1_00001022_0000-45_RCHighLow_chi2ndf5_25dx_250mom500_dir1>0flip.root");
	draw.						CompareRuns(v_filepaths, "cuts");

	return 0;
}