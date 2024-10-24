#include "Draw.h"
#include "Process.h"
#include "Includes.h"
#include "Util.h"

int main(){

	std::string type =					"cosmics";
	std::string runset =				"TreeMaker_dog1_00001022";
	std::string tag = 					"_0000-45_s0_n28000";
	// std::vector<std::string>			v_comments{"_noGRC", "_noRC_GERAM", "_noRC"};
	std::vector<std::string>			v_comments{"_Ghardware_noRC"};
    std::vector<std::string>			v_tags;
    std::vector<std::string>			v_filepath;
    std::vector<std::string>			v_filepathraw;

	std::string compfolderpath = 		"Output_PDF/" + type + "/" + runset;
	std::string runsettagname = 		runset + tag;

	std::vector<Process*> v_processes;
	Draw *p_draw = 						new Draw();

	for(std::string comment : v_comments)
	{
		std::string filename = 			runset + tag + comment;
		// data
		std::string datafolderpath = 	"ROOT_files/" + type + "/" + runset + "/" + comment;
		MakeMyDir(datafolderpath);
		std::string datafilepath = 		datafolderpath + "/" + filename + ".root";
		if(!(fopen(datafilepath.c_str(), "r")))
		{
			std::cout << "Not found: " << datafilepath << std::endl;
			int result = 				system(("scp cca9.in2p3.fr:~/public/Output_root/" + runset + "/" + comment + "/" + filename + ".root " + datafolderpath).c_str());
			if(result == 0)				std::cout << "Downloaded: " << filename << ".root" << std::endl;
			else{ 						std::cerr << "Error: scp command failed with exit code " << result << std::endl; continue;}
		}
		// draw
		std::string drawfolderpath = 	compfolderpath + "/" + comment;
		MakeMyDir(drawfolderpath);
		v_filepath.						push_back(datafilepath);
		v_processes.					push_back(new Process());
		v_processes.back()->			Datafile(comment, datafilepath, drawfolderpath, filename);
		p_draw->						Run(*v_processes.back());
	}
	if(v_comments.size() > 1) p_draw->	CompareRuns(v_processes, compfolderpath, runsettagname);
}









// int main(){

//     int run =						1;
//     int compare_runs =                  0;
//     int scan =							0;
//     int compare_scans =					0;

// 	Draw *p_process =						new Draw();
// 	p_process->runset =					"cosmics";
// 	p_process->type =						"TreeMaker_dog1_00000920";
// 	p_process->v_comments =					{"_noGRC", "_noRC_GERAM", "_noRC"};
//     std::vector<std::string>			v_tags;
//     std::vector<std::string>			v_filepath;
//     std::vector<std::string>			v_filepathraw;
    
//     if(p_process->runset == "phi"){
//         v_tags.                         insert(v_tags.end(), {
//                                         "_mu-_600MeV_x50_y-91_z-180_phi0_theta0_N5000_s0_n5000",
//                                         "_mu-_600MeV_x50_y-90_z-180_phi-10_theta0_N5000_s0_n5000",
//                                         "_mu-_600MeV_x50_y-85_z-180_phi-20_theta0_N5000_s0_n5000",
//                                         "_mu-_600MeV_x50_y-74_z-180_phi-30_theta0_N5000_s0_n5000",
//                                         "_mu-_600MeV_x50_y-64_z-180_phi-40_theta0_N5000_s0_n5000",
//                                         "_mu-_600MeV_x50_y-60_z-180_phi-45_theta0_N5000_s0_n5000",
//                                         "_mu-_600MeV_x50_y-55_z-180_phi-50_theta0_N5000_s0_n5000",
//                                         "_mu-_600MeV_x50_y-50_z-180_phi-60_theta0_N5000_s0_n5000",
//                                         "_mu-_600MeV_x50_y-43_z-180_phi-70_theta0_N5000_s0_n5000",
//                                         "_mu-_600MeV_x50_y-40_z-180_phi-80_theta0_N5000_s0_n5000",
//                                         "_mu-_600MeV_x50_y-39_z-180_phi-90_theta0_N5000_s0_n5000"});
//     }

//     else if(p_process->runset == "drift"){
//         v_tags.                         insert(v_tags.end(), {
//                                         "_mu-_600MeV_x97_y-40_z-200_phi-90_theta0_N5000_s0_n5000",
//                                         "_mu-_600MeV_x88_y-40_z-200_phi-90_theta0_N5000_s0_n5000",
//                                         "_mu-_600MeV_x78_y-40_z-200_phi-90_theta0_N5000_s0_n5000",
//                                         "_mu-_600MeV_x68_y-40_z-200_phi-90_theta0_N5000_s0_n5000",
//                                         "_mu-_600MeV_x58_y-40_z-200_phi-90_theta0_N5000_s0_n5000",
//                                         "_mu-_600MeV_x48_y-40_z-200_phi-90_theta0_N5000_s0_n5000",
//                                         "_mu-_600MeV_x38_y-40_z-200_phi-90_theta0_N5000_s0_n5000",
//                                         "_mu-_600MeV_x28_y-40_z-200_phi-90_theta0_N5000_s0_n5000",
//                                         "_mu-_600MeV_x18_y-40_z-200_phi-90_theta0_N5000_s0_n5000",
//                                         "_mu-_600MeV_x8_y-40_z-200_phi-90_theta0_N5000_s0_n5000",
//                                         "_mu-_600MeV_x1_y-40_z-200_phi-90_theta0_N5000_s0_n5000"});
//     }

//     else{
//         v_tags.                         push_back("_0000-8_s0_n23000");
//     }


// 	for(std::string comment : p_process->v_comments)
// 	{
// 		p_process->comment = comment;
// 		std::string datafolderpath = "ROOT_files/" + p_process->type + "/" + p_process->runset + "/" + p_process->comment;
// 		MakeMyDir(datafolderpath);

// 		for(std::string tag : v_tags)
// 		{   
// 			std::string filename = p_process->type + tag + comment;
// 			std::string datafilepath = datafolderpath + "/" + filename + ".root";
// 			if(!(fopen(datafilepath.c_str(), "r")))
// 			{
// 				std::cout << "Not found: " << datafilepath << std::endl;
// 				int result = system(("scp cca9.in2p3.fr:~/public/Output_root/" + p_process->runset + "/" + comment + "/" + filename + ".root " + datafolderpath).c_str());
// 				if(result == 0)	std::cout << "Downloaded: " << filename << ".root" << std::endl;
// 				else{
// 					std::cerr << "Error: scp command failed with exit code " << result << std::endl;
// 					continue;
// 				}
// 			}
// 			if(run){
//                 p_process->filename = filename;
//                 p_process->	Run(datafilepath);
//             }
// 			v_filepath.			push_back(datafilepath);
// 		}
// 	}

// 	if(scan) p_process->			Scan(v_filepath);

// 	if(compare_scans){
// 		p_process->				CompareScans(v_tags, "WF");
// 		p_process->				CompareScans(v_tags, "XP");
// 	}

// 	delete p_process;
// }