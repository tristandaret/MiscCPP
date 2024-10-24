#include "DrawOuts.h"
#include "Includes.h"
#include "Util.h"

int main(){

    int single =						1;
    int compare_runs =                  0;
    int scan =							0;
    int compare_scans =					0;

	Draw *p_draw =						new Draw();
	p_draw->scantype =					"cosmics";
	p_draw->type =						"TreeMaker_dog1_00000920";
	p_draw->v_comment =					{"_noGRC", "_noRC_GERAM", "_noRC"};
    std::vector<std::string>			v_tags;
    std::vector<std::string>			v_filepath;
    std::vector<std::string>			v_filepathraw;
    
    if(p_draw->scantype == "phi"){
        v_tags.                         insert(v_tags.end(), {
                                        "_mu-_600MeV_x50_y-91_z-180_phi0_theta0_N5000_s0_n5000",
                                        "_mu-_600MeV_x50_y-90_z-180_phi-10_theta0_N5000_s0_n5000",
                                        "_mu-_600MeV_x50_y-85_z-180_phi-20_theta0_N5000_s0_n5000",
                                        "_mu-_600MeV_x50_y-74_z-180_phi-30_theta0_N5000_s0_n5000",
                                        "_mu-_600MeV_x50_y-64_z-180_phi-40_theta0_N5000_s0_n5000",
                                        "_mu-_600MeV_x50_y-60_z-180_phi-45_theta0_N5000_s0_n5000",
                                        "_mu-_600MeV_x50_y-55_z-180_phi-50_theta0_N5000_s0_n5000",
                                        "_mu-_600MeV_x50_y-50_z-180_phi-60_theta0_N5000_s0_n5000",
                                        "_mu-_600MeV_x50_y-43_z-180_phi-70_theta0_N5000_s0_n5000",
                                        "_mu-_600MeV_x50_y-40_z-180_phi-80_theta0_N5000_s0_n5000",
                                        "_mu-_600MeV_x50_y-39_z-180_phi-90_theta0_N5000_s0_n5000"});
    }

    else if(p_draw->scantype == "drift"){
        v_tags.                         insert(v_tags.end(), {
                                        "_mu-_600MeV_x97_y-40_z-200_phi-90_theta0_N5000_s0_n5000",
                                        "_mu-_600MeV_x88_y-40_z-200_phi-90_theta0_N5000_s0_n5000",
                                        "_mu-_600MeV_x78_y-40_z-200_phi-90_theta0_N5000_s0_n5000",
                                        "_mu-_600MeV_x68_y-40_z-200_phi-90_theta0_N5000_s0_n5000",
                                        "_mu-_600MeV_x58_y-40_z-200_phi-90_theta0_N5000_s0_n5000",
                                        "_mu-_600MeV_x48_y-40_z-200_phi-90_theta0_N5000_s0_n5000",
                                        "_mu-_600MeV_x38_y-40_z-200_phi-90_theta0_N5000_s0_n5000",
                                        "_mu-_600MeV_x28_y-40_z-200_phi-90_theta0_N5000_s0_n5000",
                                        "_mu-_600MeV_x18_y-40_z-200_phi-90_theta0_N5000_s0_n5000",
                                        "_mu-_600MeV_x8_y-40_z-200_phi-90_theta0_N5000_s0_n5000",
                                        "_mu-_600MeV_x1_y-40_z-200_phi-90_theta0_N5000_s0_n5000"});
    }

    else{
        v_tags.                         push_back("_0000-8_s0_n23000");
    }


	for(std::string comment : p_draw->v_comment)
	{
		p_draw->comment = comment;
		std::string folderpath = "ROOT_files/" + p_draw->type + "/" + p_draw->scantype + "/" + p_draw->comment;
		MakeMyDir(folderpath);

		for(std::string tag : v_tags)
		{   
			std::string filename = p_draw->type + tag + comment;
			std::string filepath = folderpath + "/" + filename + ".root";
			if(!(fopen(filepath.c_str(), "r")))
			{
				std::cout << "Not found: " << filepath << std::endl;
				int result = system(("scp cca9.in2p3.fr:~/public/Output_root/" + p_draw->scantype + "/" + comment + "/" + filename + ".root " + folderpath).c_str());
				if(result == 0)	std::cout << "Downloaded: " << filename << ".root" << std::endl;
				else{
					std::cerr << "Error: scp command failed with exit code " << result << std::endl;
					continue;
				}
			}
			if(single){
                p_draw->filename = filename;
                p_draw->	Run(filepath);
            }
			v_filepath.			push_back(filepath);
		}
	}

	if(scan) p_draw->			Scan(v_filepath);

	if(compare_scans){
		p_draw->				CompareScans(v_tags, "WF");
		p_draw->				CompareScans(v_tags, "XP");
	}

	delete p_draw;
}