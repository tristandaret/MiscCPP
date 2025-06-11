#include "Draw.h"
#include "Process.h"
#include "Includes.h"
#include "Util.h"

int main()
{

   int do_run = 1;
   int do_draw = 1;
   int do_draw_comparison = 0;
   // Output structure: type -> run -> comment -> fileName

   std::vector<std::string> v_types;
   std::vector<std::string> v_runs;
   std::vector<std::string> v_comments;

   v_types.push_back("beam");
   v_runs.push_back("treeMaker_beam_may2024");
   v_comments.push_back("nd280_14.32");

   // v_types.push_back("beam");
   // v_runs.push_back("tree_mu-_40-4895MeV_sandmu");
   // v_comments.push_back("nd280_14.32");

   // v_types.push_back("cosmics");
   // v_runs.push_back("TreeMaker_dog1_00001148");
   // v_comments.push_back("nd280_14.32");

   // v_types.push_back("cosmics");
   // v_runs.push_back("tree_mu-_40-4895MeV_cosmics");
   // v_comments.push_back("BFieldON_14.32");

   std::string tag = "";

   std::vector<std::string> v_filepaths;

   Draw draw;
   draw.SetOutputComparisonFolder("Output_PDF/" + v_types[0]);

   for (int i = 0; i < (int)v_types.size(); i++) {
      std::string type = v_types[i];
      std::string run = v_runs[i];
      std::string comment = v_comments[i];

      std::string fileName = run;
      if (tag != "")
         fileName += "_" + tag;
      if (comment != "")
         fileName += "_" + comment;

      // data
      std::string inputFolderPath = "/local/home/td263283/Documents/Code/CPP/"
                                    "CPP_Projects/Read_TreeMaker/ROOT_files/" +
                                    type + "/" + run;
      MakeMyDir(inputFolderPath);
      std::string inputFilePath = inputFolderPath + "/" + fileName + ".root";
      if (!(fopen(inputFilePath.c_str(), "r"))) {
         std::cout << "Not found: " << inputFilePath << std::endl;
         int result = system(("scp cca9.in2p3.fr:~/public/Output_root/" + run + "/" +
                              comment + "/" + fileName + ".root " + inputFolderPath)
                                .c_str());
         if (result == 0)
            std::cout << "Downloaded: " << fileName << ".root" << std::endl;
         else {
            std::cerr << "Error: scp command failed with exit code " << result
                      << std::endl;
            continue;
         }
      }

      // draw
      std::string outputComparisonFolder = "Output_PDF/" + type + "/" + run;
      std::string outputPDFFolder = "Output_PDF/" + type + "/" + run + "/" + comment;
      std::string outputROOTFolder = "Output_ROOT/" + type + "/" + run + "/" + comment;
      MakeMyDir(outputPDFFolder);
      MakeMyDir(outputROOTFolder);

      Process process;
      process.SetType(type);
      process.SetRun(run);
      process.SetTag(tag);
      process.SetComment(comment);
      process.SetFileName(fileName);
      process.SetInputFile(inputFilePath);
      process.SetOutputROOTfolder(outputROOTFolder);
      process.SetCuts();
      if (do_run) {
         process.Run();
         TFile tfile(process.GetROOTOutputPath().c_str(), "RECREATE");
         process.Write();
         tfile.Close();
      }

      v_filepaths.push_back(process.GetROOTOutputPath());
      draw.SetOutputPDFFolder(outputPDFFolder);
      if (do_draw)
         draw.Run(v_filepaths.back());
   }
   if (v_comments.size() > 1 && do_draw_comparison)
      draw.Compare(v_filepaths, "runs");

   return 0;
}