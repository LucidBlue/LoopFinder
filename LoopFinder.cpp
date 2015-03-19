#include <iostream>
#include <fstream>
#include <vector>
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"

#include "singlePDB.h"
#include "cleaning.h"

namespace bfs = boost::filesystem;
using namespace std;

int main(int argc, char* const argv[])
{
	ParamData params;

	params.bind_dist = 5;
	params.min_res = 4;
	params.max_res = 8;
	params.min_perc = 80;
	params.max_dist = 9.6;
	params.len_factor = 1.75;
	params.aCarbons_only = false;
	params.loop_interface_perc = 0.5;
	params.num_partners = 2;

	bfs::path input;

	for (int i = 0; i < argc; i++)
	{
		if (string(argv[i]) == "-input")
			input = (argv[++i]);
		else if (string(argv[i]) == "-bdist")
			params.bind_dist = atof(argv[++i]);
		else if (string(argv[i]) == "-minres")
			params.min_res = atoi(argv[++i]);
		else if (string(argv[i]) == "-maxres")
			params.max_res = atoi(argv[++i]);
		else if (string(argv[i]) == "-minperc")
			params.min_perc = atof(argv[++i]);
		else if (string(argv[i]) == "-maxdist")
			params.max_dist = atof(argv[++i]);
		else if (string(argv[i]) == "-lfactor")
			params.len_factor = atof(argv[++i]);
		else if (string(argv[i]) == "-liperc")
			params.loop_interface_perc = atof(argv[++i]);
		else if (string(argv[i]) == "-nchain")
			params.num_partners = atoi(argv[++i]);
		else if (string(argv[i]) == "-acarbs")
		{
			params.aCarbons_only = true;
			++i;
		}
		
	}

	if (bfs::exists(input))
	{
		bfs::path output (input);
		output /= "results";
		bfs::create_directories(output);

		for (bfs::directory_iterator input_iter (input); 
			input_iter != bfs::directory_iterator(); input_iter++)
		{
			if (bfs::is_regular_file(*input_iter))
			{
				bfs::path output_path_strip(output);
				output_path_strip /= input_iter->path().leaf().string();

				std::vector<vector<char> > ChainDuplicates = FindDuplicates(input_iter->path());
				std::vector<char> GoodChains = RemoveDuplicates(params, ChainDuplicates);

				CleanPDB(input_iter->path(), output_path_strip, GoodChains);
				SplitPDB(output_path_strip, GoodChains);
			}
		}
		

		for (bfs::directory_iterator input_iter (output); 
			input_iter != bfs::directory_iterator(); input_iter++)
		{
			if (bfs::is_regular_file(*input_iter) && (input_iter->path().leaf().string()).find("_c.pdb") == 7)
				SinglePDB(input_iter->path(), output, params);
		}
	}
	else
	{
		std::cout << "Not a valid folder" << std::endl;
		return -1;
	}
	return 0;
}