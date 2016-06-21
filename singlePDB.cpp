#include <iostream>
#include <fstream>
#include "singlePDB.h"
#include "atomList.h"

void SinglePDB (bfs::path input, bfs::path output, ParamData params)
{
	ProteinComplex testComplex;

	std::string pdb_and_partners = input.leaf().string().substr(0,7);
	std::string filename_res = pdb_and_partners + "_res.txt";

	bfs::path file_clean = input;

	bfs::path file_res = output;
	file_res /= filename_res;

	if (testComplex.LoadPDB2(file_clean) == false)
		return;

	testComplex.AllAtomsDistCalc(params.bind_dist, params.aCarbons_only);

	testComplex.ExtractResidues();
	testComplex.PrintResidues(file_res);

	testComplex.ExtractLoops(params);

	testComplex.PrintOutput(input, output, params);
}
