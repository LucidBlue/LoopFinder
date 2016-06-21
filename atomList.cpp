/*
 *  atomList.cpp
 *  PChemAnalysis
 *
 *  Created by Brad Sheneman
 *  10/20/12
 *
 */

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "atomList.h"
#include <cmath>
#include <locale>

#include "cleaning.h"

AtomData::AtomData()
{
	atom_num = 0;
	residue_num = 0;
	chain_id = ' ';
	atom_type = "";
	atom_element = "";
	residue_type = "";

	x_cd = 0;
	y_cd = 0;
	z_cd = 0;

	occupancy = 0;
	temp_factor = 0;

	charge = 0;
	interface_atom = false;
}

AtomData::AtomData(const AtomData &source)
{
	atom_num = source.atom_num;
	residue_num = source.residue_num;
	chain_id = source.chain_id;
	atom_type = source.atom_type;
	atom_element = source.atom_element;
	residue_type = source.residue_type;

	x_cd = source.x_cd;
	y_cd = source.y_cd;
	z_cd = source.z_cd;

	occupancy = source.occupancy;
	temp_factor = source.temp_factor;

	charge = source.charge;
	interface_atom = source.interface_atom;
	interfaces = source.interfaces;
}

AtomData& AtomData::operator= (const AtomData &source)
{
	atom_num = source.atom_num;
	residue_num = source.residue_num;
	chain_id = source.chain_id;
	atom_type = source.atom_type;
	atom_element = source.atom_element;
	residue_type = source.residue_type;

	x_cd = source.x_cd;
	y_cd = source.y_cd;
	z_cd = source.z_cd;

	occupancy = source.occupancy;
	temp_factor = source.temp_factor;

	charge = source.charge;
	interface_atom = source.interface_atom;
	interfaces = source.interfaces;

	return *this;
}

bool ProteinComplex::FindChainPair(std::vector<std::string>& pair_list, std::string chain_pair)
{
	std::string pair_order_1, pair_order_2;
	std::stringstream ss;
	bool found_pair = false;

	char chain1 = chain_pair[0];
	char chain2 = chain_pair[1];

	ss << chain1 << chain2;
	ss >> pair_order_1;
	ss.clear();

	ss << chain2 << chain1;
	ss >> pair_order_2;
	ss.clear();

	for (int i = 0; i < pair_list.size(); i++)
	{
		if (pair_list[i] == pair_order_1 || pair_list[i] == pair_order_2)
			found_pair = true;
	}

	return found_pair;
}

// load a Protein Data Bank file. Each line corresponds to a single atom
bool ProteinComplex::LoadPDB2(bfs::path filename)
{
	bfs::ifstream infile(filename);
	std::string current_line;

	if(!infile.is_open())
	{
		std::cout << "File not found: " << filename.string() <<"\n";
		return false;
	}

	while (!infile.eof())
	{
		std::getline(infile, current_line);
		if (!infile.eof())
		{
			if(current_line.substr(0, 4) == "ATOM")
				break;
		}
	}

	while (!infile.eof())
	{
		if(current_line.substr(0, 4) == "ATOM")
		{
			AtomData current_atom;

			std::istringstream current_val(current_line.substr(6, 5));
			current_val >> current_atom.atom_num;
			current_val.clear();

			current_val.str(current_line.substr(12, 4));
			current_val >> current_atom.atom_type;
			current_val.clear();

			current_val.str(current_line.substr(17, 3));
			current_val >> current_atom.residue_type;
			current_val.clear();

			current_val.str(current_line.substr(21, 1));
			current_val >> current_atom.chain_id;
			current_val.clear();

			current_val.str(current_line.substr(22, 4));
			current_val >> current_atom.residue_num;
			current_val.clear();

			current_val.str(current_line.substr(30, 8));
			current_val >> current_atom.x_cd;
			current_val.clear();

			current_val.str(current_line.substr(38, 8));
			current_val >> current_atom.y_cd;
			current_val.clear();

			current_val.str(current_line.substr(46, 8));
			current_val >> current_atom.z_cd;
			current_val.clear();

			current_val.str(current_line.substr(54, 6));
			current_val >> current_atom.occupancy;
			current_val.clear();

			current_val.str(current_line.substr(60, 6));
			current_val >> current_atom.temp_factor;
			current_val.clear();

			current_val.str(current_line.substr(76, 2));
			current_val >> current_atom.atom_element;
			current_val.clear();

			InsertAtomData(current_atom);
			std::getline(infile, current_line);
		}
		else
			std::getline(infile, current_line);
	}
	infile.close();

	if (ComplexAtomData.size() <= 1)
	{
		return false;
	}
	return true;
}

// Find the corresponding chain for this atom and insert.
// Create a new chain if this is the first atom of the chain
void ProteinComplex::InsertAtomData(AtomData& atom)
{
	char current_chain = ' ';

	for (int i = 0; i < ComplexAtomData.size(); i++)
	{
		current_chain = ComplexAtomData[i][0].chain_id;
		if (current_chain == atom.chain_id)
		{
			ComplexAtomData[i].push_back(atom);
			break;
		}
	}

	if (current_chain != atom.chain_id)
	{
		std::vector<AtomData> Current_Chain;
		Current_Chain.push_back(atom);
		ComplexAtomData.push_back(Current_Chain);
	}
}

// Simple 3D distance between two atoms.
double ProteinComplex::AtomDistCalc(AtomData& atom1, AtomData& atom2)
{
	float x_1, x_2, y_1, y_2, z_1, z_2;
	x_1 = atom1.x_cd;
	x_2 = atom2.x_cd;
	y_1 = atom1.y_cd;
	y_2 = atom2.y_cd;
	z_1 = atom1.z_cd;
	z_2 = atom2.z_cd;

	return (sqrt(pow((x_2-x_1),2)+pow((y_2-y_1),2)+pow((z_2-z_1),2)));
}

// used only for debugging purposes
// void ProteinComplex::TestCalc()
// {
// 	std::cout<< ComplexAtomData[0][0].atom_num << " ";
// 	std::cout<< ComplexAtomData[0][6].atom_num << " ";
// 	std::cout<< ComplexAtomData[0][14].atom_num << " ";
// 	std::cout<< ComplexAtomData[0][18].atom_num << " " << std::endl;
// 	std::cout<< AtomDistCalc(ComplexAtomData[0][0], ComplexAtomData[0][6]) << std::endl;
// 	std::cout<< AtomDistCalc(ComplexAtomData[0][6], ComplexAtomData[0][14]) << std::endl;
// 	std::cout<< AtomDistCalc(ComplexAtomData[0][14], ComplexAtomData[0][18]) << std::endl;
// }

// calls distance calc for every pair of atoms on different chains
// and stores them in a vector of atom pairs
void ProteinComplex::AllAtomsDistCalc(double bind_distance, bool aCarbons)
{
	// iterate over possible first chains from complex
	for (std::vector<std::vector<AtomData> >::iterator ch1 = ComplexAtomData.begin(); ch1 < ComplexAtomData.end(); ch1++)
	{
		// iterate over atoms from first chain
		for (std::vector<AtomData>::iterator at1 = ch1->begin(); at1 < ch1->end(); at1++)
		{
			// iterate over possible second chains from complex
			for (std::vector<std::vector<AtomData> >::iterator ch2 = ch1 + 1; ch2 < ComplexAtomData.end(); ch2++)
			{
				// iterate over atoms from second chain
				for (std::vector<AtomData>::iterator at2 = ch2->begin(); at2 < ch2->end(); at2++)
				{
					if (AtomDistCalc(*at1, *at2) < bind_distance)
					{
						char chain1 = at1->chain_id;
						char chain2 = at2->chain_id;
						std::stringstream ss;
						std::string chain_pair;
						ss << chain1 << chain2;
						ss >> chain_pair;
						ss.clear();

						if (aCarbons)
						{
							if (!FindChainPair(at1->interfaces, chain_pair) && at1->atom_type == "CA")
								at1->interfaces.push_back(chain_pair);

							if (!FindChainPair(at2->interfaces, chain_pair) && at2->atom_type == "CA")
								at2->interfaces.push_back(chain_pair);
						}
						else{
							if (!FindChainPair(at1->interfaces, chain_pair))
								at1->interfaces.push_back(chain_pair);

							if (!FindChainPair(at2->interfaces, chain_pair))
								at2->interfaces.push_back(chain_pair);
						}
						/*
						ComplexAtomData[i][j] = atom1;
						ComplexAtomData[a][b] = atom2;
						for (int y = 0; y < atom1.interfaces.size(); y++)
							std::cout << atom1.atom_num << " " << atom1.interfaces.size() << " " << atom1.interfaces[y] << " ";

						for (int y = 0; y < atom2.interfaces.size(); y++)
							std::cout << atom2.atom_num << " " << atom2.interfaces.size() << " " << atom2.interfaces[y] << " ";
						std::cout << std::endl;
						*/

						// we have found an interface residue so increment sum for the two chains
						// but only if they are not already marked (and thus have already been counted toward sum)
						if (!at1->interface_atom)
							NumInterfaceRes[chain1]++;
						if (!at2->interface_atom)
							NumInterfaceRes[chain2]++;

						// and set interface_atom to true, regardless
						at1->interface_atom = true;
						at2->interface_atom = true;
						/*
						std::cout << at1->atom_num << " " << at1->atom_type << " " << at1->interface_atom << " ";
						if (at1->interfaces.size() > 0)
							std::cout << at1->interfaces[0];
						std::cout << std::endl;
						std::cout << at2->atom_num << " " << at2->atom_type << " " << at2->interface_atom << " ";
						if (at2->interfaces.size() > 0)
							std::cout << at2->interfaces[0];
						std::cout << std::endl;
						*/

					}
				}
			}
		}
	}
}

// extracts all the residues and marks them as interface or not,
// then inserts them with data corresponding to the residue a-carbon
void ProteinComplex::ExtractResidues()
{
	for (std::vector<std::vector<AtomData> >::iterator chainI = ComplexAtomData.begin(); chainI < ComplexAtomData.end(); chainI++)
	{
		std::vector<ResidueData> Current_Chain;
		AtomData test_res = *chainI->begin();

		bool interface = false;
		AtomData current_ACarbon;
		std::vector<std::string> res_chain_pairs;

		for (std::vector<AtomData>::iterator atomI = chainI->begin(); atomI < chainI->end(); atomI++)
		{
			AtomData current_res = *atomI;
			//populate proxy vector res_chain_pairs with any pairs not already found
			for (int a = 0; a < current_res.interfaces.size(); a++)
			{
				std::string current_pair = current_res.interfaces[a];
				if (!FindChainPair(res_chain_pairs, current_pair))
					res_chain_pairs.push_back(current_pair);
			}
			if (current_res.residue_num == test_res.residue_num)
			{
				if (current_res.atom_type == "CA")
					current_ACarbon = *atomI;
				if (current_res.interface_atom == true)
					interface = true;
				/*
				if (current_res.interfaces.size() > 0){
				for (int z = 0; z < current_res.interfaces.size(); z++)
					std::cout << current_res.interfaces[z] << " ";
				std::cout << std::endl;
				}*/
			}
			else
			{
				// otherwise, we've moved on to the next residue
				// and we have to insert the data for the previous residue
				ResidueData ResData;
				ResData.aCarbon = current_ACarbon;
				ResData.interface_res = interface;
				ResData.interfaces = res_chain_pairs;


				Current_Chain.push_back(ResData);

				/*
				if (ResData.interfaces.size() > 0){
				std::cout << "Res info: " << ResData.interfaces.size() << " ";
				for (int z = 0; z < ResData.interfaces.size(); z++)
					std::cout << ResData.interfaces[z] << " ";
				std::cout << std::endl;
				}*/

				test_res = current_res;
				interface = false;
				AtomData blank_atom;
				current_ACarbon = blank_atom;
				res_chain_pairs.clear();
				atomI--;
			}
		}
		ComplexResidues.push_back(Current_Chain);
	}
}

void ProteinComplex::PrintResidues(bfs::path filename)
{
	bfs::ofstream outfile(filename);

	for (int i = 0; i < ComplexResidues.size(); i++)
	{
		for (int j = 0; j < ComplexResidues[i].size(); j++)
		{
			outfile << std::left << std::setw(10) << ComplexResidues[i][j].aCarbon.chain_id;
			outfile << std::left << std::setw(10) << ComplexResidues[i][j].aCarbon.residue_num;
			outfile << std::left << std::setw(10) << ComplexResidues[i][j].aCarbon.residue_type;
			outfile << std::left << std::setw(10) << ComplexResidues[i][j].interface_res << std::endl;
		}
	}
}

void ProteinComplex::LoopFinder(int chain_index, ParamData params)
{
	double max_dist = params.max_dist;
	int min_res = params.min_res;
	int max_res = params.max_res;
	double min_perc = params.min_perc;
	double len_factor = params.len_factor;

	for (resVectorIter i = ComplexResidues[chain_index].begin(); i < ComplexResidues[chain_index].end(); i++)
	{
		std::vector<ResidueData> TempLoop;
		int num_res_interface = 0;
		int num_res_total = 0;
		int current_loop_length = 0;
		int first_res_num = i->aCarbon.residue_num;

		for (resVectorIter j = i; (current_loop_length < max_res) && (j < ComplexResidues[chain_index].end()); j++)
		{
			double distance;
			double percent_interface;

			ResidueData current_res = *j;
			int current_res_num = current_res.aCarbon.residue_num;
			current_loop_length = current_res_num - first_res_num + 1;

			TempLoop.push_back(current_res);
			num_res_total++;

			// average length of an amino acid is about 3.2 angstroms (this can be modified if necessary)
			// linker length shouldn't be more than about half of the loop length
			double max_linker_length = len_factor*(TempLoop.size());

			if (j->interface_res)
				num_res_interface++;

			distance = AtomDistCalc(TempLoop.front().aCarbon, TempLoop.back().aCarbon);
			percent_interface = (double(num_res_interface)/double(num_res_total))*100;

			// if loop is large enough and distance and % interface criteria are met, add to 'loops' vector
			if ((current_loop_length >= min_res)
				&& (distance <= max_dist)
				&& (distance <= max_linker_length)
				&& (percent_interface >= min_perc))
			{
				//std::cout << "nrt: " << num_res_total << " nri: " << num_res_interface << " pi: " << percent_interface << std::endl;
				TempLoop.back().distance_to_start = distance;

				LoopData This_Loop;

				std::vector<std::string> pair_duplicates;
				std::vector<std::string> loop_interactions;
				std::vector<std::pair<std::string, double> > loop_interactions_with_perc;

				// create list of duplicates, and list excluding duplicates
				for (resVectorIter loopIter = TempLoop.begin(); loopIter < TempLoop.end(); loopIter++)
				{
					for (chainPairIter pairIter = loopIter->interfaces.begin(); pairIter < loopIter->interfaces.end(); pairIter++)
					{
						pair_duplicates.push_back(*pairIter);

						if (!FindChainPair(loop_interactions, *pairIter))
							loop_interactions.push_back(*pairIter);
					}
				}

				// only update This_Loop.Interactions if the number of dupes is at least loop_length*loop_interface_perc
				for (chainPairIter pairIter = loop_interactions.begin(); pairIter < loop_interactions.end(); pairIter++)
				{
					int num_res_on_this_interface = 0;
					for (chainPairIter dupeIter = pair_duplicates.begin(); dupeIter < pair_duplicates.end(); dupeIter++)
					{
						if (*dupeIter == *pairIter)
							num_res_on_this_interface++;
					}

					double percent_res_on_interface = 100*((double)num_res_on_this_interface)/((double)TempLoop.size());
					std::pair<std::string, double> this_interface(*pairIter, percent_res_on_interface);

					loop_interactions_with_perc.push_back(this_interface);

					if (!FindChainPair(ComplexInteractions, *pairIter) && percent_res_on_interface >= params.loop_interface_perc)
						ComplexInteractions.push_back(*pairIter);
				}

				//std::cout << percent_interface << "% " << distance << std::endl;
				This_Loop.LoopResidues = TempLoop;
				This_Loop.Interactions = loop_interactions_with_perc;
				ComplexLoops.push_back(This_Loop);
			}

		}
	}
}

void ProteinComplex::ExtractLoops(ParamData params)
{
	for (int i = 0; i < ComplexResidues.size(); i++)
		LoopFinder(i, params);
}

void ProteinComplex::PrintOutput(bfs::path input_file, bfs::path output, ParamData params)
{
	double bdist = params.bind_dist;

	bfs::path output_file_loop = output;
	bfs::path output_file_PDB = output;
	output_file_loop /= "Loops.txt";
	output_file_PDB /= "cmd_line_input.txt";

	bfs::ofstream outfile_loop(output_file_loop, bfs::ofstream::app);
	bfs::ofstream outfile_PDB(output_file_PDB, bfs::ofstream::app);

	std::string output_ddg = ChangeFilename(input_file, "_ddg", "");
	std::string pdb_code = ChangeFilename(input_file, "", "");

	for (std::vector<std::string>::iterator pairIter = ComplexInteractions.begin(); pairIter < ComplexInteractions.end(); pairIter++)
	{
		char first, second;
		std::stringstream ss;

		ss << (*pairIter)[0];
		ss >> first;
		ss.clear();

		ss << (*pairIter)[1];
		ss >> second;

		std::stringstream ssfilext;

		ssfilext << "_" << first << second << "_c";
		std::string filext;
		ssfilext >> filext;


		std::string filename_clean = ChangeFilename(input_file, filext, ".pdb");

		outfile_PDB << " --pdb_filename=" << filename_clean
		<< " --partners=" << first << "_" << second
		<< " --interface_cutoff=" << bdist << " --trials=" << 20
		<< " --trial_output="<< output_ddg << std::endl;
	}

	// print all lines to loops
	for (std::vector<LoopData>::iterator loopIter = ComplexLoops.begin(); loopIter < ComplexLoops.end(); loopIter++)
	{

		int loop_length = loopIter->LoopResidues.back().aCarbon.residue_num - loopIter->LoopResidues.front().aCarbon.residue_num + 1;

		outfile_loop.precision(5);

		outfile_loop << pdb_code << " "
		<< std::left << std::setw(3) << loopIter->LoopResidues.front().aCarbon.chain_id
		<< std::left << std::setw(4) << loop_length
		<< std::left << std::setw(7) << loopIter->LoopResidues.back().distance_to_start
		<< " " << "num res: "
		<< NumInterfaceRes[loopIter->LoopResidues.front().aCarbon.chain_id] << " ";

		for (resVectorIter resIter = loopIter->LoopResidues.begin(); resIter < loopIter->LoopResidues.end(); resIter++)
		{
			ResidueData Current_Res = *resIter;
			outfile_loop << std::left << std::setw(5) << Current_Res.aCarbon.residue_type
			<< std::left << std::setw(5) << Current_Res.aCarbon.residue_num << "  ";
		}

		outfile_loop << std::endl << "INTERFACES: ";

		for (sdpVectorIter pairIter = loopIter->Interactions.begin(); pairIter < loopIter->Interactions.end(); pairIter++)
		{
			std::pair<std::string, double> Current_Pair = *pairIter;
			outfile_loop << std::left << std::setw(5) << pairIter->first
			<< std::left << std::setw(5) << pairIter->second << "  ";
		}
		outfile_loop << std::endl;
	}
	outfile_loop << std::endl;
	outfile_loop.close();
	outfile_PDB.close();
}
