#include "cleaning.h"

#include <iostream>
#include <sstream>
#include <iomanip>

std::string ChangeFilename(bfs::path input_file, std::string append, std::string extension)
{

	std::string file_name = input_file.filename().string();
	std::string file_name_new;

	file_name.erase(file_name.begin() + 4, file_name.end());
	
	file_name_new = file_name + append + extension;

	return file_name_new;
}

std::vector<std::vector<char> > FindDuplicates(bfs::path filename)
{
	bfs::ifstream infile(filename);
	
	std::vector<std::vector<char> > ChainDuplicates;

	if(!infile.is_open())
	{
		std::cout << "File not found\n";
		std::vector<std::vector<char> > empty_vec;
		return empty_vec;
	}
	
	for (int chain_iter = 0; !infile.eof(); chain_iter++)
	{
		std::string current_line;
 		std::getline(infile, current_line);

		if (current_line.find("SOURCE", 0) == 0)
			break;
		
		if (current_line.find("COMPND", 0) == 0)
		{
			if (current_line.find("CHAIN", 0) < 12)
			{
				std::vector<char> TempDuplicates;
				
				for (int iter = current_line.find("CHAIN", 0) + 6; iter < current_line.size(); iter++)
				{
					if (std::isalpha(current_line[iter], std::locale()))
						TempDuplicates.push_back(current_line[iter]);
				}
				
				ChainDuplicates.push_back(TempDuplicates);
			}
		}
	}
	infile.close();

	return ChainDuplicates;
}

std::vector<char> RemoveDuplicates(ParamData params, std::vector<std::vector<char> > ChainDuplicates)
{
	std::vector<char> Chains;

	std::vector<std::vector<char> >::iterator vecIter = ChainDuplicates.begin();
	int nchains = 0;
	while(nchains < params.num_partners && !ChainDuplicates.empty())
	{
		if(vecIter == ChainDuplicates.end())
			vecIter = ChainDuplicates.begin();
		if(vecIter->empty())
			vecIter = ChainDuplicates.erase(vecIter);
		else
		{
			Chains.push_back(*vecIter->begin());
			vecIter->erase(vecIter->begin());
			nchains++;
			vecIter++;
		}
	}

	return Chains;
}

bool IsDuplicate(char chain_ID_test, std::vector<char> Chains)
{
	bool not_in_chains = true;

	for (int i = 0; i < Chains.size(); i++)
	{
		if (chain_ID_test == Chains[i])
			not_in_chains = false;
	}
	return not_in_chains;
}

void CleanPDB(bfs::path input, bfs::path output, std::vector<char> Chains)
{
	bfs::ifstream infile(input);
	bfs::ofstream outfile(output);
	
	std::vector<char> Visited_Chains;
	

	if(!infile.is_open())
	{
		std::cout << "File not found\n";
		return;
	}

	char current_chain = '*';

	for (int i = 0; !infile.eof(); i++)
	{
		std::string current_line;
		std::getline(infile, current_line);
		

		if (current_line.find("ATOM", 0) == 0)
		{
			bool visited = false;

			for (int a = 0; a < Visited_Chains.size(); a++)
			{
				if (current_line[21] == Visited_Chains[a])
					visited = true;
			}

			if (current_chain != current_line[21] && current_chain != '*')
				Visited_Chains.push_back(current_chain);

			if (visited == false && !IsDuplicate(current_line[21], Chains))
				outfile << current_line << std::endl;

			current_chain = current_line[21];
		}
	}

	infile.close();
	outfile.close();
}

void SplitPDB(bfs::path input, std::vector<char> chains)
{
	for (int i = 0; i < chains.size(); i++)
	{
		for (int j = i+1; j < chains.size(); j++)
		{
			char chain1 = chains[i];
			char chain2 = chains[j];
			bfs::path output(input);

			std::stringstream ss;
			ss << "_" << chain1 << chain2 << "_c";
			std::string to_append;
			ss >> to_append;
			std::string output_file = ChangeFilename(output, to_append, ".pdb");
			output.remove_leaf() /= output_file;

			bfs::ifstream infile(input);
			bfs::ofstream outfile(output);

			if(!infile.is_open())
			{
				std::cout << "File not found\n";
				return;
			}

			while (!infile.eof())
			{
				std::string current_line;
				std::getline(infile, current_line);

				if (current_line.find("ATOM", 0) == 0)
				{
					if (current_line[21] == chain1 || current_line[21] == chain2)
						outfile << current_line << std::endl;
				}
			
			}
		}
	}
}