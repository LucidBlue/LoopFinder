#ifndef CLEANING_H
#define CLEANING_H

#include <string>
#include <vector>

#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"
namespace bfs = boost::filesystem;

#include "structures.h"

std::string ChangeFilename(bfs::path input_file, std::string append, std::string extension);
std::vector<std::vector<char> > FindDuplicates(bfs::path filename);
std::vector<char> RemoveDuplicates(ParamData params, std::vector<std::vector<char> > ChainDuplicates);
bool IsDuplicate(char chain_ID_test, std::vector<char> Chains);
void CleanPDB(bfs::path input, bfs::path output, std::vector<char> Chains);
void SplitPDB(bfs::path input, std::vector<char> chains);

#endif