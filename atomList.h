// 
// atomList.h
// Created by Brad Sheneman
// 10/20/2012
//

#include <string>
#include <vector>
#include <unordered_map>

#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"
namespace bfs = boost::filesystem;

#include "structures.h"

class AtomData
{
public:
	int atom_num;
	int residue_num;
	char chain_id;
	std::string atom_type;
	std::string atom_element;
	std::string residue_type;

	float x_cd;
	float y_cd;
	float z_cd;
	
	float occupancy;
	float temp_factor;
	int charge;

	bool interface_atom;
	std::vector<std::string> interfaces;
	
	AtomData();
	AtomData(const AtomData &source);
	AtomData& operator= (const AtomData &source);
};


struct AtomPairData
{
	AtomData atom1;
	AtomData atom2;
	
	double distance;
};

struct ResidueData
{
	AtomData aCarbon;
	bool interface_res;
	double distance_to_start;
	std::vector<std::string> interfaces;
};

struct LoopData
{
	std::vector<ResidueData> LoopResidues;
	std::vector<std::pair<std::string, double> > Interactions;
};

typedef std::vector<std::string>::iterator chainPairIter;

typedef std::pair<std::string, double> stringDoublePair;
typedef std::vector<stringDoublePair>::iterator sdpVectorIter;

typedef std::vector<std::vector<AtomData> >::iterator atomComplexIter;
typedef std::vector<AtomData>::iterator atomVectorIter;

typedef std::vector<std::vector<ResidueData> >::iterator resComplexIter;
typedef std::vector<ResidueData>::iterator resVectorIter;

class ProteinComplex
{
public:
	bool FindChainPair(std::vector<std::string>& pair_list, std::string chain_pair);
//  void LoadPDB(bfs::path filename);
	bool LoadPDB2(bfs::path filename);
	void InsertAtomData(AtomData& atom);
	double AtomDistCalc(AtomData& atom1, AtomData& atom2);
	void TestCalc();
	void AllAtomsDistCalc(double bind_distance, bool aCarbons);
	void ExtractResidues();
	void PrintResidues(bfs::path filename);
	void LoopFinder(int chain_index, ParamData params);
	void ExtractLoops(ParamData params);
	void PrintOutput(bfs::path input_file, bfs::path output, ParamData params);

private:
	int __num_models__;
	std::vector<std::string> ComplexInteractions;
	std::vector<std::vector<AtomData> > ComplexAtomData;
	std::vector<std::vector<ResidueData> > ComplexResidues;
	std::vector<LoopData> ComplexLoops;
	std::unordered_map<char,int> NumInterfaceRes;
};