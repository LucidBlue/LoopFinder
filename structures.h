#ifndef STRUCTURES_H
#define STRUCTURES_H

struct ParamData
{
	double bind_dist;

	int min_res;
	int max_res;
	double min_perc;
	double max_dist;
	double len_factor;
	double loop_interface_perc;
	int num_partners;
	bool aCarbons_only;
};

#endif