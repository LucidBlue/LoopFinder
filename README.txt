LoopFinder Readme
Bradley Sheneman
Last updated: 2016/06/21

****************************
Compilation:
****************************
1. Make sure you have the (filesystem + fstream) boost libraries installed

2. run the following cmd line (on Linux):

	g++ -std=c++11 -o LoopFinder loopFinder.cpp atomList.cpp cleaning.cpp singlePDB.cpp -lboost_system -lboost_filesystem

3. If that worked, you are done!

****************************
Setup and Run:
****************************

1. Create a folder of raw PDB files and note its filepath
	In the example directory, this is `/TestFiles`

2. Run LoopFinder from its location, using the parameters shown below
	(a batch file will work best for multiple runs)

3. A "results" folder will appear in your input folder, containing multiple
	intermediate files for each PDB, and a Loops.txt with data for all files

****************************
Advanced Input:
****************************

Example cmd line on Linux:
./LoopFinder -input TestFiles -bdist 6.5 -minres 4 -maxres 8 -minperc 80 -maxdist 6.2 -lfactor 1.75 -nchain 3

The flags:
-input: 	specifies the path to the location of your pdb files (relative or absolute)

-bdist: 	maximum distance (in angstroms) to be considered a bond between
			two atoms of separate residues (used when determining the interface)

-minres: 	minimum number of residues in a loop

-maxres: 	maximum number of residues in a loop

-minperc: 	minimum percentage of interface residues in a loop for it to be considered

-maxdist: 	maximum distance (in angstroms) allowed between ends of a given loop.
			(i.e. linker distance)

-lfactor:	'linker length factor' = ((average aa length) * (max fraction of total loop length)

			e.g. for an average amino acid length of 3.5 angstroms, and a max linker length of
			half the length of the loop: lfactor = (3.5) * (1/2) = 1.75

-nchain		number of binding partners in each file.
			*must be the same for all files in a single run of LoopFinder!*

****************************
Interpreting Results:
****************************

"filename.pdb" : this is identical to the original, but with all non-ATOM lines removed.
	It also has duplicates removed, according to my algorithm


"filename_<partners>_c.pdb" : contains only the lines for the two chains in 'partners'
	the number of these files created will depend on the number of chains specified above in nchains.
	e.g. if a file contains 3 chains (A, B and C) and -nchain is specified as 3, 3 files will be created:
	filename_AB_c.pdb, filename_AC_c.pdb and filename_BC_c.pdb


"filename_<partners>_res.txt" : list of all residues from the two partner chains,
	and whether they are part of the interface.


	the columns are as follows:

	chain ID, residue number, residue type, interface (0/1)


"Loops.txt" : list of all loops, with linker distance and % interface atoms

	this file does not get overwritten with subsequent runs of LoopFinder
	make sure to delete the old results folder for subsequent runs,
	or start with a fresh set of PDBs in a new location
