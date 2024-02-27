### CITATION ###

To cite RDP, please cite McCoy, K.M., Ackerman, M.E. and Grigoryan, G., 2024. A significance score for protein–protein interaction models through random docking. Protein Science, 33(2), p.e4853.



### INTRODUCTION ###

RD p-value is a statistical significance score for protein-protein interaction models. It is based on the probability of getting a model as good or better than the one provided, by randomly docking the true structure's partners together. The user must chose which relative distance metric to use to score each of the random docks and the model, which will then be used to place the model's score into the distribution of scores of the randomly docked models, generating a p-value. Currently complex RMSD, Receptor + Ligand RMSD, and DOCKQ are possible relative distance metrics. It may also be biased so that random docks must include certain binding residues and/or dock at a given angle or below.

Both RMSDs and DOCKQ require atoms to be matched between the model and true structure. This script's atom-matching process automatically figures out which atoms to pair to one another, only requiring the user to provide which chains represent each binding partner in the true structure (i.e. the chains representing each binding partner in the models do not need to be provided and will be automatically figured out).  This functionality is available even in cases of binding partners involving multiple chains, or there being multiple unique binding locations. It can even match multiple chains from the true structure to the same chain on a model, for cases in which multiple chains in the model have the same identifier. This is useful for the output of programs like SnugDock, for example, which requires all antigen chains to be labeled together as "A". It may prove useful those who want to score large numbers of models from disparate methods that number and label chains in varied ways. 

It is worth noting that our script can also compute just the underlying distance metrics (RMSDs, DOCKQ) with the -j flag. It is slightly quicker than DOCKQ.py, despite the fact that DOCKQ.py was explicitly given which chains should be paired to each other, while atom\_match\_n\_RD p-value.py was only told the chains for the docking partners in the true structure, and compared all chains to each other to automatically match them. So this script may also prove useful for those who would like a slight increase in speed, or to analyze models with chains that share the same identifier or have different chain naming between models of the same structure, when running many DOCKQ calculations en masse.





### SET UP ###

To set up RDP, please first set up Mosaist, located on GitHub at Grigoryanlab/Mosaist.

Then run RDP's makefile, providing the path to your version of g++ as CC, the path to Mosaist's include subdir as MST_INC_DIR, and the path to MST's objs dir as MST_OBJS_DIR. The version of g++ should be the same one as used by Mosaist. These can be provided as flags or filled into the makefile.

Example makefile command, with relevant paths given as flags: make CC=~/.conda/envs/myenv/bin/x86_64-conda_cos6-linux-gnu-g++ MST_INC_DIR=/path/to/Mosaist/include MST_OBJS_DIR=/path/to/Mosaist/objs

To run the backbone matching and RDP calculations, run the atom_match_n_RDP.py script. At minimum it requires the true structure (-a), the model structure (-b), the chains for the first docking partner in the crystal structure (-c1), the chains for the second docking partner in the crystal structure (-c2), the underlying distance metric for the distribution the RDP is taken from (-m), and the out path (-o). By default it uses complex RMSD as the distance metric for the Random Docking P-value. For details on how to format the input for each flag and all its other options, run the script with the -h flag. Here is an example command: python3 atom_match_n_RDP.py -a /path/to/true/structure/1ahw.pdb -b /path/to/model/structure/1ahw_0000.pdb -a1 'A,B' -a2 'C' -m r -o out/path/testRDP.txt 
