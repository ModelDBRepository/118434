This is the readme for the matlab scripts associated with the paper:

Kulvicius T, Tamosiunaite M, Ainge J, Dudchenko P, Worgotter F (2008)
Odor supported place cell model and goal navigation in rodents.
J Comput Neurosci 25:481-500


navPFdevelv03.m - place field model implemented by using feed-forward
                   network with winner takes all learning
                   algorithm.  Place cells are formed from visual
                   and olfactory input.

createCellsAll.m - creates place fields. This script is called by
                   "navPFdevelv03.m".

plotPFM.m - plots place field maps. This function is called by
                   "navPFdevelv03.m".

navPFQLv14.m - goal navigation based on palce cells and Q-learning
                   with function aproximation.

navUrineBasedv04.m - goal navigation based on self generated odor
                   marks (no place fields).

navPFQLUv07.m - goal navigation using combined navigation algorithm
                   (Q-learning based on place cells + self-marking
                   navigation)

These files were supplied by Tomas Kulvicius.
