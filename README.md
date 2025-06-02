# README for MATLAB Scripts

Kulvicius T, Tamosiunaite M, Ainge J, Dudchenko P, Worgotter F (2008)
Odor supported place cell model and goal navigation in rodents.
*J Comput Neurosci* 25:481-500

## MATLAB Scripts

- **navPFdevelv03.m**
  Place field model implemented by using feed-forward network with winner takes all learning algorithm. Place cells are formed from visual and olfactory input.

- **createCellsAll.m**
  Creates place fields. This script is called by *navPFdevelv03.m*.

- **plotPFM.m**
  Plots place field maps. This function is called by *navPFdevelv03.m*.

- **navPFQLv14.m**
  Goal navigation based on place cells and Q-learning with function approximation.

- **navUrineBasedv04.m**
  Goal navigation based on self-generated odor marks (no place fields).

- **navPFQLUv07.m**
  Goal navigation using combined navigation algorithm (Q-learning based on place cells + self-marking navigation).

These files were supplied by Tomas Kulvicius.

---

2025-06-02: Converted README to Markdown.