## What is ***disorder*** ?
- ***disorder*** is an open source software designed for generating irreducible site-occupancy configurations (i.e., symmetrically inequivalent disordered crystal structures), which can be used for disordered doping, including substitution doping and vacancy doping.
- The ***disorder*** code works for arbitrary parent cells with any supercell expansion matrix, and for any number of atomic types with arbitrary stoichiometry. Most important, a linear scale of run time with the number of irreducible configurations is achieved, which is the best possible scaling for this type of problem.


## Features
- Build supercell form arbitrary input cell with any supercell expansion matrix.
- Search space group operations of the supercell and identify its point group symbol.
- Construct the equivalent atomic matrix of the atomic positions of the doping site.
- Generate irreducible configurations for any number of atomic types with arbitrary stoichiometry.
- Count the degeneracy (the number of equivalent configurations) of each irreducible configuration.
- Output the the crystal structure file of each irreducible configuration.

## License
- Copyright © 2020 Ji-Chun Lian
- ***disorder*** is licensed under the **MIT License**, See the **LICENSE** file for license rights and limitations.

## Author & Contact
- Ji-Chun Lian (jichunlian@hnu.edu.cn), Department of Applied Physics, School of Physics and Electronics, Hunan University, Changsha, China
- If you have any questions, suggestions, and problems regarding ***disorder***, please feel free to contact me.

## How to Cite ?
- Please cite the following article when you use ***disorder***：\
[1] J.-C. Lian, H.-Y. Wu, W.-Q. Huang, W. Hu, and G.-F. Huang, [Phys. Rev. B **102**, 134209 (2020)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.134209).
- If you want to understand the algorithm implemented in ***disorder***, you can read this reference carefully. After several updates, as an aside, the latest version of ***disorder*** possesses better time and space complexity than the original algorithm described in this reference, but the key concept is consistent.

## Future
- The code used for generating Special Quasirandom Structures (SQS) is under developing.
