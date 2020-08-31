# 1. About
## 1.1. What is ***`disorder`*** ?
- ***`disorder`*** is an open source software designed for generating irreducible site-occupancy configurations (i.e., symmetrically inequivalent disordered structures), which can be used for disordered doping, including substitution doping and vacancy doping.

## 1.2. Features
- Build a supercell form arbitrary input cell with any supercell expansion matrix.
- Search space group operations of the supercell and identify its point group symbol.
- Construct the equivalent atomic matrix of the atomic positions for disordered doping.
- Generate irreducible configurations for any number of atomic types with arbitrary stoichiometry.
- Count the degeneracy (the number of equivalent configurations) of each irreducible configuration.
- Output the **POSCAR** file (the input file of **`VASP`**) of each irreducible configuration.

## 1.3. License
- Copyright © 2020 Ji-Chun Lian
- ***`disorder`*** is licensed under the **MIT License**, See the **LICENSE** file for license rights and limitations.

## 1.4. Author & Contact
- Ji-Chun Lian (jichunlian@hnu.edu.cn), Department of Applied Physics, School of Physics and Electronics, Hunan University, Changsha, China
- If you have any questions, suggestions, and problems regarding ***`disorder`***, please feel free to contact me.


# 2. Compiling
- The default compiler is ifort.
- Please modify the Makefile in the **src/** folder, if you want to use another compiler.
- Make compile the executables into the **bin/** folder:

```
make disorder     # Only disorder  is compiled
make supercell    # Only supercell is compiled
make              # Both disorder and supercell are compiled
make clean        # rm *.mod *.o
```

# 3. Input Files

### 3.1. SPOSCAR (Crystal structure file)
- The format of **SPOSCAR** is consistent with the **POSCAR** of **`VASP.5.x`**.
- **SPOSCAR** should be a supercell, because ***`disorder`*** does not expand the cell inputted from SPOSCAR.
- We provide a util (i.e., ***`supercell`***) for expanding the small cell.
- The non-diagonal supercell expansion matrix is also supported in ***`supercell`***.
- More detailed information about ***`supercell`*** can be seen at the end of this document.
- Any other software can also be used for cell expansion, as long as the format of **SPOSCAR** is correct.


### 3.2. INDSOD (Running control file)
- **INDSOD** contains **11** parameters to control the running of ***`disorder`***.
- Its format and the meaning of all parameters are shown below:

```
&input                              # The title of namelist
  nsub = integer                    # nsub-nary substitution (Default 2)
  subs = integer,integer,...        # Atomic numbers, nsub integers (No default)
  symb = character,character,...    # Atomic symbols, nsub integers (No default)
  site = integer                    # Atomic species to be substituted (Default 1)
  prec = real                       # Symmetry precision (Default 1D-4)
  fast = logical                    # Whether to turn on Fast Mode (Semi-automatic)
  lpro = logical                    # Whether to display the progress bar (Default .false.)
  lpos = logical                    # Whether to output the POSCAR-x... files (Default .false.)
  leqa = logical                    # Whether to output the EQAMAT file (Default .false.)
  lspg = logical                    # Whether to output the SPGMAT file (Default .true.)
  lcfg = logical                    # Whether to output the CFGMAT file (Default .true.)
/                                   # Terminator
```
>**Note:**  
> 1. **&input** and **/** are indispensable and immutable.  
> 2. The order of the above parameters is irrelevant.  
> 3. The exclamation mark **!** can be used for annotations.  
> 4. The atomic symbol of **Kw** represents the vacancy.



# 4. Output Files

### 4.1. SPGMAT (Space group operations file)
### 4.2. EQAMAT (Equivalent atomic matrix file)
### 4.3. CFGMAT (Irreducible atomic configurations file)
### 4.4. POSCAR (Irreducible disordered structure files)



# 5. Running
- We provide two examples for testing (i.e., SnxPb1-xTe and TiO2-VO).
- We use the first example to introduce the running of ***`disorder`***.
- The second example is used to introduce the running of ***`supercell`***.

## 5.1. Program ***`disorder`***

```
cd examples/1_SnxPb1-xTe/
../../bin/disorder
```

## 5.2. Program ***`supercell`***

```
cd examples/2_TiO2-VO/
```

- **If the small cell is an unit cell**

```
cp POSCAR_unit POSCAR
../../bin/supercell
  3  0  0
  0  3  0
  0  0  1
```
- **If the small cell is a primitive cell**

```
cp POSCAR_prim POSCAR
../../bin/supercell
  3  0  3
  0  3  3
 -1 -1  0
```
