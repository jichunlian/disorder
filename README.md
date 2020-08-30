# Compiling
- The default compiler is ifort.
- Please modify the Makefile in the src/ folder, if you want to use another compiler.
- Make compile the executables into the bin/ folder:
> make disorder  # Only disorder  is compiled\
> make supercell # Only supercell is compiled\
> make           # Both disorder and supercell are compiled\
> make clean     # rm *.mod *.o



# Input Files

## SPOSCAR (Crystal structure file)
- The format of SPOSCAR is consistent with the POSCAR of VASP.5.x.
- SPOSCAR should be a supercell, because 'disorder' does not expand the cell inputted from SPOSCAR.
- We provide a util (i.e., 'supercell') for expanding the small cell.
- The non-diagonal supercell expansion matrix is also supported in 'supercell'.
- More detailed information about 'supercell' can be seen at the end of this document.
- Any other software can also be used for cell expansion, as long as the format of SPOCAR is correct.


## INDSOD (Running control file)
- INDSOD contains 11 parameters to control the running of 'disorder'.
- Its format and the meaning of all parameters are shown below:

&input                              # The title of namelist\
  nsub = integer                    # nsub-nary substitution (Default 2)\
  subs = integer,integer,...        # Atomic numbers, nsub integers (No default)\
  symb = character,character,...    # Atomic symbols, nsub integers (No default)\
  site = integer                    # Atomic species to be substituted (Default 1)\
  prec = real                       # Symmetry precision (Default 1D-4)\
  fast = logical                    # Whether to turn on Fast Mode (Semi-automatic)\
  lpro = logical                    # Whether to display the progress bar (Default .false.)\
  lpos = logical                    # Whether to output the POSCAR-x... files (Default .false.)\
  leqa = logical                    # Whether to output the EQAMAT file (Default .false.)\
  lspg = logical                    # Whether to output the SPGMAT file (Default .true.)\
  lcfg = logical                    # Whether to output the CFGMAT file (Default .true.)\
/                                   # Terminator\

- '&input' and '/' are indispensable and immutable.
- The order of the above parameters is irrelevant.
- The exclamation mark '!' can be used for annotations.
- The atomic symbol of 'Kw' represents the vacancy.



# Output Files

## SPGMAT (Space group operations file)
## EQAMAT (Equivalent atomic matrix file)
## CFGMAT (Irreducible atomic configurations file)
## POSCAR-x... (Irreducible disordered structure files)



# Running
- We provide two examples for testing (i.e., SnxPb1-xTe and TiO2-VO).
- We use the first example to introduce the running of 'disorder'.
- The second example is used to introduce the running of 'supercell'.

## disorder

> cd examples/1_SnxPb1-xTe/\
> ../../bin/disorder


## supercell

> cd examples/2_TiO2-VO/

### If the small cell is an unit cell
> cp POSCAR_unit POSCAR\
> ../../bin/supercell\
  3  0  0\
  0  3  0\
  0  0  1

### If the small cell is a primitive cell
> cp POSCAR_prim POSCAR\
> ../../bin/supercell\
  3  0  3\
  0  3  3\
 -1 -1  0
