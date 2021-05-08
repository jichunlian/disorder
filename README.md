## This branch is a beta version of ***disorder***
### New tags
- ***cmax-tag:*** Specifies the maximum number of configurations allowed to output. When ***cmax*** = 0 or >= ***N*** (***N*** is the number of irreducible configurations), all irreducible configurations will be output, otherwise only ***cmax*** of them will be output randomly. (**Default:** ***cmax*** **= 0**)
- ***lsep-tag:*** Specifies whether the structure of every irreducible configuration is write to the files **POSCAR-x...** (**TRUE**) or whether it is merged (**FALSE**) to a single file **POSCAR**. (**Default:** ***lsep*** **= .true.**)
