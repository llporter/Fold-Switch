Fold switchers search
Loren & Lauren
9.27.2016

1. Lauren processed the PDB sequence & DSSP output files down to almost the best we can do (missing just 128 chains from 56 PDB files- in Appendix I).
2. Got Òall_good_pdbs_092716.txtÓ file from her. Renamed it ss.txt.
3. Parsed it into the two separate files:
PDBSEARCHER_ParseSeqs -in ss.txt -out ss.seqs -ss ss.ss
4. Format the .seqs file for searching:
formatdb -p T -i ss.seqs
5. Divide up the directories:
PDBSEARCHER_DivideUpDirectories -in ss.seqs -size 1000
This made 298 directories out of the 297,968 PDB domains. It also makes a file called ÒPDBlist.txtÓ with all the domains in it.
6. Create mega-script that has all the individual calls to find fold-switchers:
PDBSEARCHER_MakeFindSwitcherScriptDirectories -in PDBlist.txt Ðsize 1000 -out bigscript.com
7. Divide up bigscript.com:
PDBSEARCHER_DivideUpFindSwitcherScript -in bigscript.com -size 1000
This made 298 scripts out of the 297,968 PDB domains.
8. Make a qsub script:
PDBSEARCHER_MakeQsubScript Ðsize 298
9. Execute the qsub script (after chmod u+x):
./qsubscript.com

-muts (number of mutations)






