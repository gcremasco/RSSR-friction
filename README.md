SPATIAL  FOUR-BAR  MECHANISM OF RSSR TYPE WITH DYNAMIC FRICTION. POSITION 
AND FORCE ANALYSIS

Version: Beta1

Matlab  algorithm implementing the analytic solutions of the RSSR linkage 
with dynamic Coulomb friction. Official implementation used on paper: 

Cremasco  Coelho,  G.,  “Exact  solutions for the spatial  four-bar  with 
dynamic  friction:  Position,  force,  Painlevé paradox,  and singularity 
asymptotics”, Mech. Mach. Theory 142 (2025), 106163.

Implemented in Matlab R2018b.

Remark: 

- The  primary  goal  of this script was to produce numeric evidence that
  the solution on the paper is correct.  It  was  used to produce many of
  the figures shown there.
- The current algorithm was tested for some RSSR configurations. Particu-
  larly,  the inputs originally adopted here will reproduce results shown
  in  the  paper,  which  have  shown  agreement  with  ADAMS   multibody
  simulations. Nonetheless, NOT  ALL  RSSR  CONFIGURATIONS WERE TESTED in 
  this computer implementation.  In  other  words, more testing should be
  done to guarantee it works as a general-purpose RSSR solver.
- If  you  identify  errors,  please  notify the author. The email can be 
  found on the paper.
- A huge part of this m-file contains raw data from ADAMS simulations for
  the  RSSR  mechanism simulated in the article. This was used to compare
  exact  and  ADAMS  solutions.  Unless  you want these results, for some 
  reason, you can delete all the final lines indicated by the end of this
  script.

 Guilherme Cremasco Coelho
 sept-10th-2026
