
These are the matlab codes used in the simulation study of the 

Bai, J. "Panel data models with interactive fixed effects" Econometrica, 77, 1229-1279, 2009


The main program is called Mul_interactNew.m, which is still in the simulation mode. You need to replace the data generating process with actual data.  But first you may want to replicate Table 4 in my paper (in the online supplement of the Econometrica paper).

The rest five programs (*.m) are subroutines. There are also three output files when running the codes.

Also included is a utility program for generating the tables, called  Multable.m.  To replicate Table 4, set maxiter=100  (or 1000). When the computation is finished, just run

Multable(BETA), it will produce Table 4 in my paper.

Let me know if you have any problems by sending me an email at jb3064@columbia.edu

Jushan Bai
