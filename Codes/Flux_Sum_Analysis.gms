*************************************************************
****************** FLUX SUM ANALYSIS ************************
***************  MIP + TSV output (updated) ******************
*************************************************************

$INLINECOM /*  */
$offdigit
$eolCom  #
$onempty

* solver settings
option mip = cplex;
*option mip = baron;
*$onecho > baron.opt
*Threads 16
*$offecho

* ---- Quiet run / clean logs
option decimals = 8;
option solprint = off;
option sysout   = off;
option limrow   = 0, limcol = 0;

*************************************************************
******************* Sets ************************************
*************************************************************

SETS
    i           "metabolites"
$include "metabolites.txt"

    SM(i)       "selected metabolites"
*$include "single_metabolite_test.txt"*
$include "NEW_selected_metabolites.txt"

    j           "reactions"
$include "reactions.txt"
;

alias(i,i1);

SET active(i) "currently selected metabolite";

*************************************************************
******************* Parameters ******************************
*************************************************************

PARAMETERS
    S(i,j)      "stoichiometric matrix"
$include "sij.txt"

    v_max(j)    "upper bounds"
$include "ANOTHER_upper_bound_80.txt"

    v_min(j)    "lower bounds"
$include "ANOTHER_lower_bound_80.txt"
;

*************************************************************
******************* Scalars *********************************
*************************************************************

SCALARS
    M   "big-M constant for abs linearization" /100000/;
;

*************************************************************
******************* Variables *******************************
*************************************************************

FREE VARIABLES
    v(j)        "reaction flux"
    Z           "objective";

POSITIVE VARIABLES
    f_plus(i,j) "positive part of S*v"
    f_minus(i,j)"negative part of S*v";

BINARY VARIABLES
    i_plus(i,j)
    i_minus(i,j);

*************************************************************
******************* Equations *******************************
*************************************************************

EQUATIONS
    objective
    mass_balance(i)
    rearrangement_1(i,j)
    rearrangement_2(i,j)
    rearrangement_3(i,j)
    rearrangement_4(i,j)
    lower_bound(j)
    upper_bound(j)
    photorespiration
    leafmalateshunt
    shoot_CO2_out
    biomass_enforcement_1
    biomass_enforcement_2
    biomass_enforcement_3
    biomass_enforcement_4
;

*************************************************************
******************* Objective *******************************
*************************************************************

objective..
    Z =e= 0.5 * sum((i,j)$active(i), f_plus(i,j) + f_minus(i,j));

*************************************************************
******************* Mass Balance ****************************
*************************************************************

mass_balance(i)..
    sum(j, S(i,j) * v(j)) =e= 0;

*************************************************************
******** Absolute Value Linearization (MIP) ******************
*************************************************************
* Optional sparsity filter: add $(S(i,j)<>0) if your S is sparse.

rearrangement_1(i,j)$active(i)..
    S(i,j) * v(j) =e= f_plus(i,j) - f_minus(i,j);

rearrangement_2(i,j)$active(i)..
    f_plus(i,j) =l= i_plus(i,j) * M;

rearrangement_3(i,j)$active(i)..
    f_minus(i,j) =l= i_minus(i,j) * M;

rearrangement_4(i,j)$active(i)..
    i_plus(i,j) + i_minus(i,j) =e= 1;

*************************************************************
******************* Flux Bounds *****************************
*************************************************************

lower_bound(j)..
    v_min(j) =l= v(j);

upper_bound(j)..
    v(j) =l= v_max(j);

*************************************************************
************** Additional Biological Constraints ************
*************************************************************

photorespiration..
    v('R03140[B,p]') =e= 0.0003 * v('R00024[B,p]');

leafmalateshunt..
    v('R00342[M,c]') =e= 0;

shoot_CO2_out..
    v('OUT_C00011[S,c]') =g= 0;

biomass_enforcement_1..
    v('Root_Biomass[R]')  =g= 2.885;

biomass_enforcement_2..
    v('Shoot_Biomass[S]') =g= 0.721;

biomass_enforcement_3..
    v('Seed_Biomass[K]')  =g= 2.164;

biomass_enforcement_4..
    v('Leaf_Biomass[L]')  =g= 1.442;

*************************************************************
******************* Model ***********************************
*************************************************************

Model flux_sum_analysis /all/;

flux_sum_analysis.optfile=1;

*************************************************************
******************* TSV Output ******************************
*************************************************************

File results /Nucl_ROOT_Palmitoleic_FSA_RESULTS_debug.tsv/;
Put results;

Put
    "Metabolite"        @25,
    "MaxFluxSum"        @40,
    "MinFluxSum"        @55,
    "MaxSolveStat"      @70,
    "MaxModelStat"      @85,
    "MinSolveStat"      @100,
    "MinModelStat"      @115,
    "MaxRootBio"        @130,
    "MaxShootBio"       @145,
    "MaxSeedBio"        @160,
    "MaxLeafBio"        @175,
    "MinRootBio"        @190,
    "MinShootBio"       @205,
    "MinSeedBio"        @220,
    "MinLeafBio"        @235 /;


*************************************************************
*********************** Iteration ***************************
*************************************************************

scalar maxZ, minZ, maxSS, maxMS, minSS, minMS;
scalar
    maxRootBio, maxShootBio, maxSeedBio, maxLeafBio,
    minRootBio, minShootBio, minSeedBio, minLeafBio;

scalar nSM, k;
nSM = card(SM);
k   = 0;

loop(i1$SM(i1),
	k = k + 1;
    put_utility 'log' / 'FSA progress: ', k:0:0, '/', nSM:0:0, '  ', i1.tl:0:0;

    active(i)  = no;
    active(i1) = yes;

* Maximize flux-sum (MIP)*
    solve flux_sum_analysis using mip maximizing Z;
    maxZ  = Z.l;
    maxSS = flux_sum_analysis.solvestat;
    maxMS = flux_sum_analysis.modelstat;
	maxRootBio  = v.l('Root_Biomass[R]');
	maxShootBio = v.l('Shoot_Biomass[S]');
	maxSeedBio  = v.l('Seed_Biomass[K]');
	maxLeafBio  = v.l('Leaf_Biomass[L]');


* Minimize flux-sum (MIP)*
    solve flux_sum_analysis using mip minimizing Z;
    minZ  = Z.l;
    minSS = flux_sum_analysis.solvestat;
    minMS = flux_sum_analysis.modelstat;
	minRootBio  = v.l('Root_Biomass[R]');
	minShootBio = v.l('Shoot_Biomass[S]');
	minSeedBio  = v.l('Seed_Biomass[K]');
	minLeafBio  = v.l('Leaf_Biomass[L]');


    Put
		i1.tl            @25,
		maxZ:0:8         @40,
		minZ:0:8         @55,
		maxSS:0:0        @70,
		maxMS:0:0        @85,
		minSS:0:0        @100,
		minMS:0:0        @115,
		maxRootBio:0:8   @130,
		maxShootBio:0:8  @145,
		maxSeedBio:0:8   @160,
		maxLeafBio:0:8   @175,
		minRootBio:0:8   @190,
		minShootBio:0:8  @205,
		minSeedBio:0:8   @220,
		minLeafBio:0:8   @235 /;



);

PutClose;
