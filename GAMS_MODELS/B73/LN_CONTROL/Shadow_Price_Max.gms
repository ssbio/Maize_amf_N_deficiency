*************************************************************
* Shadow Price Analysis for Whole_Maize
*************************************************************
*****************Nicholas Ampimah****************************
*************************************************************

$INLINECOM /*  */
$offdigit

OPTIONS

        lp       = cplex
        limrow   = 10000
        limcol   = 10000
        optCR    = 1E-9
        optCA    = 0.0
        iterlim  = 100000
        decimals = 8
        reslim   = 100000
        work     = 50000000
;

*************************************************************
**********************Defining Sets**************************

SETS

    i       Set of metabolites

$include "metabolites.txt"

    j       Set of reactions

$include "reactions.txt"

;

*************************************************************
********************Defining Parameters**********************

PARAMETERS

    S(i,j)          Stoichiometric matrix

$include "sij.txt"

    v_max(j)        Maximum reaction flux

$include "B73_LN_CONTROL_upper_bound.txt"

    v_min(j)        Minimum reaction flux

$include "B73_LN_CONTROL_lower_bound.txt"

    
    shadow_max(i)   Shadow price at maximum objective
;

*************************************************************
********************Defining Equations***********************

EQUATIONS

    objective               Objective function
    mass_balance1(i)         Steady-state mass balance
    lower_bound(j)           Lower reaction bounds
    upper_bound(j)           Upper reaction bounds
    photorespiration         Setting photorespiration rate
    leafmalateshunt          Shutting down leaf malate shunt
    shoot_CO2_out            Shoot CO2 export constraint
    biomass_enforcement_1    Root biomass enforcement
    biomass_enforcement_2    Shoot biomass enforcement
    biomass_enforcement_3    Seed biomass enforcement
    biomass_enforcement_4    Leaf biomass enforcement
;

*************************************************************
********************Defining Variables***********************

FREE VARIABLES

    v(j)    "Reaction flux"
    obj     "Objective value"
;

*************************************************************
*********************Model Equations*************************

objective..

    obj =e= v('Seed_Biomass[K]') + v('Root_Biomass[R]') + v('Shoot_Biomass[S]') + v('Leaf_Biomass[L]');


mass_balance1(i)..

    sum(j, S(i,j) * v(j)) =e= 0;


lower_bound(j)..

    v_min(j) =l= v(j);


upper_bound(j)..

    v(j) =l= v_max(j);


photorespiration..

    v('R03140[B,p]') =e=
        0.0003 * v('R00024[B,p]');


leafmalateshunt..

    v('R00342[M,c]') =e= 0;


shoot_CO2_out..

    v('OUT_C00011[S,c]') =g= 0;


biomass_enforcement_1..

    v('Root_Biomass[R]') =g= 4.188;


biomass_enforcement_2..

    v('Shoot_Biomass[S]') =g= 1.560;


biomass_enforcement_3..

    v('Seed_Biomass[K]') =g= 0;


biomass_enforcement_4..

    v('Leaf_Biomass[L]') =g= 4.722;

*************************************************************
*********************Defining Model**************************

MODEL Whole_Maize /all/;

Whole_Maize.optfile   = 1;
Whole_Maize.holdfixed = 1;

*************************************************************

********Shadow Prices at Maximum Objective******************

SOLVE Whole_Maize USING LP MAXIMIZING obj;

ABORT$(Whole_Maize.solvestat <> 1)
    "Solver failed during objective maximization",
    Whole_Maize.solvestat,
    Whole_Maize.modelstat;

ABORT$((Whole_Maize.modelstat <> 1) AND
       (Whole_Maize.modelstat <> 2))
    "No optimal solution found during objective maximization",
    Whole_Maize.modelstat;

shadow_max(i) = mass_balance1.m(i);

*************************************************************
****************Shadow-Price Output File*********************

FILE RESULTS /B73_LN_CONTROL_shadow_prices.txt/;

PUT RESULTS;

PUT
    "Metabolite":30,
    "Shadow_price_max":20
/;

LOOP(i,

    PUT
        i.tl:30,
        shadow_max(i):20:8
    /;
);

PUTCLOSE RESULTS;

*************************************************************
*********************Display Results*************************

DISPLAY shadow_max;

*************************************************************