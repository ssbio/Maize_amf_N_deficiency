*************************************************************
*FVA model for Whole Maize
*************************************************************
*****************Nicholas Ampimah************************
*************************************************************
$INLINECOM /*  */
$offdigit

OPTIONS

        limrow = 10000
        limcol = 10000
        optCR = 1E-9
        optCA = 0.0
        iterlim = 100000
        decimals = 8
        reslim = 100000
        work = 50000000
;

*********Defining Sets**************************************
SETS

	i					set of metabolites

$include "metabolites.txt"	

	j					set of reactions

$include "reactions.txt"


;

alias(j,j1);

*************************************************************

***********Defining Parameters*******************************
PARAMETERS

	S(i,j)					stoichiometric matrix

$include "sij.txt"

	v_max(j)				maximum flux of v(j)
	
$include "LO3_LN_CONTROL_upper_bound.txt"

	v_min(j)				minimum flux of v(j)

$include "LO3_LN_CONTROL_lower_bound.txt"

	c(j)
;

**************************************************************

*********Defining Equations***********************************
EQUATIONS

	objective				objective function
	mass_balance1(i)			steady state mass balance
	
	lower_bound(j)				lower bounds on reactions
	upper_bound(j)				upper bounds on reactions
	photorespiration			setting maximum rate of photorespiration
	leafmalateshunt				shutting down leaf malate shunt
	shoot_CO2_out				shoot can survive without photosynthesis
	biomass_enforcement_1(j)	biomass enforcement_1
	biomass_enforcement_2(j)	biomass enforcement_2
	biomass_enforcement_3(j)	biomass enforcement_3
	biomass_enforcement_4(j)	biomass enforcement_4
;
**************************************************************

*********Defining Variables***********************************
FREE VARIABLES

	v(j)					reaction flux
	z					objective value
;

****************************************************************

***************Defining Model***********************************

objective..			z =e= sum(j,c(j)*v(j));

mass_balance1(i)..		sum(j, S(i,j) * v(j)) =e= 0;


lower_bound(j)..		v_min(j) =l= v(j);

upper_bound(j)..		v(j) =l= v_max(j);

photorespiration..
    v('R03140[B,p]') =e= 0.0003 * v('R00024[B,p]');

leafmalateshunt..
    v('R00342[M,c]') =e= 0;

shoot_CO2_out..
    v('OUT_C00011[S,c]') =g= 0;

biomass_enforcement_1(j)..
    v('Root_Biomass[R]')  =g= 2.131;

biomass_enforcement_2(j)..
    v('Shoot_Biomass[S]') =g= 0.532;

biomass_enforcement_3(j)..
    v('Seed_Biomass[K]')  =g= 1.598;

biomass_enforcement_4(j)..
    v('Leaf_Biomass[L]')  =g= 1.065;


;

******************************************************************

**********Solving Model*********************

Model fva /all/;

fva.optfile = 1;

fva.holdfixed = 1;

FILE RESULTS /LO3_LN_CONTROL_FVA.txt/;

PUT RESULTS;



PUT "reaction     min_rate    max_rate"/;

LOOP(j1,

      c(j) = 0;
      c(j1) = 1;
      
      PUT j1.tl:0:100;

      solve fva using lp MINIMIZING z;
      
      PUT z.l:20:5;
      
      solve fva using lp MAXIMIZING z;
      
      PUT z.l:20:5/;
      
);

PUTCLOSE;
**********************************************   

