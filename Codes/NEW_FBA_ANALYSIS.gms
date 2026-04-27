*************************************************************
*Basic FBA model for Maize_Whole_Plant
*************************************************************
*****************Nicholas Amppimah************************
*************************************************************
$INLINECOM /*  */
$offdigit

OPTIONS

        decimals = 8
	lp = cplex
;

*********Defining Sets**************************************
SETS

	i					set of metabolites

$include "metabolites.txt"	

	j					set of reactions

$include "reactions.txt"


;
*************************************************************

***********Defining Parameters*******************************
PARAMETERS

	S(i,j)					stoichiometric matrix

$include "sij.txt"

	v_max(j)				maximum flux of v(j)
	
$include "ANOTHER_upper_bound_80.txt"

	v_min(j)				minimum flux of v(j)

$include "ANOTHER_lower_bound_80.txt"




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
	obj					objective value
;

****************************************************************

***************Defining Model***********************************

objective..			obj =e= v('Seed_Biomass[K]') + v('Root_Biomass[R]') + v('Shoot_Biomass[S]') +  v('Leaf_Biomass[L]') ;



mass_balance1(i)..		sum(j, S(i,j) * v(j)) =e= 0;


lower_bound(j)..		v_min(j) =l= v(j);

upper_bound(j)..		v(j) =l= v_max(j);

photorespiration..		v('R03140[B,p]') =e= 0.0003 * v('R00024[B,p]');

leafmalateshunt..		v('R00342[M,c]') =e= 0;

shoot_CO2_out..			v('OUT_C00011[S,c]') =g= 0;

biomass_enforcement_1(j)..	v('Root_Biomass[R]') =g= 0.4* obj;

biomass_enforcement_2(j)..  v('Shoot_Biomass[S]') =g= 0.1* obj;

biomass_enforcement_3(j)..  v('Seed_Biomass[K]') =g= 0.3* obj;

biomass_enforcement_4(j)..  v('Leaf_Biomass[L]') =g= 0.2*obj;


***************Set bounds***********************************


Model arabidopsis_root /all/;

******************************************************************

**********Solving Model*********************

Solve arabidopsis_root using lp maximizing obj;

arabidopsis_root.holdfixed = 1;


********************************************
****************Output File*****************

FILE RESULTS /FBA_Results.txt/;

PUT RESULTS;

PUT "The max Biomass value is : " obj.l:0:8//;

loop(j,
	put j.tl:0:30, "    ", v.l(j):0:8/;
);

PUTCLOSE;
**********************************************