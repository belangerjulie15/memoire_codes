# memoire_codes
Codes used for the thesis: Portfolio optimisation applied for variable annuities


This file contains many R codes:

-Memoire_Codes : which contains basic formulae and functions (which we can find in the package Optimisation. Power. Utility), and codes used to plot the graphs in the thesis.

-Memoire_Codes1 : which contains functions used to plot relativity between variable (i.e. ptf value wrt to S_t, xi wrt S_t, ..). It also contains functions that Optimisation.Power.Utility contains (such as, find_lambda_optimal)

-Memoire_Codes2 : which contains codes to plot relations between variables whit and without taking into account fees. For example, Xi_t wrt to S_t (with and without fees). Also contains Histogramme Portefeuille optimal: dynamique et forme fermée:  to compare the results obtained while similating directly the optimal ptf at maturity and with rebalancing at each simulation node. 

-Memoire_Codes3.3: which contains the basic function while optimising a portfolio. It can return, the average of the utility of the mutual funds, the average of the utility of the final portfolio, the probability that the guaranty is exercised, the martingale verification value. 

-Memoire_Codes4 : which contains excatly the same elements of Memoire_Codes3.3 but with a parallelized code (for Windows and Mac).

-Memoire_Codes5 : which contains the functions used to find the Optimal Fee with the 4 methods (Méthode Martingale, MM Bornée [0,1], proportion cte -au choix- et par simulation du portfeuille optimal final directement).

-Memoire_Codes6 : which contains all the functions used for mutual funds optimisation for MANY insured.
