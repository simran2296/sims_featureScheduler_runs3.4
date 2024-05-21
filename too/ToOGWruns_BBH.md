**Simulate BBH ToOs**

**SCHEMA (different from earlier simulations)**:

1) simulate 3xNtriggers where Ntrigger is the number of recommended triggers as all-sky events.
2) follow up any of the area that is visible, regardless of the location of the trigger
3) simulate
- 1 run with 3 x 3/4 x Ntrigger
- 1 run with 3 x 4/5 x Ntrigger
- 1 run with 3 x Ntrigger
- 1 run with 3 x 6/5 x Ntrigger
- 1 run with 3 x 5/4 x Ntrigger
- 1 run with 3 x Ntrigger with 30sec exposures (e.g. where 120sec exposures are asked simulate 4x30s) QUESTION: can this run use the same triggers as the 3 x Ntrigger run so that we can measure the impact of exposure length aside from other stochastic changes?


**Ntriggers = 7** triggers

All after June 2027 and before January 2030

Case A) 9 triggers with 10 sq deg area -> yields 3 observable

Case B) 6 triggers with 20 sq deg area -> yields 2 observable

Case C) 6 triggers with 30 sq deg area -> yields 2 observable


**A/B/C** (7)

if bright: 
            filters = (u)g(r)i 

else: 
            filters = riz

* Night 0,2,7,9,39:
  
            1x[filters] x 1 pass 30sec
