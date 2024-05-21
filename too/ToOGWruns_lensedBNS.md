**Simulate lensed BNS ToOs**

**SCHEMA (different from earlier simulations and different from the other GW cases)**:

1) simulate Ntriggers with center in the visible sky.
2) follow up any of the area that is visible, regardless of the location of the trigger
3) simulate
- 1 run with 3 x Ntrigger
- 1 run with 3 x Ntrigger with 30sec exposures (e.g. where 120sec exposures are asked simulate 4x30s) QUESTION: can this run use the same triggers as the 3 x Ntrigger run so that we can measure the impact of exposure length aside from other stochastic changes?


**Ntriggers = 2** triggers

All after June 2027 and before January 2030

Case A) 1 trigger 900 sq deg in observable sky

Case B) 1 trigger 15 sq deg in observable sky

**A:** (1)

* Night >= 0:
  
          [g] x 1 pass 30sec
          [r] x 1 pass 90sec

**B:** (1)
  * Night >= 0: _QUESTION: How do we deal with this? clearly we cannot do a 1500 seconds integration, we need to break it up. I think break it up in 150sec exposures
  
          10x[g] + 10x[r] x 1 pass 150sec 

