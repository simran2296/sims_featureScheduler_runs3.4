**Simulate BBH and NBH ToOs** (for lensed BNS and BBH see https://github.com/fedhere/_scoc_/blob/main/ToOGWruns_lensedBNS.md and https://github.com/fedhere/_scoc_/blob/main/ToOGWruns_BBH.md)

**SCHEMA (different from earlier simulations)**: 
1) simulate 3xNtriggers where Ntrigger is the number of recommended triggers as all-sky events.
2) follow up any of the area that is visible, regardless of the location of the trigger
3) simulate 
- 1 run with 3 x 1/2 x Ntrigger
- 1 run with 3 x 3/4 x Ntrigger
- 1 run with 3 x Ntrigger
- 1 run with 3 x 3/2 x Ntrigger
- 1 run with 3 x 2 x Ntrigger
- 1 run with 3 x Ntrigger with 30sec exposures (e.g. where 120sec exposures are asked simulate 4x30s) _QUESTION: can this run use the same triggers as the 3 x Ntrigger run so that we can measure the impact of exposure length aside from other stochastic changes?_

**Ntrigger = 22** triggers

All after June 2027 and before January 2030 simulate trigger time at random over 24h.

Case A) 9 triggers with 30 sq deg area -> should yield 3 visible. 

Case B) 15 triggers with 50 sq deg area -> should yield 5 observable. 

Case C) 24 triggers with 100 sq deg area -> should yield 8 observable. 

Case D) 6 triggers with 250 sq deg area -> should yield 4 observable. 

Case E) 12 triggers with 500 sq deg area -> should yield 2 observable. 



define Tw observability window as the time from the trigger to the end of visibility of a field in the area (_QUESTION: what is the max airmass we observe at??_)

**A:** (3x3)
*  Night 0:
      
        if Tw> 4.7h:

           [(u)griz(y)] x 3 passes 120sec
   
        if 3.2h < Tw < 4.7h:
   
           [(u)griz(y)] x 2 passes 120sec
   
        if Tw < 3.2h:
   
           [(u)griz(y)] x 1 passes 120sec
* Night 1,2:

        [(u)griz(y)] x 1 pass 120sec
* Night 3 only 1/3 of the times:

        [(u)griz(y)] x 1 pass 120sec
  

**B/C:** (3x13)
* Night 0:

      if Tw> 4.7h:

        [gri] x 3 passes 120sec

      else if 3.2h < Tw < 4.7h:

        [gri] x 3 passes 120sec

      else if Tw < 3.2h:

        [gri] x 1 pass 120sec + [4xr] x 1 pass 120sec
* Night 1,2:

      [ri] x 1 pass 180sec

* Night 3 only 1/3 of the times:

      [ri] x 1 pass 180sec
      

**D/E:** (3x6)
* Night 0:

        [gr] x 1 pass
  
* Night 1,2:
  
        [gr] x 1 pass 120sec
 
* Night 3 only 1/3 of the times:
  
        [gr] x 1 pass 120sec


