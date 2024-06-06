**Simulate BBH and NBH ToOs** (for lensed BNS and BBH see [https://github.com/fedhere/_scoc_/blob/main/ToOGWruns_lensedBNS.md](https://github.com/fedhere/_scoc_/blob/main/ToOGWruns_lensedBNS.md) and [https://github.com/lsst-sims/sims_featureScheduler_runs3.4/blob/main/too/ToOGWruns_BBH.md](https://github.com/lsst-sims/sims_featureScheduler_runs3.4/blob/main/too/ToOGWruns_BBH.md))

**Definition**:

**N_O** = N Observable Triggers is the number of expected number of O5 all-sky events that will meet the selection cuts established in [https://lssttooworkshop.github.io/images/Rubin_2024_ToO_workshop_final_report.pdf](https://lssttooworkshop.github.io/images/Rubin_2024_ToO_workshop_final_report.pdf).

**N_t** = N Triggers is the number of all-sky alerts (produced by LKV or any other discovery survey) that meet our selection cuts, regardless of location

**Night 0**: the first night after the trigger when the sky location is observable and our follow-up starts. _NOTE_: we should start with night 0 follow up as soon as possible and no later than 48 hours from the time of the LKV alert (could be same as night _of_ trigger)

**t0**: the time when the trigger is received by Rubin

**Tw**: observability window as the time from the trigger to the end of visibility of a field in the area 


**SCHEMA**: 

1) simulate N_t = 3xN_O
   
3) follow up any of the area that is visible, regardless of the location of the trigger

4) simulate 
- 1 run with N=1/2 x N_t (downscope)
- 1 run with N=3/4 x N_t (downscope)
- 1 run with N=N_t
- 1 run with N=3/2 x N_t (upscope)
- 1 run with N=2 x N_t (upscope)
- 1 run with N_t with 30sec exposures (e.g. where 120sec exposures are asked simulate 4x30s) _NOTE_ : n this run should use the same triggers as the N=Ntrigger run so that we can measure the impact of exposure length aside from other stochastic changes

**N_t = 66** triggers simulated

All after June 2027 and before January 2030 simulate trigger time at random over 24h.

Case A) N_t = 9 with 30 sq deg area -> should yield N_O=3. 

Case B) N_t = 15 with 50 sq deg area -> should yield N_O=5. 

Case C) N_t = 24 with 100 sq deg area -> should yield N_O=8. 

Case D) N_t = 6 with 250 sq deg area -> should yield N_O=4. 

Case E) N_t = 12  with 500 sq deg area -> should yield  N_O=2. 

_Expected yield: N_O=22_

**Observing plans**

Initiate the observing sequence ASAP after trigger but no later than 48 hours after trigger

**A:** (3x3)
*  Night 0:
      
        if Tw > 4.7h:

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


