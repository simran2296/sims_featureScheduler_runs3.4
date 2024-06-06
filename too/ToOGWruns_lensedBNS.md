**Simulate lensed BNS ToOs**


**Definition**:

**N_O** = N Observable Triggers is the number of expected number of O5 all-sky events that will meet the selection cuts established in [https://lssttooworkshop.github.io/images/Rubin_2024_ToO_workshop_final_report.pdf](https://lssttooworkshop.github.io/images/Rubin_2024_ToO_workshop_final_report.pdf).

**N_t** = N Triggers is the number of all-sky alerts (produced by LKV or any other discovery survey) that meet our selection cuts, regardless of location

**Night 0**: the first night after the trigger when the sky location is observable and our follow-up starts. _NOTE_: we should start with night 0 follow up as soon as possible and no later than 48 hours from the time of the LKV alert (could be same as night _of_ trigger)

**t0**: the time when the trigger is received by Rubin

**Tw**: observability window as the time from the trigger to the end of visibility of a field in the area 


**SCHEMA**: 

1) simulate N_t = 3xN_O

3) follow up if location of the trigger is visible

3) simulate
- 1 run with N_t
- 1 run with N_t with 30sec exposures (e.g. where 120sec exposures are asked simulate 4x30s) QUESTION: can this run use the same triggers as the 3 x Ntrigger run so that we can measure the impact of exposure length aside from other stochastic changes?


**N_t=6** triggers

All after June 2027 and before January 2030

Case A) N_t=3 trigger 900 sq deg in observable sky

Case B) N_t=3trigger 15 sq deg in observable sky

**A:** (1)

* Night >= 0:
  
          [g] x 1 pass 30sec
          [r] x 1 pass 90sec

**B:** (1)
  * Night >= 0: 
  
          10x[g] + 10x[r] x 1 pass 150sec 

