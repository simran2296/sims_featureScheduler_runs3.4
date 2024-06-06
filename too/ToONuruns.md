

**Simulate neutrino ToO**
only observe if Galactic Latitude > |10 deg|

**Definition**:

**N_O** = N Observable Triggers is the number of expected number of O5 all-sky events that will meet the selection cuts established in [https://lssttooworkshop.github.io/images/Rubin_2024_ToO_workshop_final_report.pdf](https://lssttooworkshop.github.io/images/Rubin_2024_ToO_workshop_final_report.pdf).

**N_t** = N Triggers is the number of all-sky alerts (produced by LKV or any other discovery survey) that meet our selection cuts, regardless of location

**Night 0**: the first night after the trigger when the sky location is observable and our follow-up starts. _NOTE_: we should start with night 0 follow up as soon as possible and no later than 48 hours from the time of the LKV alert (could be same as night _of_ trigger)

**t0**: the time when the trigger is received by Rubin

**Tw**: observability window as the time from the trigger to the end of visibility of a field in the area 


**SCHEMA :**

- simulate N_t=4xN_O
- follow up any of the trigger location is visible **and Galactic Latitude >|10 deg|** _Note: the factor should be 4 accounts for visibility (x3) and additional constraint latitude >|10 deg|)


**N_t** = 160 (16/year) all sky triggers

**Area** = 7 sq deg (one pointing)

* Night 0: 

        [g] @ t=t0 1 pass 120sec 

        [r] @ t=t0 +> 15 minutes 1 pass 30sec

        [z] within the night 1 pass 30sec

* Night 1: 

        [g] @ t=t0 1 pass 120s at some point in the night

        [r] @ t=t0 1 pass 30sec right after

* 6 <= Night <= 8:

        [g] x 1 pass 30s at some point in the night

        [r] x 1 pass 30s right after

        [z] x 1 pass 30s right after

* 0 <= Night <= 60: first night with u-band filter after t0
    
        [u] x 1 pass 30 s

_Expected Yield N_0=40_
