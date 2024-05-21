

**Simulate neutrino ToO**
only observe if Galactic Latitude > |10 deg|

**SCHEMA (different from earlier simulations):**

simulate 4xNtriggers where Ntrigger is the number of recommended triggers as all-sky events.
follow up any of the area that is visible **and Galactic Latitude >|10 deg|** _QUESTION: I think the factor should be 4 not three if the constrain is latitude >|10 deg|)


**Ntrigger** = 40 (4/year) triggers

**Area** = 7 sq deg (one pointing)
t0 = time of trigger

* Night 0: 

        [g] @ t=t0 120sec

        [r] @ t=t0 + 15 minutes

        [z] within the night

* Night 1: 

        [g] @ t=t0 120s

        [r] @ t=t0 + 15 minutes

* 6 <= Night <= 8:

        [g] @ t=t0 120sec

        [r] @ t=t0 + 15 minutes

* 0 <= Night <= 60: first night with u-band filter after t0
    
        [u]
