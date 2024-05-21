**Simulate neutrino ToO**
only simulate Galactic Latitude > 10 deg

**N** = 40 (4/year)

**Area** = 7 sq deg (one pointing)
t0 = time of trigger

Night 0: 

*  4xg @ t=t0
  
*  1xr @ t=t0 + 15 minutes
  
*  1xz within the night

Night 1: 
  
*  4xg @ t=t0
  
*  1xr @ t=t0 + 15 minutes

6 <= Night <= 8:
  
*  4xg @ t=t0
  
*  1xr @ t=t0 + 15 minutes

0 <= Night <= 60: first night with u-band filter after t0
*   1xu

**Simulate SSO ToO**

**N** = 300 (30/year)

**Area**: 10 sq deg

Night 0:
*   if not twilight:

      2xr x 3 passes separated by 33 minutes with dither between the two 30s observations
*   else:

    2xz @ 15sec x 3 passes separated by 33 minutes with dither between the two 15s observations

