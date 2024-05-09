**Simulate BBH and NBH ToOs**

_N total = 22_

All after June 2027 and before January 2030 _should we simulate the yield or the population??_

Case A) 9 triggers with 30 sq deg area -> yealds 3 visible. Only cover if center in footprint

Case B) 15 triggers with 50 sq deg area -> yields 5 observable. Only cover if center in footprint

Case C) 24 triggers with 100 sq deg area -> yelds 8 observable. Only cover if center in footprint

Case D) 6 triggers with 250 sq deg area -> yelds 4 observable. Only cover if center in footprint

Case E) 12 triggers with 500 sq deg area -> yelds 2 observable. Only cover if center in footprint

Simulate Tw observability window of area of interest randomly

**A:**
*  Night 0:
      
        if Tw> 4.7h:

           4x(u)griz(y)) x 3 passes
   
        if 3.2h < Tw < 4.7h:
   
           4x(u)griz(y)) x 2 passes
   
        if Tw < 3.2h:
   
           4x(u)griz(y) x 1 passes
* Night 1,2:

        4x(u)griz(y) x 1 pass
* Night 3 only 1/3 of the times:

        4x(u)griz(y) x 1 pass
  

**B/C:**
* Night 0:

      if Tw> 4.7h:

        (4xgri) x 3 passes

      else if 3.2h < Tw < 4.7h:

        (4xgri) x 3 passes

      else if Tw < 3.2h:

        (4xgri) + (4xr)
* Night 1,2:

      6xri x 1 pass

* Night 3 only 1/3 of the times:

      6xri x 1 pass
      

**D/E:**
* Night 0:

        1xgr x 1 pass
  
* Night 1,2:
  
        4xgr 30s x 1 pass
 
* Night 3 only 1/3 of the times:
  
        4xgr x 1 pass

**Simulate lensed BNS ToOs**

_N total = 2_

All after June 2027 and before January 2030

Case A) 1 trigger 900 sq deg in observable sky

Case B) 1 trigger 15 sq deg in observable sky

**A:** 

* Night >= 0:
  
    1xg + 3xr x 1 pass

**B:**
  * Night >= 0:
  
    50xg + 50xr x 1 pass

**Simulate BBH ToOs**

_N total = 7_

All after June 2027 and before January 2030

Case A) 9 triggers with 10 sq deg area -> yields 3 observable

Case B) 6 triggers with 20 sq deg area -> yields 2 observable

Case C) 6 triggers with 30 sq deg area -> yields 2 observable


**A/B/C**

if bright: filters = (u)g(r)i 

else: filters = riz

* Night 0,2,7,9,39:
  
      1xfilters x 1 pass
