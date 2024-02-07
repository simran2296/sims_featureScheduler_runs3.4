
too.py has initial update of old ToO simulations. Checks different ToO rates and two follow up strategies.

~~Assumes all ToO events have a 6.5 degree radius search area.~~

Only attempt to follow ToO events that are within our regular footprint. Give up if the event is over 3 days old.


**TBD:** 
**- check footprint size breakdown**
**- check p(t) function**
**- check fraction that should go to completion** 
------

For ToO workshop in March 2024, we would like to simulate:

**Nevents/year**: 
- 8
- 16
- 32
- 64 


**Footprint (same for all sims)**: 
- 50% 100 sqdeg
- 35% 50 sqdeg
- 15% 20 sqdeg


**Strategy**: 

SET 1:
A) 
   NIGHT 0: _g+z_ 1 image per filter over full footprint t0+1 hour after trigger or at night start 
   REPEATS:   t0+2 t0+4 hours from trigger
   NIGHT 1: 6x_g_ + 6x_z_ 
B) 
   NIGHT 0: _g+i+z_ 1 image per filter over full footprint t0+1 hour after trigger or at night start 
   REPEATS:   t0+2 t0+4 hours from trigger
   NIGHT 1: 6x_g_ + 6x_z_ 
C) 
   NIGHT 0: _g+r+i+z_ 1 image per filter over full footprint t0+1 hour after trigger or at night start 
   REPEATS:   t0+2 t0+4 hours from trigger
   NIGHT 1: 6x_g_ + 6x_z_ 
D) 
   NIGHT 0: _(u)+g+r+i+z+(y)_ 1 image per filter over full footprint t0+1 hour after trigger or at night start 
   REPEATS:   t0+2 t0+4 hours from trigger
   NIGHT 1: 6x_g_ + 6x_z_ 

SET 2 A, B, C, D, same as set 1 but with NIGHT 2:  6x_g_ + 6x_z_


**ALL SETS INTERRUPTS**
20% of the ToOs should be taken to completion
80% of the ToOs shoudl be interrupted at time t0+t with t drawn from p(t)~t**2 between t0 and t0+100 hours









