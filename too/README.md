
too.py has initial update of old ToO simulations. Checks different ToO rates and two follow up strategies. Rates of 5, 10, 100 events per year. 
followups of:
1) gz, gz 1 hour later
4) gz, gz 1 hour later, gz 2 hours later, gz 4 hours later.
(note z is always available, one of u and y will not be available)

Also run the ignore variant that doesn't try to count the ToO observations as part of
the regular survey.

~~Assumes all ToO events have a 6.5 degree radius search area (132 sq deg).~~

Only attempt to follow ToO events that are within our regular footprint. Give up if the event is over 3 days old.







