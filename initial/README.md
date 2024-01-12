An initial port of v3.3 up to v3.4. Refactoring to add documentation and make it easier to drop into
`rubin_scheduler.scheduler.example_scheduler`

* Updated to `IntRounded` to ensure code repeats cross-platform
* Patched small bug that lead to rare longer than needed slewtime.
* update to N good seeing basis function that makes a minor change

To do:
* Swap the zenith mask with the more general shadow mask basis function.
