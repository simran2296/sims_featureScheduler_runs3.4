An initial port of v3.3 up to v3.4. Refactoring to add documentation and make it easier to drop into
`rubin_scheduler.scheduler.example_scheduler`

changes since v3.3

* Updated to `IntRounded` to ensure code repeats cross-platform
* Patched small bug that lead to rare longer than actual reported slewtime
* update to N good seeing basis function that makes a minor change in execution
* Swap the zenith mask with the more general shadow mask basis function
* update to scripted surveys and DDF preschuling so DDF sequences don't start right before sunrise twilight
* Update to traveling salesman solver to ensure cross-platform repeatability.


notes on MAF runtime
real    46m56.927s
real    135m42.152s
real    313m58.593s
real    869m7.487s

So the DDFs and Science would be good to break up and run more parallel.
