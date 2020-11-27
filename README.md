# Cheops Observability Interpolator

Until ESA makes a useable Cheops observability web tool, here is what I use to figure out if a target is within Cheops' cone of shame.
This uses the interpolated maps of [ESA Cheops Sky](https://www.cosmos.esa.int/web/cheops/the-cheops-sky) to estimate the time a target is observable, what the max efficiency is, and what the combined observability time is (all assuming a minimum threshold of 50% efficiency).

Check out [Example.ipynb](Example.ipynb) for a quick example.

I caution anyone - *do not use this tool to decide exactly which targets to observe and when*. For final target selection and scheduling, the official virtual machine-based feasibility checker should be used. But to figure out which of your longlist of candidates are likely to be observable with Cheops, this tools should work fine.
