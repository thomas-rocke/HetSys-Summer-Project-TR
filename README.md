# Summer Project Repo - Tom Rocke

This repository hosts the summer project reproducable result, calculating the formation energy of a P interstitial point defect in an InP bulk.

The repo assumes an existing and up-to-date install of the quippy Python module - installation instructions can be found in "Installing_QUIP.md"


"Reproducible_Result.ipynb" is the main component under peer review. "GAP_Interfaces.py" may also be of interest to current/future GAP users, as it contains both a wrapper for how we call GAP models ("GapQUIP" function), but also a Python class implementing a portion of GAP ("GapPy") which may be useful to understanding the internals of GAP prediction.