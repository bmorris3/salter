Documentation
=============

This is the documentation for ``salter``.

Source Code
-----------
`...is available on GitHub <https://github.com/bmorris3/salter/>`_.

Installation
------------

Install via pip (via GitHub)::

    pip install https://github.com/bmorris3/salter/archive/master.zip

Dependencies
------------
- numpy
- scipy
- matplotlib
- astropy
- astroquery
- batman-package
- kplr
- h5py

Prep your local data caches
---------------------------

1. First open and run the notebook `choose_targets.ipynb`.

2. In the repo, run this to construct and cache the database of kepler light curves (runtime = ~1 hr on conference wifi)::

    python -c "from salter import cache_light_curves; cache_light_curves()"

3. Run this to cache a local copy of the joined tables from the NASA Exoplanet Archive and the Exoplanet Orbit Database::

    python -c "from salter import cache_planets_table; cache_planets_table()"

.. automodapi:: salter
