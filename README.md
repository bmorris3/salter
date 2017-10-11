# salter
Stellar Active Latitudes with Transiting Exoplanet Residuals

### The Problem
Do other stars have active latitudes, as we seen on [the Sun](https://en.wikipedia.org/wiki/Spörer%27s_law)? If so, it's something very basic and important related to the dynamo and convection properties.

One of the best (only?) other determinations of this is from analysis of *Kepler* data for Hat-P-11 by [Morris et al.](https://arxiv.org/abs/1708.02583) This system is very unique, in both geometry and brightness.

How else can we constrain the starspots/activity as a function of latitude?!

### Idea
Let's use ensemble transit data. If we use transits with many different impact parameters, we'll be sensitive to a wide range of latitudes. Now the game is measuring the in- vs. out-of-transit scatter as a function of latitude due to starspot occultations.

### Prep your local data cache
1. First open and run the notebook `choose_targets.ipynb`. 

2. In the repo, run this to construct and cache the database of kepler light curves (runtime = ~1 hr on conference wifi)
```
python -c "from salter.cache import cache_light_curves; cache_light_curves()"
```
3. Run this to cache a local copy of the joined tables from the NASA Exoplanet Archive and the Exoplanet Orbit Database:
```
python -c "from salter.cache import cache_planets_table; cache_planets_table()"
```

### Example notebook: 

Open up the notebook `show_lc.ipynb` to see what the package can do.

### Dependencies

* numpy 
* matplotlib
* scipy
* h5py, hdf5
* astropy
* batman
* kplr
* astroquery
