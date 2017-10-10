# salter
Stellar Active Latitudes with Transiting Exoplanet Residuals

## The Problem
Do other stars have active latitudes, as we seen on [the Sun](https://en.wikipedia.org/wiki/Spörer%27s_law)? If so, it's something very basic and important related to the dynamo and convection properties.

One of the best (only?) other determinations of this is from analysis of *Kepler* data for Hat-P-11 by [Morris et al.](https://arxiv.org/abs/1708.02583) This system is very unique, in both geometry and brightness.

How else can we constrain the starspots/activity as a function of latitude?!

## Idea
Let's use ensemble transit data. If we use transits with many different impact parameters, we'll be sensitive to a wide range of latitudes. Now the game is measuring the in- vs. out-of-transit scatter as a function of latitude due to starspot occultations.

## Prep your local data cache
In the repo, run this to construct and cache the database of kepler light curves (runtime = ~1 hr on conference wifi)
```
python -c "from salter.cache import cache_light_curves; cache_light_curves()"
```
Run this to cache a local copy of the joined tables from the NASA Exoplanet Archive and the Exoplanet Orbit Database:
```
python -c "from salter.cache import cache_joined_table; cache_joined_table()"
```

## Dependencies

* numpy 
* matplotlib
* h5py, hdf5
* astropy
* batman
* kplr
