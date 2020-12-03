# Spatial Map Generator
This repo contains the code used to produce approximated maps of hit-angle distributions of photon emission in liquid argon.

# Requirements
To run this macros, you need to install `ROOT` on your pc.

The code has been tested with `ROOT 6.20` on a machine running `Ubuntu 18.04.5 LTS`.

# How to run
To generate the Spatial Maps, you have to use the macro `createMap.cpp`.

## Default configuration
From the current directory, run the following command:

```
$>root
root [0] .L createMap.cpp
root [1] createMap()
```

or directly:

```
$>root createMap.cpp
```

## Basic customization
The easiest parameters to configure are the following:
- `nOpticsPerPoint`: number of optical photons emitted for each point
- `minR`: starting radius (in mm) for sampling
- `maxR`: ending radius (in mm) for sampling
- `angleBins`: number of angle bins for the produced maps

From the current directory, run the following command:
```
$>root
root [0] .L createMap.cpp
root [1]
root [2] Int_t nOpticsPerPoint = 5000;     // choose the number of photons sampled in each point
root [3] Double_t minR = 300.0;            // the simulation will start sampling photons from R=300
root [4] Double_t maxR = 500.0;            // until reaching R=500
root [5] Int_t angleBins = 100;            // the produced map have 100 bins for the angle-axis
root [6]
root [7] createMap(nOpticsPerPoint, minR, maxR, angleBins)
```

## Advanced customization: 
You can also customize further parameters by using two structs: `SimConfig` for the simulation parameters, `LGND200Geom` for the simulated geometry.

### Simulation parameters
Using a global shared variable named `conf` you can change the following simulation parameters.
|        Parameter           | Description |
| -------------------------- | ----------- |
| `attenuationLength`        | attenuation length of photons in the material (*default l=500.0*)  |
| `shroudCaptureProb`        | geometric acceptance of fibers volume compared with the ideal cylinder (*default pr=0.54*) |
| `shiftAngleForXYSampling`  | angular slice (*in rad*) from which produce photons by symmetry assumptions (*default ang=0.22 rad*) |

For example, you could run:
```
$>root
$>.L createMap.cpp
$>
$>conf.attenuationLength = 100.0;          // choose attenuation length
$>conf.shroudCaptureProb = 0.75;           // choose the prob of acceptance by the shrouds
$>conf.shiftAngleForXYSampling = 0.22;     // choose the angle shifting from which sample photons
$>
$>Int_t nOpticsPerPoint=10000;
$>createMap(nOpticsPerPoint)
```

### Geometry
Using the global shared variable `geom`, you can change the parameters of the simulated geometry.

The following parameters concern the definition of the fiber shrouds:
|        Parameter           | Description |
| -------------------------- | ----------- |
| `innerShroudRadius`        | distance of the inner shroud from the origin (*default r=175 mm*) |
| `outerShroudRadius`        | distance of the outer shroud from the origin (*default r=295 mm*) |
| `topZShroud`  | highest Z coordinate for both the shroud (*default z=+845 mm*) |
| `bottomZShroud`  | lowest Z coordinate for both the shroud (*default z=-845 mm*)|

The following parameters concern the definition of the Germanium detector:
|        Parameter           | Description |
| -------------------------- | ----------- |
| `radiusGeStrings`  | distance of the Germanium strings from the origin (*default r=235 mm*) |
| `topZGeStrings`    | highest Z coordinate of the Germanium strings (*default r=425 mm*) |
| `bottomZGeStrings` | lowest Z coordinate of the Germanium strings (*default z=-425 mm*) |
| `radGeCrystal`     | radius of each single Germanium crystal (*default r=40 mm*)|
| `nGeStrings`       | number of Germanium strings which compose the detector (*default n=14*)|