# Spatial Map Generator
This repo contains the code used to produce approximated maps of hit-angle distributions of photon emission in liquid argon.

![alt text](doc/inner_map.png?raw=true      "Inner Shroud Map")
![alt text](doc/outer_map.png?raw=true      "Outer Shroud Map")
![alt text](doc/pr_inner_det.png?raw=true   "Prob. Inner Detection Map")

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

For example, suppose you want to focus on a reduced portion of liquid argon, with radius between `300` and `500`.
Also, you would like to have a better resolution of the map, increasing the photons per point to `5000` and the angle bins to `250`.
Then, you can run the following commands:
```
$>root
root [0] .L createMap.cpp
root [1]
root [2] Int_t nOpticsPerPoint = 5000;    // increase resolution by increasing nr photons per point
root [3] Double_t minR = 300.0;           // focus on smaller region, starting from radius 300
root [4] Double_t maxR = 500.0;           // and ending at radius 500
root [5] Int_t angleBins = 250;           // increase resolution of angles, setting 250 bins
root [6]
root [7] createMap(nOpticsPerPoint, minR, maxR, angleBins)
```

## Advanced customization: 
You can also customize further parameters by using two structs: `SimConfig` for the simulation parameters, `LGND200Geom` for the simulated geometry.

### Simulation parameters
Using a global shared variable struct `conf` you can change the following simulation parameters.
|        Parameter           | Description |
| -------------------------- | ----------- |
| `attenuationLength`        | attenuation length of photons in the material (*default l=500.0*)  |
| `shroudCaptureProb`        | geometric acceptance of fibers volume compared with the ideal cylinder (*default pr=0.54*) |
| `shiftAngleForXYSampling`  | angular slice (*in rad*) from which produce photons by symmetry assumptions (*default ang=0.22 rad*) |

For example, suppose you want to generate a spatial map with an attenuation length `l=100`.
Also, you want to have a better approximation by increasing the number of photons per point to `10000`.
Then, you can run the following commands:
```
$>root
root [0] .L createMap.cpp
root [1] 
root [2] conf.attenuationLength = 100.0;    // set att length to 100
root [3] Int_t nOpticsPerPoint=10000;       // set nr photons per point to 10K
root [4] 
root [5] createMap(nOpticsPerPoint)
```

### Geometry
Using the global shared struct `geom`, you can change the parameters of the simulated geometry.

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

For example, assume you want to create a spatial map for a novel detector geometry, composed by `20` strings and smaller crystals with `r=10 mm`. Then, you can run the following commands:
```
$>root
root [0] .L createMap.cpp
root [1] 
root [2] geom.nGeStrings = 20;      // set the number of Ge strings to 20
root [3] geom.radGeCrystal = 10;    // set the cristal radius to 10 mm
root [4] 
root [5] createMap()
```