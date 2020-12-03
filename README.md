# Spatial Map Generator
This repo contains the code used to produce approximated maps of hit-angle distributions of photon emission in liquid argon.

# Requirements
To run this macros, you need to install ROOT on your pc.

The code has been tested with ROOT 6.20 on a machine running Ubuntu 18.04.5 LTS.

# How to run
To generate the Spatial Maps, you have to use the macro `createMap.cpp`.

## Default configuration
From the current directory, run the following command:

```
$>root
$>.L createMap.cpp
$>createMap()
```

or directly:

```
$>root createMap.cpp
```

## Basic customization
The first customizable parameters are the following:
- `nOpticsPerPoint`: number of optical photons emitted for each point
- `minR`: starting radius (in mm) for sampling
- `maxR`: ending radius (in mm) for sampling
- `angleBins`: number of angle bins for the produced maps

From the current directory, run the following command:
```
$>root
$>.L createMap.cpp
$>
$>Int_t nOpticsPerPoint = 5000;     // choose the number of photons sampled in each point
$>Double_t minR = 300.0;            // the simulation will start sampling photons from R=300
$>Double_t maxR = 500.0;            // until reaching R=500
$>Int_t angleBins = 100;            // the produced map have 100 bins for the angle-axis
$>
$>createMap(Int_t nOpticsPerPoint = 1000, Double_t minR = 0, Double_t maxR = 710, Int_t angleBins = 100)
```

## Advanced customization: 
You can also customize further parameters by using two structs: `SimConfig` for the simulation parameters, `LGND200Geom` for the simulated geometry.

### Simulation parameters
The following simulation parameters can be customized:
- `SimConfig.attenuationLength`: the attenuation length of photons in the material, default l=500.0 in liquid argon
- `SimConfig.shroudCaptureProb`: the geometric acceptance of fibers volume compared with the ideal cylinder modeled as probability, default pr=0.54 in LGND200
- `SimConfig.shiftAngleForXYSampling`: sampling angular slice from which produce photons, by symmetry assumptions it is 0.22 in LGND200

### Geometry
The following parameters of the geometry can be customized:
The fiber shrouds:
- `LGND200Geom.innerShroudRadius`: distance of the inner shroud from the origin, by default r=175mm in LGND200
- `LGND200Geom.outerShroudRadius`: distance of the outer shroud from the origin, by default r=295mm in LGND200
- `LGND200Geom.topZShroud`: highest Z coordinate for both the shroud, by default z=+845mm in LGND200
- `LGND200Geom.bottomZShroud`: lowest Z coordinate for both the shroud, by default z=-845mm in LGND200
The Germanium detector:
- `LGND200Geom.radiusGeStrings`: distance of the Germanium strings from the origin, by default r=235mm in LGND200
- `LGND200Geom.topZGeStrings`: highest Z coordinate of the Germanium strings, by default z=+425mm in LGND200
- `LGND200Geom.bottomZGeStrings`: lowest Z coordinate of the Germanium strings, by default z=-425mm in LGND200
- `LGND200Geom.radGeCrystal`: radius of each single Germanium crystal, by default r=40mm in LGND200
- `LGND200Geom.nGeStrings`: number of Germanium strings which compose the detector, by default n=14 in LGND200