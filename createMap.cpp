#include <iostream>
#include <TEntryList.h>
#include <TParameter.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TBranch.h>
#include <TChain.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TMath.h>
#include <TRandom.h>
#include <array>
#include <vector>
#include <set>
#include <assert.h>
#include "src/GeometryLib.cpp"

#define INNER_REGION -1
#define MIDDLE_REGION 0
#define OUTER_REGION +1

#define INNER_SHROUD 0
#define OUTER_SHROUD 1
#define GE_VOLUME -2

using namespace std;

struct LGND200Geom{
    // Fiber shrouds
    const double innerShroudRadius = +175.0;
    const double outerShroudRadius = +295.0;
    const double topZShroud        = +845.0;
    const double bottomZShroud     = -845.0;
    // Ge Detector
    const double topZGeStrings     = +425.0;
    const double bottomZGeStrings  = -425.0;
    const double radiusGeStrings   = +235.0;
    const double radGeCrystal      =  +40.0;
    const int nGeStrings = 14;
    vector<double> geCenters_x;
    vector<double> geCenters_y;
};

struct SimConfig{
    double attenuationLength       = 500.0;
    double shroudCaptureProb       =  0.54;     // geometric coverage of fiber shrouds (fibers / ideal cylinder)
    double shiftAngleForXYSampling =  0.22;	    // symmetry assumptions, e.g. 0.22 = 2 * PI / nGeStrings / 2, for nGeStrings=14
};

// Rnd Generator
int seed = 123456789;
TRandom * rnd = new TRandom(seed);
LGND200Geom geom;
SimConfig conf;

void initializePositionGeCrystals(){
    // place the Ge strings at the same distance, over the radius defined in the geometry
    cout << "[Info] Initialization Ge Crystals...\n";
    Double_t geAngle = 0;
    for(Int_t ige=0; ige<geom.nGeStrings; ige++){
        // compute center of Ge crystal
        Double_t x0 = geom.radiusGeStrings * cos(geAngle);
        Double_t y0 = geom.radiusGeStrings * sin(geAngle);
        // append coordinates
        geom.geCenters_x.push_back(x0);
        geom.geCenters_y.push_back(y0);
        // increment angle shift for next string
        geAngle += 2 * TMath::Pi() / geom.nGeStrings;
        // debug
        cout << "\tCrystal " << ige+1 << ":\t";
        cout << fixed << setprecision(2) << "x: " << x0 << ",\ty: " << y0 << ",\tr: " << sqrt(x0*x0+y0*y0) << endl;
    }
}

struct Point getHitOnGermanium(Point point, Point dirPoint, Bool_t enableGe){
    // Idea: the method to find an intersection  assumes the circle to be centered in 0,0
    // We can translate the points to have each Ge crystal centered
    Double_t distance = INF;
    struct Point closestGeHit = Point();
    if(!enableGe)
        return closestGeHit;	// If Ge disable, return always NoHit!
    // We aim to average the result over a slice, then sample the rotation angle
    Double_t original_x1 = point.x, original_y1 = point.y;
    Double_t original_x2 = dirPoint.x, original_y2 = dirPoint.y;
    for(int i=0; i<geom.geCenters_x.size(); i++){
        // Shift the point coordinates w.r.t. the center of crystal
        point.x = point.x - geom.geCenters_x[i];
        point.y = point.y - geom.geCenters_y[i];
        dirPoint.x = dirPoint.x - geom.geCenters_x[i];
        dirPoint.y = dirPoint.y - geom.geCenters_y[i];
       
        vector<Point> geHits = computeIntersection(point, dirPoint, geom.radGeCrystal);
        if(geHits.size() > 0){
            Point geHit = geHits[0];    // closest hit on Ge crystal
            geHit.x += geom.geCenters_x[i];
            geHit.y += geom.geCenters_y[i];
            // Note: no shift on Z
            if(geHit.distance < distance){
                distance = geHit.distance;
                closestGeHit = geHit;
            }
        }
        // Restore original coord
        point.x = original_x1;
        point.y = original_y1;
        dirPoint.x = original_x2;
        dirPoint.y = original_y2;
    }
    return closestGeHit;
}


Bool_t checkIfEnableGe(Point prodPoint, Point dirPoint){
    // Compute intersection with the Ge volume, i.e. the cylinder that includes the Ge Strings
    double geVolumeRadius = geom.radiusGeStrings + geom.radGeCrystal;
    vector<Point> geHits = computeIntersection(prodPoint, dirPoint, geVolumeRadius);
    if(geHits.size()==0)
        return false;	// No interesection at all, no difference between the two maps
    // If there is an interesection, look at where it occurs in Z
    // Use parametrization of line: line = P + t D, where P is a vector (point), D direction vector
    Point geHit = geHits[0];    // consider only the first it on Ge crystal
    Double_t t = (geHit.x - prodPoint.x) / (dirPoint.x-prodPoint.x);	// From param: t = (x-x1)/dx
    Double_t z_int = prodPoint.z + t * (dirPoint.z-prodPoint.z);   		// z = z1 + t * dz, where z1 point, dz direction vector
    if(geom.bottomZGeStrings < z_int && z_int < geom.topZGeStrings)
        return true;	// Interesection with Z in [BOTTOMGE,TOPGE], enable crystals
    return false;		// Intersection with Z outside the Ge volume, disable crystals
}

struct Point generateRandomPointInAngleSlice(Double_t radius, Double_t minZ=geom.bottomZShroud, Double_t maxZ=geom.topZShroud){
    Double_t angleShift = conf.shiftAngleForXYSampling;
    Double_t phi = -angleShift + rnd->Rndm() * (2*angleShift);  // Random Phi in angle
    struct Point point = Point(radius * cos(phi),
                               radius * sin(phi),
                               minZ + rnd->Rndm() * maxZ); // Random Z position
    return point;
}

pair<struct Point, Double_t> generatePhotonEmission(struct Point point){
    Double_t theta = rnd->Rndm() * TMath::Pi();
    Double_t phi = rnd->Rndm() * 2 * TMath::Pi();
    // generate direction point
    Point dirPoint = Point(point.x + 1 * sin(theta) * cos(phi),
                           point.y + 1 * sin(theta) * sin(phi),
                           point.z + 1 * cos(theta));
    Double_t lenTrajectory = 0;	// distance to hit the plane delimited by upper or lower boundary
    if(theta <= TMath::Pi() / 2){
        lenTrajectory = (geom.topZShroud-point.z)/cos(theta);	    // trajectory towards top Z
    }else{
        lenTrajectory = (geom.bottomZShroud-point.z)/cos(theta);    // trajectory towards bottom Z
    }
    return make_pair(dirPoint, lenTrajectory);
}

vector<struct Point> computeAllHits(struct Point prodPoint, struct Point dirPoint, Double_t lenTrajectory, Bool_t enableGe){
    // Possible scenarios: outer hit (close), inner hit (close), inner hit (far), outer hit(far).
    // assert: if hit shroud, always the outer first
    // assert: if 2nd hit, check if no Ge before
    vector<struct Point> hits;
    // 2 cases:
    struct Point geHit;
    geHit = getHitOnGermanium(prodPoint, dirPoint, enableGe);
    // First, check for germanium hit
    Double_t distanceToGeHit = INF;
    if(geHit.isDefined()){
        distanceToGeHit = geHit.distance;
    }
    for(auto shroud : {INNER_SHROUD, OUTER_SHROUD}) {
        Int_t shroudRadius;
        if(shroud==INNER_SHROUD)
            shroudRadius = geom.innerShroudRadius;
        else
            shroudRadius = geom.outerShroudRadius;
        vector<Point> shroudHits = computeIntersection(prodPoint, dirPoint, shroudRadius, lenTrajectory);
        for(auto &hit : shroudHits){
            if(hit.isDefined() && prodPoint.checkSameDirection(dirPoint.x, dirPoint.y, hit) && hit.distance<distanceToGeHit) {    // Discard no hit on outer shroud
                hit.shroud = shroud;
                hits.push_back(hit);
            }
        }
    }
    std::sort(hits.begin(), hits.end());    // sort ascending distance
    return hits;
}

void runToyOpticsFromPoint(Double_t radius, Int_t nOptics, TH1D * &prInnerD, TH2D* &innerMap, TH2D* &outerMap){
    cout << "\rRunning at R=" << radius << std::flush;
    // Flag sample Ge
    Bool_t enableGe;
    // Loop on angles
    Double_t kInnerHits = 0, kOuterHits = 0;
    // Get coordinate
    struct Point prodPoint = Point();
    Int_t kGe = 0;
    while(nOptics>0){
        // Get random point (y1 according to angle shifting, z1 random) 
        prodPoint = generateRandomPointInAngleSlice(radius);
        pair<struct Point, Double_t> photonEmission = generatePhotonEmission(prodPoint);
        struct Point dirPoint = photonEmission.first;
        Double_t lenTrajectory = photonEmission.second;
        enableGe = checkIfEnableGe(prodPoint, dirPoint);
        assert(lenTrajectory>=0);        
        // Compute the possible hits
        Int_t pointPlacement = getPointPlacement(prodPoint.x, prodPoint.y, geom.innerShroudRadius, geom.outerShroudRadius);
        vector<Point> hits = computeAllHits(prodPoint, dirPoint, lenTrajectory, enableGe);        
        assert(pointPlacement!=INNER_REGION || hits.size()<=2);     //INNER_REGION => nr hits <=2
        assert(pointPlacement!=MIDDLE_REGION || hits.size()<=3);    //MIDDLE_REGION => nr hits <=3
        assert(pointPlacement!=OUTER_REGION || hits.size()<=4);     //OUTER_REGION => nr hits <=4

        if(hits.size() > 0) {
            vector <Double_t> probabilities, norm_probabilities;
            Int_t i = 0;
            Double_t NORM_FACTOR = 0;
            for (auto hit : hits) {
                // prob: (1-shroudCaptPr)^(i-1) * shroudCaptPr * Exp(-dist/L)
                i++;    // avoid i=0
                Double_t pr = pow(1 - conf.shroudCaptureProb, i - 1) * conf.shroudCaptureProb * TMath::Exp(-hit.distance / conf.attenuationLength);
                probabilities.push_back(pr);
                NORM_FACTOR += pr;
            }
            for (auto pr : probabilities) {
                norm_probabilities.push_back(pr / NORM_FACTOR);
            }
            // Sampling the hit
            Double_t rndChoice = rnd->Rndm();
            Int_t id_currentHit = 0;
            while (rndChoice > norm_probabilities[id_currentHit]) {
                rndChoice -= norm_probabilities[id_currentHit];
                id_currentHit++;
            }
            assert(id_currentHit < hits.size());
            Point sampledHit = hits[id_currentHit];
            // Fill map with the sampled hit
            assert(sampledHit.isDefined());
            assert(sampledHit.shroud == INNER_SHROUD || sampledHit.shroud == OUTER_SHROUD);
            double angle = atan2(sampledHit.y, sampledHit.x);
            if (sampledHit.shroud == INNER_SHROUD) {
                innerMap->Fill(radius, angle);
                kInnerHits += 1;
            } else if (sampledHit.shroud == OUTER_SHROUD) {
                outerMap->Fill(radius, angle);
                kOuterHits += 1;
            }           
            // Update counters and debug    
            nOptics--;
        if (enableGe)
                kGe++;
        }
    }
    if(kInnerHits + kOuterHits > 0){
        prInnerD->Fill(radius, kInnerHits/(kInnerHits+kOuterHits));
    }
}

void createMap(Int_t nOpticsPerPoint=1000, Double_t minR=0, Double_t maxR = 710, Int_t angleBins=100){
    // initialize Ge Crystals
    initializePositionGeCrystals();
    // definition of bins
    Double_t rbins = maxR - minR + 1;
    Int_t minAngle = -ceil(TMath::Pi()), maxAngle = +ceil(TMath::Pi());
    // create output file and maps
    TString outfile;
        outfile.Form("SpatialMap_%dR_%dAngleSlices_%dops_AttLen%f_NewSampling_Last.root", (int)rbins, angleBins, nOpticsPerPoint, conf.attenuationLength);
        TFile * file = new TFile(outfile, "RECREATE");
    TH1D * prInnerD = new TH1D("PrInnerDet", "Pr ~ Fract. Inner/(Inner+Outer) Detections", rbins, minR, maxR);
    TH2D * innerMap = new TH2D("InnerMap", "R-Angle Inner Shroud Map", rbins, minR, maxR, angleBins, minAngle, maxAngle);
    TH2D * outerMap = new TH2D("OuterMap", "R-Angle Outer Shroud Map", rbins, minR, maxR, angleBins, minAngle, maxAngle);
    // iterate on radius and sample optical photons
    cout << "[Info] Sampling photon hits...\n";
    for(Double_t x=minR; x<=maxR; x += 1){
        runToyOpticsFromPoint(x, nOpticsPerPoint, prInnerD, innerMap, outerMap);
    }
    // write maps to output file
    prInnerD->Write();
    innerMap->Write();
    outerMap->Write();
    file->Close();
    cout << "\n[Info] Done.\n";
}
