//TODO: license??

#include "cameras/perspective.h"
#include "cameras/distortion.h"
#include "paramset.h"
#include <iostream> // TODO:remove
#include <vector>
#include <cmath>
#include <unordered_map>
#include "sampling.h"
#include <cassert>

#define NO_MODEL "MODEL_NOT_FOUND"
#define POLY_DEGREE 5

namespace pbrt {

  std::unordered_map<std::string, int> DistortionCamera::numCoeffsForModel = { {"poly3lensfun", 1},
                                                                               {"poly5lensfun", 2},
                                                                               {"ptlens", 3} };

  DistortionCamera::DistortionCamera(const AnimatedTransform &CameraToWorld,
                                     const Bounds2f &screenWindow, Float shutterOpen,
                                     Float shutterClose, Float lensRadius, Float focalDistance,
                                     Float fov, Film *film, const Medium *medium,
                                     std::string distortionModel, DistortionCamera::coeffVec coeffs,
                                     int centerOffsetX, int centerOffsetY)
    : ProjectiveCamera(CameraToWorld, Perspective(fov, 1e-2f, 1000.f),
                       screenWindow, shutterOpen, shutterClose, lensRadius,
                       focalDistance, film, medium),
      distortionModel(distortionModel), coeffs(coeffs),
      centerOffsetX(centerOffsetX), centerOffsetY(centerOffsetY), fittedCoeffs(){

        /* ------------ Validate given model ---------------------------*/
        if (distortionModel == NO_MODEL) {
          Error("No model for lense distortion given. Rendering will be done without distortion, equivalent to PerspectiveCamera.");
          fittedCoeffs = coeffVec({0, 1});
        }
        else {
          auto modelIter = numCoeffsForModel.find(distortionModel);
          if (modelIter == numCoeffsForModel.end()) {
            Error("Model %s is unsupported (did you make a typo?). Rendering will be done without distortion, equivalent to PerspectiveCamera.",
                  distortionModel.c_str());
            fittedCoeffs = coeffVec({0, 1});
          }
          else {
            int requiredNum = modelIter->second;
            int givenNum = coeffs.size();
            if (givenNum > requiredNum) {
              Error("Model %s requires %d coefficients, but %d provided. I will use only the first %d.", 
                    distortionModel.c_str(), requiredNum, givenNum, requiredNum);
              auto first = coeffs.begin();
              auto last = coeffs.begin() + requiredNum;
              coeffs = coeffVec(first, last);
            }
            else if (givenNum < requiredNum) {
              Error("Model %s requires %d coefficients, but only %d provided. I will assume the missing %d to be zero.",
                    distortionModel.c_str(), requiredNum, givenNum, requiredNum - givenNum);
              coeffVec zeros = coeffVec(requiredNum - givenNum, .0);
              coeffs.insert(coeffs.end(), zeros.begin(), zeros.end());
            }
            assert(coeffs.size() == requiredNum);
            fittedCoeffs = InvertDistortion(coeffs, POLY_DEGREE);
          }
        }

        /* -------- Define transformations for ray generation --------*/
        Float xRes = film->fullResolution.x;
        Float yRes = film->fullResolution.y;
        Float cornerRadius = sqrt(pow(xRes, 2) + pow(yRes, 2)) / 2.;
        std::cout << "Corner radius: " << cornerRadius << std::endl;
        NormalizeToCornerRadius = Scale(1. / cornerRadius, 1. / cornerRadius, 1.);
        Denormalize = Inverse(NormalizeToCornerRadius);

        Point3f p = Point3f(xRes, yRes, 0);
        std::cout << "Ecke normalisiert" << NormalizeToCornerRadius(p) << std::endl;

        /* ---- transform image center ----*/
        std::cout << "Center offset: " << centerOffsetX << ", " << centerOffsetY << std::endl;
        imageCenterNormalized = NormalizeToCornerRadius(Point3f(xRes / 2., yRes / 2., 0));
        imageCenterNormalized = imageCenterNormalized + Vector3f(centerOffsetX, centerOffsetY, 0);
        std::cout << "Image Center normalised: " << imageCenterNormalized << std::endl;

        Vector3f rad = Vector3f(imageCenterNormalized) - Vector3f(NormalizeToCornerRadius(p));
        std::cout << "Normalised distance from image center to Ecke: " << rad.Length() << std::endl;
    }

  DistortionCamera::coeffVec DistortionCamera::InvertDistortion(DistortionCamera::coeffVec coeffs,
                                                                int polyDegree) {
    // create x vector for sampling distortion values
    // and y vector with sampled values
    int sampleSize = 1000;     // what is a reasonable value here?
    Float scale(sampleSize);
    std::vector<Float> x(sampleSize);
    std::vector<Float> y(sampleSize);

    // get pointer to model function
    Float (*modelFunc)(const Float, const coeffVec);
    if (distortionModel == "poly3lensfun")
      modelFunc = ModelPoly3LensFun;
    else if (distortionModel == "poly5lensfun")
      modelFunc = ModelPoly5LensFun;
    else if (distortionModel == "ptlens")
      modelFunc = ModelPTLens;
    else
      Error("Model %s not supported. this should have been caught in the constructor!", distortionModel.c_str());

    // fill sample vectors according to the model used
    for (int i = 0; i < sampleSize; i++) {
      x[i] = i / scale;
      y[i] = (*modelFunc)(x[i], coeffs);
    }

    coeffVec polyCoeffs = fitPolyCoeffs(y, x, polyDegree);
    return polyCoeffs;
  }

  // the fitted radial model is applied here
  Point3f DistortionCamera::CalculateRayStartpoint(const CameraSample &sample) const {
    Float centerX = imageCenterNormalized.x;
    Float centerY = imageCenterNormalized.y;
    Point3f pNormalized = NormalizeToCornerRadius(Point3f(sample.pFilm.x, sample.pFilm.y, 0));
    Float radius = sqrt(pow(pNormalized.x - centerX, 2) + pow(pNormalized.y - centerY, 2));
    Float rNew = evalPolynomial(fittedCoeffs, std::vector<Float>({radius}))[0];
    Float rRatio = rNew / radius;
    return Denormalize(Point3f(pNormalized.x * rRatio + centerX * (1 - rRatio),
                               pNormalized.y * rRatio + centerY * (1 - rRatio), 0));
  }


  Float DistortionCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
    ProfilePhase prof(Prof::GenerateCameraRay);
    Point3f pCamera = RasterToCamera(CalculateRayStartpoint(sample));
    *ray = Ray(Point3f(0,0,0), Normalize(Vector3f(pCamera)));

    // Modify ray for depth of field
    if (lensRadius > 0) {
        // Sample point on lens
        Point2f pLens = lensRadius * ConcentricSampleDisk(sample.pLens);

        // Compute point on plane of focus
        Float ft = focalDistance / ray->d.z;
        Point3f pFocus = (*ray)(ft);

        // Update ray for effect of lens
        ray->o = Point3f(pLens.x, pLens.y, 0);
        ray->d = Normalize(pFocus - ray->o);
    }
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    ray->medium = medium;
    *ray = CameraToWorld(*ray);
    return 1;
  }

  DistortionCamera *CreateDistortionCamera(const ParamSet &params,
                                            const AnimatedTransform &cam2world,
                                            Film *film, const Medium *medium) {
    // for now extract the parameters for the distortion camera from the 
    // paramset, verify that they are there and go on to create the usual
    // perspective camera
    Float centerOffsetX = params.FindOneFloat("centerx", 0.0);
    Float centerOffsetY = params.FindOneFloat("centery", 0.0);
    std::string model = params.FindOneString("model", NO_MODEL);
    std::cout << "Model in Create: " << model << std::endl;
    int n;
    const Float *coeffPtr =  params.FindFloat("coefficients", &n);
    DistortionCamera::coeffVec coeffs;
    for (int i = 0; i < n; i++)
      coeffs.push_back(*(coeffPtr + i));
    std::cout << "Distortion Camera with " << model << " and " << n << " coeffs : ";
    for (DistortionCamera::coeffVec::iterator it = coeffs.begin(); it != coeffs.end(); it++)
      std::cout << *it << " ";
    std::cout << std::endl;

    // ------------- equivalent things from PerspectiveCamera ----------------
    Float shutteropen = params.FindOneFloat("shutteropen", 0.f);
    Float shutterclose = params.FindOneFloat("shutterclose", 1.f);
    if (shutterclose < shutteropen) {
        Warning("Shutter close time [%f] < shutter open [%f].  Swapping them.",
                shutterclose, shutteropen);
        std::swap(shutterclose, shutteropen);
    }
    Float lensradius = params.FindOneFloat("lensradius", 0.f);
    Float focaldistance = params.FindOneFloat("focaldistance", 1e6);
    Float frame = params.FindOneFloat(
        "frameaspectratio",
        Float(film->fullResolution.x) / Float(film->fullResolution.y));
    Bounds2f screen;
    if (frame > 1.f) {
        screen.pMin.x = -frame;
        screen.pMax.x = frame;
        screen.pMin.y = -1.f;
        screen.pMax.y = 1.f;
    } else {
        screen.pMin.x = -1.f;
        screen.pMax.x = 1.f;
        screen.pMin.y = -1.f / frame;
        screen.pMax.y = 1.f / frame;
    }
    int swi;
    const Float *sw = params.FindFloat("screenwindow", &swi);
    if (sw) {
        if (swi == 4) {
            screen.pMin.x = sw[0];
            screen.pMax.x = sw[1];
            screen.pMin.y = sw[2];
            screen.pMax.y = sw[3];
        } else
            Error("\"screenwindow\" should have four values");
    }
    Float fov = params.FindOneFloat("fov", 90.);
    Float halffov = params.FindOneFloat("halffov", -1.f);
    if (halffov > 0.f)
        // hack for structure synth, which exports half of the full fov
        fov = 2.f * halffov;
    // ------------- equivalent things from PerspectiveCamera end -------------

    return new DistortionCamera(cam2world, screen, shutteropen, shutterclose,
                                lensradius, focaldistance, fov, film, medium,
                                model, coeffs, centerOffsetX, centerOffsetY);
  }

}
