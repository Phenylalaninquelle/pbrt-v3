#include "cameras/perspective.h"
#include "cameras/distortion.h"
#include "paramset.h"
#include <vector>
#include <unordered_map>
#include "sampling.h"

#define NO_MODEL "MODEL_NOT_FOUND"
#define POLY_DEGREE 5

namespace pbrt {

  // add new models and their number of coefficients here
  std::unordered_map<std::string, int> DistortionCamera::numCoeffsForModel = { {"poly3lensfun", 1},
                                                                               {"poly5lensfun", 2},
                                                                               {"ptlens", 3} };
  DistortionCamera::DistortionCamera(const AnimatedTransform &CameraToWorld,
                                     const Bounds2f &screenWindow, Float shutterOpen,
                                     Float shutterClose, Float lensRadius, Float focalDistance,
                                     Float fov, Film *film, const Medium *medium,
                                     std::string distortionModel, DistortionCamera::coeffVec coeffs,
                                     Float centerOffsetX, Float centerOffsetY)
    : ProjectiveCamera(CameraToWorld, Perspective(fov, 1e-2f, 1000.f),
                       screenWindow, shutterOpen, shutterClose, lensRadius,
                       focalDistance, film, medium),
      distortionModel(distortionModel), coeffs(coeffs),
      centerOffsetX(centerOffsetX), centerOffsetY(centerOffsetY), fittedCoeffs(){

        /* 
         * check for valid input of model and params
         * then set or compute the fitted coefficients accordingly
         */
        const coeffVec PERSPECTIVE_COEFFS = coeffVec({0, 1});
        if (distortionModel == NO_MODEL) {
          Error("No model for lense distortion given. Rendering will be done without distortion, equivalent to PerspectiveCamera.");
          fittedCoeffs = PERSPECTIVE_COEFFS;
        }
        else {
          auto modelIter = numCoeffsForModel.find(distortionModel);
          if (modelIter == numCoeffsForModel.end()) {
            Error("Model %s is unsupported (did you make a typo?). Rendering will be done without distortion, equivalent to PerspectiveCamera.",
                  distortionModel.c_str());
            fittedCoeffs = PERSPECTIVE_COEFFS;
          }
          else {
            int requiredNum = modelIter->second;
            int givenNum = coeffs.size();
            if (givenNum > requiredNum) {
              Error("Model %s requires %d coefficients, but %d provided. I will use only the first %d.", 
                    distortionModel.c_str(), requiredNum, givenNum, requiredNum);
              coeffs = coeffVec(coeffs.begin(), coeffs.begin() + requiredNum);
            }
            else if (givenNum < requiredNum) {
              Error("Model %s requires %d coefficients, but only %d provided. I will assume the missing %d to be zero.",
                    distortionModel.c_str(), requiredNum, givenNum, requiredNum - givenNum);
              coeffVec zeros = coeffVec(requiredNum - givenNum, 0);
              coeffs.insert(coeffs.end(), zeros.begin(), zeros.end());
            }
            fittedCoeffs = InvertDistortion(coeffs, POLY_DEGREE);
          }
        }
        SetupImageNormalization(film->fullResolution.x, film->fullResolution.y, centerOffsetX, centerOffsetY);
    }

  void DistortionCamera::SetupImageNormalization(Float xRes, Float yRes,
                                                 Float centerOffsetX, Float centerOffsetY) {
        /* ---- transform image center to pixel coordinates ----*/
        /*
         * Image center in Lensfun database is given as an offset from the ideal
         * center of the viewframe. The values are normalised  so that the smaller
         * image dimension has the value 2 (yay, another coordinate system). 
         * So, to transform from this offset coordinates to Pixelcoordinates, we 
         * need to scale the offset by yRes/2 (assuming, that yRes < xRes). 
         */
        Float centerOffsetScaleFactor = (yRes < xRes) ? yRes/2 : xRes/2;
        Transform OffsetToPixelCoordinates = Scale(centerOffsetScaleFactor, centerOffsetScaleFactor, 1.);
        Point3f CenterOffsetPixels = OffsetToPixelCoordinates(Point3f(centerOffsetX, centerOffsetY, 0));

        // find biggest center to corner distance based on sign of the offsets:
        Point3f farthestCorner = Point3f(0, 0, 1);
        if (CenterOffsetPixels.x < 0)
          farthestCorner.x = xRes;
        if (CenterOffsetPixels.y < 0)
          farthestCorner.y = yRes;

        Point3f shiftedCenter = Point3f(xRes / 2 + CenterOffsetPixels.x,
                                        yRes / 2 + CenterOffsetPixels.y, 1);
        Float cornerRadius = (Vector3f(farthestCorner) - Vector3f(shiftedCenter)).Length();
        NormalizeToCornerRadius = Scale(1. / cornerRadius, 1. / cornerRadius, 1.);
        Denormalize = Inverse(NormalizeToCornerRadius);
        imageCenterNormalized = NormalizeToCornerRadius(shiftedCenter);
  }

  DistortionCamera::coeffVec DistortionCamera::InvertDistortion(DistortionCamera::coeffVec coeffs,
                                                                int polyDegree) const {
    // create x vector for sampling distortion values
    // and y vector with sampled values
    int sampleSize = 1000;
    Float scale(sampleSize);
    coeffVec x(sampleSize);
    coeffVec y(sampleSize);

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

  Point3f DistortionCamera::CalculateRayStartingpoint(const CameraSample &sample) const {
    Float centerX = imageCenterNormalized.x;
    Float centerY = imageCenterNormalized.y;
    Point3f pNormalized = NormalizeToCornerRadius(Point3f(sample.pFilm.x, sample.pFilm.y, 0));
    Float radius = sqrt(pow(pNormalized.x - centerX, 2) + pow(pNormalized.y - centerY, 2));
    Float rNew = evalPolynomial(fittedCoeffs, coeffVec({radius}))[0];
    Float rRatio = rNew / radius;
    return Denormalize(Point3f(pNormalized.x * rRatio + centerX * (1 - rRatio),
                               pNormalized.y * rRatio + centerY * (1 - rRatio), 0));
  }


  Float DistortionCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
    ProfilePhase prof(Prof::GenerateCameraRay);
    Point3f pCamera = RasterToCamera(CalculateRayStartingpoint(sample));
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
    Float centerOffsetX = params.FindOneFloat("centerx", 0.0);
    Float centerOffsetY = params.FindOneFloat("centery", 0.0);
    std::string model = params.FindOneString("model", NO_MODEL);
    int n;
    DistortionCamera::coeffVec coeffs;
    const Float *coeffPtr =  params.FindFloat("coefficients", &n);
    if (coeffPtr) {
      for (int i = 0; i < n; i++)
        coeffs.push_back(*(coeffPtr + i));
    }
    else {
      coeffs = DistortionCamera::coeffVec({0, 1});
    }

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

    return new DistortionCamera(cam2world, screen, shutteropen, shutterclose,
                                lensradius, focaldistance, fov, film, medium,
                                model, coeffs, centerOffsetX, centerOffsetY);
  }

}
