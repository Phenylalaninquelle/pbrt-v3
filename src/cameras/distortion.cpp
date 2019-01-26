//TODO: license??

#include "cameras/perspective.h"
#include "cameras/distortion.h"
#include "paramset.h"
#include <iostream> // TODO:remove
#include <vector>
#include <cmath>
#include "sampling.h"

#define NO_MODEL "MODEL_NOT_FOUND"

namespace pbrt {

  std::unordered_map<std::string, int> DistortionCamera::num_coeffs_for_model = 
    { {"poly3lensfun", 1}, };

  std::unordered_set<std::string> DistortionCamera::supported_models = 
    {"poly3lensfun",};

  DistortionCamera::DistortionCamera(const AnimatedTransform &CameraToWorld,
                                     const Bounds2f &screenWindow, Float shutterOpen,
                                     Float shutterClose, Float lensRadius, Float focalDistance,
                                     Float fov, Film *film, const Medium *medium,
                                     std::string distortion_model, DistortionCamera::coeffVec coeffs)
    : ProjectiveCamera(CameraToWorld, Perspective(fov, 1e-2f, 1000.f),
                       screenWindow, shutterOpen, shutterClose, lensRadius,
                       focalDistance, film, medium),
      distortion_model(distortion_model),
      coeffs(coeffs),
      fitted_coeffs(){

        /* -------- Define transformations for ray generation --------*/
        //TODO: this normalisation makes no sense?
        Float xRes = film->fullResolution.x;
        Float yRes = film->fullResolution.y;
        //if (xRes >= yRes)
          //RasterToNDC = Scale(1. / xRes, 1. / xRes, 1.);
        //else
          //RasterToNDC = Scale(1. / yRes, 1. / yRes, 1.);
        RasterToNDC = Scale(1. / xRes, 1. / yRes, 1.);
        NDCToRaster = Inverse(RasterToNDC);
        /* -----------------------------------------------------------*/

        /* ------------ Validate given model ---------------------------*/
        // if no model is given either stop or do something smart
        if (distortion_model == NO_MODEL) {
          Error("No Model given.");
        }
        else {
          // check if given model is supported
          if (supported_models.find(distortion_model) == supported_models.end()) {
            Error("Model %s is unsupported. Abort.", distortion_model.c_str());
            // TODO: how do you stop this thing?
          }
        }
        /* -----------------------------------------------------------*/
        // TODO: remove hardcoded polnomial degree
        fitted_coeffs = InvertDistortion(distortion_model, coeffs, 5);
        //fitted_coeffs = coeffVec({0, 1});
    }

  DistortionCamera::coeffVec DistortionCamera::InvertDistortion(std::string distortion_model,
                                                                DistortionCamera::coeffVec coeffs,
                                                                int poly_degree) {
    // create x vector for sampling distortion values
    // and y vector with sampled values
    int sample_size = 100;
    
    std::vector<Float> x(sample_size);
    std::vector<Float> y(sample_size);
    // fill sample vectors according to the model used
    if (distortion_model == "poly3lensfun") {
      Float k = coeffs[0];
      Float scale(sample_size);
      for (int i = 0; i < sample_size; i++) {
        x[i] = i / scale;
        y[i] = ModelPoly3LensFun(x[i], k);
      }
    }

    coeffVec poly_coeffs = fit_poly_coeffs(y, x, poly_degree);
    return poly_coeffs;
  }

  // the fitted radial model is applied here
  Point3f DistortionCamera::CalculateRayStartpoint(const CameraSample &sample) const {
    Point3f pFilm = Point3f(sample.pFilm.x, sample.pFilm.y, 0);
    Point3f pNDC = RasterToNDC(pFilm);
    Float xCenter = .5, yCenter = .5;
    Float radius = sqrt(pow(pNDC.x - xCenter, 2) + pow(pNDC.y - yCenter, 2));
    Float r_new = eval_polynomial(fitted_coeffs, std::vector<Float>({radius}))[0];
    Float r_ratio = r_new / radius;
    return NDCToRaster(Point3f(pNDC.x * r_ratio + xCenter * (1 - r_ratio),
                               pNDC.y * r_ratio + yCenter * (1 - r_ratio), 0));
  }

  Float DistortionCamera::ModelPoly3LensFun(const Float radius, const Float k) const {
    return radius * (1 - k + k * pow(radius, 2));
  }

  // Compute the distorted position of a Point on the film in raster coordinates
  //Point3f DistortionCamera::ModelPoly3LensFun(const Point3f pFilm, const Float k) const {
    //Point3f pNDC = RasterToNDC(pFilm);
    //Float xCenter = .5, yCenter = .5;
    //Float radius = sqrt(pow(pNDC.x - xCenter, 2) + pow(pNDC.y - yCenter, 2));
    //Float K = 1 - k + k * pow(radius, 2);
    //return NDCToRaster(Point3f(xCenter * (1 - K) + pNDC.x * K,
                               //yCenter * (1 - K) + pNDC.y * K,
                               //0));
  //}

  Float DistortionCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
    ProfilePhase prof(Prof::GenerateCameraRay);
    Point3f pCamera = RasterToCamera(CalculateRayStartpoint(sample));
    //Point3f pCamera = Point3f(1,1,0);
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
                                model, coeffs);
  }

}
