// TODO: license?


#if defined (_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CAMERAS_DISTORTION_H
#define PBRT_CAMERAS_DISTORTION_H

#include "pbrt.h"
#include "camera.h"
#include "perspective.h"
#include "film.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>

namespace pbrt {

  class DistortionCamera : public ProjectiveCamera {
    
    //TODO: static variable with the suported models
    static std::unordered_set<std::string> supported_models;
    static std::unordered_map<std::string, int> num_coeffs_for_model;
    //std::unordered_map<std::string, void>function_for_model;
    public:
      typedef std::vector<Float> coeffVec;
      DistortionCamera(const AnimatedTransform &CameraToWorld,
                       const Bounds2f &screenWindow, Float shutterOpen,
                       Float shutterClose, Float lensRadius, Float focalDistance,
                       Float fov, Film *film, const Medium *medium,
                       std::string distortion_model, coeffVec coeffs);
      Float GenerateRay(const CameraSample &sample, Ray *ray) const;
    private:
      Point3f ApplyDistortionModel(const CameraSample &sample) const;
      Point3f ModelPoly3LensFun(const Point3f pFilm, const Float k) const;
      std::string distortion_model;
      coeffVec coeffs;
      Transform RasterToNDC, NDCToRaster;

  };
  // for now the distortionCamera factory creates a perspective camera 
  // for testing purposes
  DistortionCamera *CreateDistortionCamera(const ParamSet &params,
                                            const AnimatedTransform &cam2world,
                                            Film *film, const Medium *medium);
}

#endif // PBRT_CAMERAS_DISTORTION_H
