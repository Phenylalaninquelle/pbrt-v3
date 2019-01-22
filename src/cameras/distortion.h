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


namespace pbrt {

  class DistortionCamera : public ProjectiveCamera {
    
    public:
      typedef std::vector<Float> coeffVec;
      DistortionCamera(const AnimatedTransform &CameraToWorld,
                       const Bounds2f &screenWindow, Float shutterOpen,
                       Float shutterClose, Float lensRadius, Float focalDistance,
                       Float fov, Film *film, const Medium *medium,
                       std::string distortion_model, coeffVec coeffs);
      Float GenerateRay(const CameraSample &sample, Ray *ray) const;
      //Float GenerateRayDifferential(const CameraSample &sample,
                                    //RayDifferential *ray) const;
    private:
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
