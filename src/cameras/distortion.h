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


namespace pbrt {
  // for now the distortionCamera factory creates a perspective camera 
  // for testing purposes
  PerspectiveCamera *CreateDistortionCamera(const ParamSet &params,
                                            const AnimatedTransform &cam2world,
                                            Film *film, const Medium *medium);
}

#endif // PBRT_CAMERAS_DISTORTION_H
