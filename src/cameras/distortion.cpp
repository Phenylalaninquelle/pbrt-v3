//TODO: license??

#include "cameras/perspective.h"
#include "cameras/distortion.h"
#include "paramset.h"
#include <iostream> // TODO:remove


namespace pbrt {

  PerspectiveCamera *CreateDistortionCamera(const ParamSet &params,
                                            const AnimatedTransform &cam2world,
                                            Film *film, const Medium *medium) {
    // for now extract the parameters for the distortion camera from the 
    // paramset, verify that they are there and go on to create the usual
    // perspective camera
    auto model = params.FindOneString("model", "model_not_found");
    int n;
    params.FindFloat("coefficients", &n);
    std::cout << "Distortion Camera with " << model << " and " << n << " coeffs" << std::endl;
    return CreatePerspectiveCamera(params, cam2world, film, medium);
  }

}
