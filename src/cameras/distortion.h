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
#include <cassert>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>


namespace pbrt {

  class DistortionCamera : public ProjectiveCamera {

    // NOTE: when adding a new distortion model, add it and the number of 
    // its coefficients here
    static std::unordered_map<std::string, int> num_coeffs_for_model;
    public:
      typedef std::vector<Float> coeffVec;
      DistortionCamera(const AnimatedTransform &CameraToWorld,
                       const Bounds2f &screenWindow, Float shutterOpen,
                       Float shutterClose, Float lensRadius, Float focalDistance,
                       Float fov, Film *film, const Medium *medium,
                       std::string distortion_model, coeffVec coeffs);
      Float GenerateRay(const CameraSample &sample, Ray *ray) const;
    private:
      coeffVec InvertDistortion(coeffVec coeffs, int poly_degree);
      Point3f CalculateRayStartpoint(const CameraSample& sample) const;
      std::string distortion_model;
      coeffVec coeffs;
      Transform RasterToNDC, NDCToRaster;
      coeffVec fitted_coeffs;
  };

  DistortionCamera *CreateDistortionCamera(const ParamSet &params,
                                            const AnimatedTransform &cam2world,
                                            Film *film, const Medium *medium);

  template<typename T> std::vector<T> fit_poly_coeffs(const std::vector<T>& x,
                                                      const std::vector<T>& y,
                                                      int degree) {
    using namespace boost::numeric::ublas;

    assert(x.size() == y.size());

    degree++;

    // construct system of equations
    int n = x.size();
    matrix<T> vandermonde_matrix(n, degree);
    matrix<T> Y(n, 1);

    for (int i = 0; i < n; i++)
      Y(i, 0) = y[i];
    for (int row = 0; row < n; row++) {
      T val = 1.0f;
      for (int col = 0; col < degree; col++) {
        vandermonde_matrix(row, col) = val;
        val *= x[row];
      }
    }

    matrix<T> vandermonde_transposed(trans(vandermonde_matrix));
    matrix<T> vandermonde_product(prec_prod(vandermonde_transposed, vandermonde_matrix));
    matrix<T> solved_coeffs(prec_prod(vandermonde_transposed, Y));

    // solve equations
    permutation_matrix<int> perm(vandermonde_product.size1());
    const int singular = lu_factorize(vandermonde_product, perm);
    assert(singular == 0);
    lu_substitute(vandermonde_product, perm, solved_coeffs);

    return std::vector<T>(solved_coeffs.data().begin(), solved_coeffs.data().end());
  }

  template<typename T> std::vector<T> eval_polynomial(const std::vector<T>& coeffs,
                                                      const std::vector<T>& x) {
    std::vector<T> result(x.size());
    for (unsigned int i = 0; i < x.size(); i++) {
      T x_tmp = 1;
      T y_tmp = 0;
      for (unsigned int j = 0; j < coeffs.size(); j++) {
        y_tmp += coeffs[j] * x_tmp;
        x_tmp *= x[i];
      }
      result[i] = y_tmp;
    }
    return result;
  }

  inline Float ModelPoly3LensFun(const Float radius, const DistortionCamera::coeffVec coeffs) {
    Warning("Poly3");
    return radius * (1 - coeffs[0] + coeffs[0] * pow(radius, 2));
  }

  inline Float ModelPoly5LensFun(const Float radius, const DistortionCamera::coeffVec coeffs) {
    Warning("Poly5");
    return radius * (1 + coeffs[0] * pow(radius, 2) + coeffs[1] * pow(radius, 4));
  }

  inline Float ModelPTLens(const Float radius, const DistortionCamera::coeffVec coeffs) {
    Warning("Ptlens");
    return radius * (coeffs[0] * pow(radius, 3) + coeffs[1] * pow(radius, 2), + coeffs[2] * radius + 1 - coeffs[0] - coeffs[1] - coeffs[2]);
  }
}

#endif // PBRT_CAMERAS_DISTORTION_H
