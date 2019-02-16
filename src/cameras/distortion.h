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
    public:
      static std::unordered_map<std::string, int> numCoeffsForModel;
      typedef std::vector<Float> coeffVec;
      DistortionCamera(const AnimatedTransform &CameraToWorld,
                       const Bounds2f &screenWindow, Float shutterOpen,
                       Float shutterClose, Float lensRadius, Float focalDistance,
                       Float fov, Film *film, const Medium *medium,
                       std::string distortionModel, coeffVec coeffs,
                       Float centerOffsetX, Float centerOffsetY);
      Float GenerateRay(const CameraSample &sample, Ray *ray) const;
    private:
      coeffVec InvertDistortion(coeffVec coeffs, int poly_degree);
      Point3f CalculateRayStartpoint(const CameraSample& sample) const;
      void SetupImageNormalization(Float xRes, Float yRes, Float centerOffsetX, Float centerOffsetY);
      std::string distortionModel;
      coeffVec coeffs;
      Transform NormalizeToCornerRadius, Denormalize;
      Point3f imageCenterNormalized;
      coeffVec fittedCoeffs;
      int centerOffsetX, centerOffsetY;
  };

  DistortionCamera *CreateDistortionCamera(const ParamSet &params,
                                            const AnimatedTransform &cam2world,
                                            Film *film, const Medium *medium);

  template<typename T> std::vector<T> fitPolyCoeffs(const std::vector<T>& x,
                                                      const std::vector<T>& y,
                                                      int degree) {
    using namespace boost::numeric::ublas;

    assert(x.size() == y.size());

    degree++;

    // construct system of equations
    int n = x.size();
    matrix<T> vandMatrix(n, degree);
    matrix<T> Y(n, 1);

    for (int i = 0; i < n; i++)
      Y(i, 0) = y[i];
    for (int row = 0; row < n; row++) {
      T val = 1.0f;
      for (int col = 0; col < degree; col++) {
        vandMatrix(row, col) = val;
        val *= x[row];
      }
    }

    matrix<T> vandTransposed(trans(vandMatrix));
    matrix<T> vandProd(prec_prod(vandTransposed, vandMatrix));
    matrix<T> solvedCoeffs(prec_prod(vandTransposed, Y));

    // solve equations
    permutation_matrix<int> perm(vandProd.size1());
    const int singular = lu_factorize(vandProd, perm);
    assert(singular == 0);
    lu_substitute(vandProd, perm, solvedCoeffs);

    return std::vector<T>(solvedCoeffs.data().begin(), solvedCoeffs.data().end());
  }

  template<typename T> std::vector<T> evalPolynomial(const std::vector<T>& coeffs,
                                                      const std::vector<T>& x) {
    std::vector<T> result(x.size());
    for (unsigned int i = 0; i < x.size(); i++) {
      T xTmp = 1;
      T yTmp = 0;
      for (unsigned int j = 0; j < coeffs.size(); j++) {
        yTmp += coeffs[j] * xTmp;
        xTmp *= x[i];
      }
      result[i] = yTmp;
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
