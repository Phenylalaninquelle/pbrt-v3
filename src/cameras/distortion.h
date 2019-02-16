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
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>

namespace pbrt {

  class DistortionCamera : public ProjectiveCamera {
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
      coeffVec InvertDistortion(coeffVec coeffs, int poly_degree) const;
      Point3f CalculateRayStartingpoint(const CameraSample& sample) const;
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
    namespace ublas = boost::numeric::ublas;

    if (x.size() != y.size())
      throw "x and y of unequal size";

    degree++;

    // construct system of equations
    int n = x.size();
    ublas::matrix<T> vandMatrix(n, degree);
    ublas::matrix<T> Y(n, 1);

    for (int i = 0; i < n; i++)
      Y(i, 0) = y[i];
    for (int row = 0; row < n; row++) {
      T val = 1.0f;
      for (int col = 0; col < degree; col++) {
        vandMatrix(row, col) = val;
        val *= x[row];
      }
    }

    ublas::matrix<T> vandTransposed(ublas::trans(vandMatrix));
    ublas::matrix<T> vandProd(prec_prod(vandTransposed, vandMatrix));
    ublas::matrix<T> solvedCoeffs(ublas::prec_prod(vandTransposed, Y));

    // solve equations
    ublas::permutation_matrix<int> perm(vandProd.size1());
    const int isSingular = ublas::lu_factorize(vandProd, perm);
    if (isSingular)
      throw "System matrix is singular. Did you pass only unique values in x?";
    ublas::lu_substitute(vandProd, perm, solvedCoeffs);

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
    return radius * (1 - coeffs[0] + coeffs[0] * pow(radius, 2));
  }

  inline Float ModelPoly5LensFun(const Float radius, const DistortionCamera::coeffVec coeffs) {
    return radius * (1 + coeffs[0] * pow(radius, 2) + coeffs[1] * pow(radius, 4));
  }

  inline Float ModelPTLens(const Float radius, const DistortionCamera::coeffVec coeffs) {
    return radius * (coeffs[0] * pow(radius, 3) + coeffs[1] * pow(radius, 2), + coeffs[2] * radius + 1 - coeffs[0] - coeffs[1] - coeffs[2]);
  }
}

#endif // PBRT_CAMERAS_DISTORTION_H
