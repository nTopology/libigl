// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "random_points_on_mesh.h"
#include "doublearea.h"
#include "cumsum.h"
#include "histc.h"
#include <iostream>
#include <cassert>

template <typename DerivedV, typename DerivedF, typename DerivedB, typename DerivedFI>
IGL_INLINE void igl::random_points_on_mesh(
  const int n,
  const Eigen::MatrixBase<DerivedV > & V,
  const Eigen::MatrixBase<DerivedF > & F,
  Eigen::PlainObjectBase<DerivedB > & B,
  Eigen::PlainObjectBase<DerivedFI > & FI)
{
  using namespace Eigen;
  using namespace std;
  typedef typename DerivedV::Scalar Scalar;
  typedef Matrix<Scalar,Dynamic,1> VectorXs;
  VectorXs A;
  doublearea(V,F,A);
  // Should be traingle mesh. Although Turk's method 1 generalizes...
  assert(F.cols() == 3);
  VectorXs C;
  VectorXs A0(A.size()+1);
  A0(0) = 0;
  A0.bottomRightCorner(A.size(),1) = A;
  // Even faster would be to use the "Alias Table Method"
  cumsum(A0,1,C);
  const Scalar Cmax = C(C.size()-1);
  assert(Cmax > 0 && "Total surface area should be positive");
  // Why is this more accurate than `C /= C(C.size()-1)` ?
  for(int i = 0;i<C.size();i++) { C(i) = C(i)/Cmax; }
  const VectorXs R = (VectorXs::Random(n,1).array() + 1.)/2.;
  assert(R.minCoeff() >= 0);
  assert(R.maxCoeff() <= 1);
  histc(R,C,FI);
  FI = FI.array().min(int(F.rows()) - 1); // fix the bin when R(i) == 1 exactly
  const VectorXs S = (VectorXs::Random(n,1).array() + 1.)/2.;
  const VectorXs T = (VectorXs::Random(n,1).array() + 1.)/2.;
  B.resize(n,3);
  B.col(0) = 1.-T.array().sqrt();
  B.col(1) = (1.-S.array()) * T.array().sqrt();
  B.col(2) = S.array() * T.array().sqrt();
}

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedB,
  typename DerivedFI,
  typename DerivedX>
IGL_INLINE void igl::random_points_on_mesh(
  const int n,
  const Eigen::MatrixBase<DerivedV > & V,
  const Eigen::MatrixBase<DerivedF > & F,
  Eigen::PlainObjectBase<DerivedB > & B,
  Eigen::PlainObjectBase<DerivedFI > & FI,
  Eigen::PlainObjectBase<DerivedX> & X)
{
  random_points_on_mesh(n,V,F,B,FI);
  X = DerivedX::Zero(B.rows(),V.cols());
  for(int x = 0;x<B.rows();x++)
  {
    for(int b = 0;b<B.cols();b++)
    {
        auto fi = FI(x);
        auto f = F(fi, b);
        auto bval = B(x, b);
        X.row(x) += bval*V.row(f);
    }
  }
}

template <typename DerivedV, typename DerivedF, typename ScalarB, typename DerivedFI>
IGL_INLINE void igl::random_points_on_mesh(
  const int n,
  const Eigen::MatrixBase<DerivedV > & V,
  const Eigen::MatrixBase<DerivedF > & F,
  Eigen::SparseMatrix<ScalarB > & B,
  Eigen::PlainObjectBase<DerivedFI > & FI)
{
  using namespace Eigen;
  using namespace std;
  Matrix<ScalarB,Dynamic,3> BC;
  random_points_on_mesh(n,V,F,BC,FI);
  vector<Triplet<ScalarB> > BIJV;
  BIJV.reserve(n*3);
  for(int s = 0;s<n;s++)
  {
    for(int c = 0;c<3;c++)
    {
      assert(FI(s) < F.rows());
      assert(FI(s) >= 0);
      const int v = F(FI(s),c);
      BIJV.push_back(Triplet<ScalarB>(s,v,BC(s,c)));
    }
  }
  B.resize(n,V.rows());
  B.reserve(n*3);
  B.setFromTriplets(BIJV.begin(),BIJV.end());
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template void igl::random_points_on_mesh<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 3, 1, -1, 3> >(int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&);
template void igl::random_points_on_mesh<Eigen::Matrix<float, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<float, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<float, -1, 3, 1, -1, 3> >(int, Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> >&);
template void igl::random_points_on_mesh<Eigen::Matrix<float, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, float, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(int, Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::SparseMatrix<float, 0, int>&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::random_points_on_mesh<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::random_points_on_mesh<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::SparseMatrix<double, 0, int>&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::random_points_on_mesh<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(int, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::SparseMatrix<double, 0, int>&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
