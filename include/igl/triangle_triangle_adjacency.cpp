// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "triangle_triangle_adjacency.h"
#include "all_edges.h"
#include "unique_simplices.h"
#include "unique_edge_map.h"
#include <algorithm>
#include <iostream>

// Extract the face adjacencies
template <typename DerivedF, typename TTT_type, typename DerivedTT>
IGL_INLINE void igl::triangle_triangle_adjacency_extractTT(
  const Eigen::PlainObjectBase<DerivedF>& F,
  std::vector<std::vector<TTT_type> >& TTT,
  Eigen::PlainObjectBase<DerivedTT>& TT,NTInterrupter* mInterrupter)
{
  TT.setConstant((int)(F.rows()),F.cols(),-1);

  for(int i=1;i<(int)TTT.size();++i)
  {
    if(NTInterrupter::wasInterrupted(mInterrupter))return;
    std::vector<int>& r1 = TTT[i-1];
    std::vector<int>& r2 = TTT[i];
    if ((r1[0] == r2[0]) && (r1[1] == r2[1]))
    {
      TT(r1[2],r1[3]) = r2[2];
      TT(r2[2],r2[3]) = r1[2];
    }
  }
}

template <typename DerivedF, typename DerivedTT>
IGL_INLINE void igl::triangle_triangle_adjacency(
  const Eigen::PlainObjectBase<DerivedF>& F,
  Eigen::PlainObjectBase<DerivedTT>& TT,NTInterrupter* mInterrupter)
{
  DerivedTT TTi;
  return triangle_triangle_adjacency(F,TT,TTi,mInterrupter);
}

template <typename DerivedF, typename TTT_type>
IGL_INLINE void igl::triangle_triangle_adjacency_preprocess(
  const Eigen::PlainObjectBase<DerivedF>& F,
  std::vector<std::vector<TTT_type> >& TTT,NTInterrupter* mInterrupter)
{
  for(int f=0;f<F.rows();++f){
    if(NTInterrupter::wasInterrupted(mInterrupter))return;
    for (int i=0;i<F.cols();++i)
    {
      // v1 v2 f ei
      if(NTInterrupter::wasInterrupted(mInterrupter))return;
      int v1 = F(f,i);
      int v2 = F(f,(i+1)%F.cols());
      if (v1 > v2) std::swap(v1,v2);
      std::vector<int> r(4);
      r[0] = v1; r[1] = v2;
      r[2] = f;  r[3] = i;
      TTT.push_back(r);
    }
  }
  std::sort(TTT.begin(),TTT.end());
}

// Extract the face adjacencies indices (needed for fast traversal)
template <typename DerivedF, typename TTT_type, typename DerivedTTi>
IGL_INLINE void igl::triangle_triangle_adjacency_extractTTi(
  const Eigen::PlainObjectBase<DerivedF>& F,
  std::vector<std::vector<TTT_type> >& TTT,
  Eigen::PlainObjectBase<DerivedTTi>& TTi,NTInterrupter* mInterrupter)
{
  TTi.setConstant((int)(F.rows()),F.cols(),-1);

  for(int i=1;i<(int)TTT.size();++i)
  {
    if(NTInterrupter::wasInterrupted(mInterrupter))return;
    std::vector<int>& r1 = TTT[i-1];
    std::vector<int>& r2 = TTT[i];
    if ((r1[0] == r2[0]) && (r1[1] == r2[1]))
    {
      TTi(r1[2],r1[3]) = r2[3];
      TTi(r2[2],r2[3]) = r1[3];
    }
  }
}

// Compute triangle-triangle adjacency with indices
template <typename DerivedF, typename DerivedTT, typename DerivedTTi>
IGL_INLINE void igl::triangle_triangle_adjacency(
  const Eigen::PlainObjectBase<DerivedF>& F,
  Eigen::PlainObjectBase<DerivedTT>& TT,
  Eigen::PlainObjectBase<DerivedTTi>& TTi,NTInterrupter* mInterrupter)
{
  std::vector<std::vector<int> > TTT;
  triangle_triangle_adjacency_preprocess(F,TTT,mInterrupter);
  triangle_triangle_adjacency_extractTT(F,TTT,TT,mInterrupter);
  triangle_triangle_adjacency_extractTTi(F,TTT,TTi,mInterrupter);
}

template <
  typename DerivedF,
  typename TTIndex,
  typename TTiIndex>
  IGL_INLINE void igl::triangle_triangle_adjacency(
    const Eigen::PlainObjectBase<DerivedF> & F,
    std::vector<std::vector<std::vector<TTIndex> > > & TT,
    std::vector<std::vector<std::vector<TTiIndex> > > & TTi,NTInterrupter* mInterrupter)
{
  return triangle_triangle_adjacency(F,true,TT,TTi,mInterrupter);
}

template <
  typename DerivedF,
  typename TTIndex>
  IGL_INLINE void igl::triangle_triangle_adjacency(
    const Eigen::PlainObjectBase<DerivedF> & F,
    std::vector<std::vector<std::vector<TTIndex> > > & TT,NTInterrupter* mInterrupter)
{
  std::vector<std::vector<std::vector<TTIndex> > > not_used;
  return triangle_triangle_adjacency(F,false,TT,not_used,mInterrupter);
}

template <
  typename DerivedF,
  typename TTIndex,
  typename TTiIndex>
  IGL_INLINE void igl::triangle_triangle_adjacency(
    const Eigen::PlainObjectBase<DerivedF> & F,
    const bool construct_TTi,
    std::vector<std::vector<std::vector<TTIndex> > > & TT,
    std::vector<std::vector<std::vector<TTiIndex> > > & TTi,NTInterrupter* mInterrupter)
{
  using namespace Eigen;
  using namespace std;
  assert(F.cols() == 3 && "Faces must be triangles");
  // number of faces
  typedef typename DerivedF::Index Index;
  typedef Matrix<typename DerivedF::Scalar,Dynamic,2> MatrixX2I;
  typedef Matrix<typename DerivedF::Index,Dynamic,1> VectorXI;
  MatrixX2I E,uE;
  VectorXI EMAP;
  vector<vector<Index> > uE2E;
  if(NTInterrupter::startSection(mInterrupter, 0.5)) return;
  unique_edge_map(F,E,uE,EMAP,uE2E, mInterrupter);
  if(NTInterrupter::endSection(mInterrupter)) return;
  if(NTInterrupter::startSection(mInterrupter, 0.5)) return;
  triangle_triangle_adjacency(E,EMAP,uE2E,construct_TTi,TT,TTi,mInterrupter);
  if(NTInterrupter::endSection(mInterrupter)) return;
}

template <
  typename DerivedE,
  typename DerivedEMAP,
  typename uE2EType,
  typename TTIndex,
  typename TTiIndex>
  IGL_INLINE void igl::triangle_triangle_adjacency(
    const Eigen::PlainObjectBase<DerivedE> & E,
    const Eigen::PlainObjectBase<DerivedEMAP> & EMAP,
    const std::vector<std::vector<uE2EType> > & uE2E,
    const bool construct_TTi,
    std::vector<std::vector<std::vector<TTIndex> > > & TT,
    std::vector<std::vector<std::vector<TTiIndex> > > & TTi,NTInterrupter* mInterrupter)
{
  using namespace std;
  using namespace Eigen;
  typedef typename DerivedE::Index Index;
  const size_t m = E.rows()/3;
  assert((size_t)E.rows() == m*3 && "E should come from list of triangles.");
  // E2E[i] --> {j,k,...} means face edge i corresponds to other faces edges j
  // and k
  TT.resize (m,vector<vector<TTIndex> >(3));
  if(construct_TTi)
  {
    TTi.resize(m,vector<vector<TTiIndex> >(3));
  }

  // No race conditions because TT*[f][c]'s are in bijection with e's
  // Minimum number of iterms per openmp thread
  //const size_t num_e = E.rows();
# ifndef IGL_OMP_MIN_VALUE
#   define IGL_OMP_MIN_VALUE 1000
# endif
# pragma omp parallel for if (m>IGL_OMP_MIN_VALUE)
  // Slightly better memory access than loop over E
  for(Index f = 0;f<(Index)m;f++)
  {
    if(NTInterrupter::wasInterrupted(mInterrupter, double(f)/double(m)))return;
    for(Index c = 0;c<3;c++)
    {
      if(NTInterrupter::wasInterrupted(mInterrupter))return;
      const Index e = f + m*c;
      //const Index c = e/m;
      const vector<uE2EType> & N = uE2E[EMAP(e)];
      for(const auto & ne : N)
      {
        const Index nf = ne%m;
        // don't add self
        if(nf != f)
        {
          TT[f][c].push_back(nf);
          if(construct_TTi)
          {
            const Index nc = ne/m;
            TTi[f][c].push_back(nc);
          }
        }
      }
    }
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
// generated by autoexplicit.sh
template void igl::triangle_triangle_adjacency<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
// generated by autoexplicit.sh
template void igl::triangle_triangle_adjacency<Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&);
template void igl::triangle_triangle_adjacency<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::triangle_triangle_adjacency<Eigen::Matrix<int, -1, 2, 0, -1, 2>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, long, long, long>(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 2, 0, -1, 2> > const&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> > const&, std::vector<std::vector<long, std::allocator<long> >, std::allocator<std::vector<long, std::allocator<long> > > > const&, bool, std::vector<std::vector<std::vector<long, std::allocator<long> >, std::allocator<std::vector<long, std::allocator<long> > > >, std::allocator<std::vector<std::vector<long, std::allocator<long> >, std::allocator<std::vector<long, std::allocator<long> > > > > >&, std::vector<std::vector<std::vector<long, std::allocator<long> >, std::allocator<std::vector<long, std::allocator<long> > > >, std::allocator<std::vector<std::vector<long, std::allocator<long> >, std::allocator<std::vector<long, std::allocator<long> > > > > >&);
template void igl::triangle_triangle_adjacency<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, unsigned long, int, int>(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, std::vector<std::vector<unsigned long, std::allocator<unsigned long> >, std::allocator<std::vector<unsigned long, std::allocator<unsigned long> > > > const&, bool, std::vector<std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > >, std::allocator<std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > > >&, std::vector<std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > >, std::allocator<std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > > >&);
#endif
