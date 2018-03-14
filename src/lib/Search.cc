#include <nabo/nabo.h>
#include <Eigen/Core>

extern "C" {
  int SearchNeighbors(const double **Position, int *i_pair, int *j_pair, double *rij2, const int NumberOfFluidParticles, const int NumberOfParticles, int ParticleBoundary, const int NumberOfNeighbors, const double dist)
  {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> M;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> q;
    Eigen::MatrixXi indices;
    Eigen::MatrixXd dists;

    M.resize(2, NumberOfParticles);
    q.resize(2, NumberOfFluidParticles);
    indices.resize(NumberOfNeighbors, q.cols());
    dists.resize(NumberOfNeighbors, q.cols());
    for (int s = 0; s < NumberOfParticles; s++) {
      M(0, s) = Position[s][0];
      M(1, s) = Position[s][1];
      // printf("%.5lf %.5lf\n", Position[s][0], Position[s][1]);
    }

    for (int s = 0; s < NumberOfFluidParticles; s++) {
      q(0, s) = Position[s][0];
      q(1, s) = Position[s][1];
    }

    Nabo::NNSearchD* nns = Nabo::NNSearchD::createKDTreeTreeHeap(M);
    nns->knn(q, indices, dists, NumberOfNeighbors);
    int niac = 0;

    for (int s = 0; s < NumberOfFluidParticles; s++) {
      for (int k = 0; k < NumberOfNeighbors; k++) {
        if ((indices(k, s) > s) && (dists(k, s) < dist) && (indices(k, s) < ParticleBoundary)) {
          i_pair[niac] = s;
          j_pair[niac] = indices(k, s);
          rij2[niac] = dists(k, s);
          niac++;
        }
      }
    }

    delete nns;
    return niac;
  }
}
