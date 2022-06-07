#include <iomanip>
#include <armadillo>

#include "TwoBodyOperator.hh"
#include "TwoBodySpace.hh"
#include "HartreeFock.hh"
#include "ModelSpace.hh"
#include "Operator.hh"
#include "Orbits.hh"

Monopole::~Monopole()
{}

Monopole::Monopole()
{}

Monopole::Monopole(Operator& H)
  : H(&H), modelspace(H.GetModelSpace())
{
  Orbits& orbits = modelspace->GetOrbits();
  int norbs = orbits.GetNumberOrbits();
  for (int i1=0; i1<norbs; i1++) {
    Orbit& o1 = orbits.GetOrbit(i1);
    for (int i3=0; i3<norbs; i3++) {
      Orbit& o3 = orbits.GetOrbit(i3);
      if (o1.k != o3.k) continue;
      for (int i2=0; i2<norbs; i2++) {
        for (int i4=0; i4<norbs; i4++) {
          Orbit& o2 = orbits.GetOrbit(i2);
          Orbit& o4 = orbits.GetOrbit(i4);
          if (o2.k != o4.k) continue;
          if (o1.e2 + o2.e2 != o3.e2 + o4.e2) continue;
          idx_to_ijkl.push_back({i1,i2,i3,i4});
          double norm = 1.0;
          if (i1==i2) norm *= sqrt(2);
          if (i3==i4) norm *= sqrt(2);
          double v = 0.0;
          for (int J=std::abs(o1.j2-o2.j2)/2; J<=(o1.j2+o2.j2)/2; J++){
            if (i1==i2 and J%2==1) continue;
            if (i3==i4 and J%2==1) continue;
            v += (2*J+1) * H.Get2BME_J(J,J,o1,o2,o3,o4);
          }
          v *= norm / (o1.j2+1);
          vmon2.push_back(v);
        }
      }
    }
  }
}

void Monopole::Print() {
  for (int idx=0; idx<idx_to_ijkl.size(); idx++){
    std::cout << "i =" << std::setw(4) << idx_to_ijkl[idx][0]
      << ", j =" << std::setw(4) << idx_to_ijkl[idx][1]
      << ", k =" << std::setw(4) << idx_to_ijkl[idx][2]
      << ", l =" << std::setw(4) << idx_to_ijkl[idx][3]
      << ", Vmon =" << std::setw(12) << vmon2[idx] << std::endl;
  }
}

HartreeFock::~HartreeFock()
{}

HartreeFock::HartreeFock(Operator& H, std::map<int,double> holes)
  : H(H), modelspace(H.GetModelSpace()), holes(holes)
{
  if(not H.orthonormalized){
    std::cout << " OrthoNormalize the Hamiltonian first!" << std::endl;
    exit(0);
  }
  monopole = Monopole(H);
  Orbits orbits = modelspace->GetOrbits();
  int norbs = orbits.GetNumberOrbits();
  C = arma::mat(norbs, norbs, arma::fill::eye);
  rho = arma::mat(norbs, norbs, arma::fill::zeros);
  F = arma::mat(norbs, norbs, arma::fill::zeros);
  V = arma::mat(norbs, norbs, arma::fill::zeros);
  SPEs = arma::vec(norbs, arma::fill::zeros);
  S = H.S;

  r = UpdateFock();
  DiagonalizeFock();
  UpdateDensityMatrix();
  CalcEnergy();
}

void HartreeFock::Solve() {
  for (int n_iter=0; n_iter<100; n_iter++) {
    r = UpdateFock(n_iter);
    DiagonalizeFock();
    UpdateDensityMatrix();
    CalcEnergy();
    PrintStatus(n_iter);
    if (r < 1.e-8) break;
  }
}

double HartreeFock::UpdateFock(int n_itr)
{
  arma::mat Fock_old = F;
  Orbits& orbits = modelspace->GetOrbits();
  int norbs = orbits.GetNumberOrbits();
  V = arma::mat(norbs, norbs, arma::fill::zeros);
  if (n_itr != -1) {
    for (int idx=0; idx<monopole.idx_to_ijkl.size(); idx++) {
      int i = monopole.idx_to_ijkl[idx][0];
      int j = monopole.idx_to_ijkl[idx][1];
      int k = monopole.idx_to_ijkl[idx][2];
      int l = monopole.idx_to_ijkl[idx][3];
      V(i,k) += monopole.vmon2[idx] * rho(j,l);
    }
  }
  F = H.OneBody + V;
  double diff = 0.0;
  for (int i=0; i<norbs; i++) {
    for (int j=0; j<norbs; j++) {
      diff += sqrt( pow(( F(i,j) - Fock_old(i,j) ), 2) );
    }
  }
  return diff;
}

void HartreeFock::DiagonalizeFock() {
  OneBodySpace& obs = modelspace->GetOneBodySpace();
  for (int ich=0; ich< obs.GetNumberChannels(); ich++) {
    std::vector<unsigned long long> tmp(obs.channels[ich].begin(), obs.channels[ich].end());
    arma::uvec sub_idx(std::vector<unsigned long long>(tmp.begin(), tmp.end()));
    arma::mat Fch = F(sub_idx, sub_idx);
    arma::vec eig;
    arma::mat vec;
    arma::eig_sym(eig, vec, Fch);
    SPEs(sub_idx) = eig;
    C(sub_idx, sub_idx) = vec;
    //ElectronStates = GetElectronStates();
  }
}

void HartreeFock::UpdateDensityMatrix() {
  OneBodySpace& obs = modelspace->GetOneBodySpace();
  Orbits& orbits = modelspace->GetOrbits();
  int norbs = orbits.GetNumberOrbits();
  arma::mat tmp(norbs, norbs, arma::fill::zeros);
  std::map<int,int> orbit_idx_to_spe_idx;
  for (int ich=0; ich<obs.GetNumberChannels(); ich++) {
    std::vector<unsigned long long> tmp(obs.channels[ich].begin(), obs.channels[ich].end());
    arma::uvec sub_idx(std::vector<unsigned long long>(tmp.begin(), tmp.end()));

    int n, e2;
    //vec SPEsCh(SPEs.n_cols, fill::zeros);
    int n_ch = obs.channels[ich].size();
    for (int i=0; i<n_ch; i++) {
      int idx_spe = obs.channels[ich][i];
      if (i < n_ch/2 ) {
        n = i;
        e2=-1;
      } else {
        n = i - floor(n_ch/2);
        e2= 1;
      }
      Orbit& o = orbits.GetOrbit(idx_spe);
      int idx = orbits.GetOrbitIndex(n+o.l, o.l, o.j2, e2);
      orbit_idx_to_spe_idx[idx] = idx_spe;
    }
  }
  for ( auto hole_idx : holes) {
    int occ = holes.at(hole_idx.first);
    int spe_idx = orbit_idx_to_spe_idx.at(hole_idx.first);
    tmp(spe_idx, spe_idx) = occ;
  }
  rho = C * tmp * C.t();
}

void HartreeFock::CalcEnergy() {
  Orbits& orbits = modelspace->GetOrbits();
  int norbs = orbits.GetNumberOrbits();
  E1 = 0.0;
  E2 = 0.0;
  for (int i=0; i<norbs; i++) {
    for (int j=0; j<norbs; j++) {
      Orbit& oi = orbits.GetOrbit(i);
      E1 += H.OneBody(i,j) * rho(i,j) * (oi.j2+1);
      E2 += V(i,j) * rho(i,j) * 0.5 * (oi.j2+1);
    }
  }
  EHF = E1+E2;
}

void HartreeFock::PrintStatus(int n_iter) {
  std::cout << "n iter:" << std::setw(4) << n_iter
    << ", OneBody: " << std::fixed << std::setw(12) << E1
    << ", TwoBody: " << std::fixed << std::setw(12) << E2
    << ", EHF: " << std::fixed << std::setw(12) << EHF
    << std::endl;
  //std::cout << "\n";
  //printf("nth iteration: %d, HF energy: %3f ", n_iter, EHF);
  //if (detail) {
  //  std::cout << "\n" << "\n";
  //  F.print("Fock Matrix");
  //  std::cout << "\n";
  //  SPEs.t().print("SPEs");
  //  std::cout << "\n";
  //  C.print("Coeffs");
  //  std::cout << "\n";
  //  rho.print("Density Matrix");
  //  std::cout << "\n";
  //}
}


//
//    HartreeFock::HartreeFock(Operator Ham1, map <int, double> holes1) {
//        // En; // Not sure if this belongs here
//        Ham = Ham1;
//        holes = holes1;
//        modelspace = Ham.modelspace;
//        monopole = Monopole(Ham);
//        monopole.set_monopole2();
//
//        Orbits orbs = modelspace.orbits;
//        int norbs = orbs.get_num_orbits();
//        S = Ham1.S;
//        C.zeros(norbs, norbs);
//        rho.zeros(norbs, norbs);
//        F.zeros(norbs, norbs);
//        V.zeros(norbs, norbs);
//        SPEs.zeros(norbs);
//
//        r = UpdateFock();
//        DiagonalizeFock();
//        UpdateDensityMatrix();
//        CalcEnergy();
//    };
//
//// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//
//    void Monopole::set_monopole2(){
//        Orbits orbs = modelspace.orbits;
//        int norbs = orbs.get_num_orbits();
//        for (int i1=0; i1<norbs; i1++) {
//            Orbit o1 = orbs.get_orbit(i1);
//            for (int i3=0; i3<norbs; i3++) {
//                Orbit o3 = orbs.get_orbit(i3);
//                if (o1.l != o3.k) {continue;}
//                for (int i2=0; i2<norbs; i2++) {
//                    for (int i4=0; i4<norbs; i4++) {
//                        Orbit o2 = orbs.get_orbit(i2);
//                        Orbit o4 = orbs.get_orbit(i4);
//                        if (o2.k != o4.k) {continue;}
//                        if (o1.e + o2.e != o3.e + o4.e) {continue;}
//                        idx_to_ijkl.push_back({i1,i2,i3,i4});
//                        double norm = 1.0;
//                        if (i1==i2) {norm *= sqrt(2);}
//                        if (i3==i4) {norm *= sqrt(2);}
//                        double v = 0.0;
//                        for (int J : boost::irange( floor(abs(o1.j-o2.j)/2),  floor((o1.j+o2.j)/2)+1 ) ) {
//                            if (i1==i2 && J%2==1) {continue;}
//                            if (i3==i4 && J%2==1) {continue;}
//                            v += (2*J+1) * Ham.two.get_2bme_orbitsJ(o1,o2,o3,o4,J,J);
//                        v *= norm / (o1.j+1);
//                        v2.push_back(v);
//                        }
//                    }
//                }
//            }
//        }
//    }
//
//    void Monopole::print_monopole2() {
//        for (int idx=0; idx<idx_to_ijkl.size(); idx++) {
//            printf("%d %d %d %d", idx_to_ijkl[idx][0], idx_to_ijkl[idx][1], idx_to_ijkl[idx][2], idx_to_ijkl[idx][3]);
//            // cout << fixed << setw(2) << idx_to_ijkl[idx] << v2[idx] << endl;
//            // idx_to_ijkl[idx] contains a vector so need to find a more efficient method
//            // need to add v2 but not sure about length of vector
//        }
//    }
//
//
//
//    double HartreeFock::UpdateFock(int n_iter) {
//        Mat<double> Fock_old = F;
//        Orbits orbs = modelspace.orbits;
//        int norbs = orbs.get_num_orbits();
//        Mat<double> V(norbs, norbs, fill::zeros);
//        if (n_iter != -1) {
//            for (int idx=0; idx<monopole.idx_to_ijkl.size(); idx++) {
//                int i = monopole.idx_to_ijkl[idx][0];
//                int j = monopole.idx_to_ijkl[idx][1];
//                int k = monopole.idx_to_ijkl[idx][2];
//                int l = monopole.idx_to_ijkl[idx][3];
//                V(i,k) += monopole.v2[idx] * rho[j,l];
//            }
//            for (int i=0; i<norbs; i++) {
//                for (int j=0; i<i+j; j++) {
//                    V(j,i) = V(i,j);
//                }
//            }
//        }
//        F = Ham.one + V;
//        double r1;
//        for (int i=0; i<norbs; i++) {
//            for (int j=0; j<norbs; j++) {
//                r1 += sqrt( pow(( F(i,j) - Fock_old(i,j) ), 2) );
//            }
//        }
//        return r1;
//    }
//
//    void HartreeFock::DiagonalizeFock() {
//        OneBodySpace one_body_space = modelspace.one;
//        for (int ich=0; ich< one_body_space.number_channels; ich++) {
//            vector <int> filter_idx = one_body_space.channels[ich];
//
//            Mat<double> Fch, Sch;
//            Fch.zeros(size(F));
//            Sch.zeros(size(S));
//
//            for (auto i : filter_idx) {
//                for (auto j : filter_idx) {
//                    Fch(i,j) = F(i,j);
//                    Sch(i,j) = S(i,j);
//                }
//            }
//
//            // Fch.print("Fch");
//            // Sch.print("Sch");
//
//            cx_vec eigval;
//            cx_mat eigvec_col;
//            eig_pair( eigval, eigvec_col, Fch, Sch );
//            cx_mat eigvec_row = eigvec_col.t();
//            // uvec idxs = arma::sort_index(real(eigval), "ascend");
//
//            for (int i=0; i<eigvec_col.n_rows; i++) {
//
//                // int idx = idxs(i);
//                complex<double> norm = arma::as_scalar(eigvec_row.row(i) * Sch * eigvec_col.col(i));
//                eigvec_col.col(i) = -1*eigvec_col.col(i)/norm;
//            }
//
//            for (auto i : filter_idx) {
//                for (auto j : filter_idx) {
//                    C(i,j) = real(eigvec_col(i,j));
//                }
//            }
//
//            for (int i=0; i<filter_idx.size(); i++) { SPEs(i) =  real( eigval(i) ); }
//
//            C.print("C");
//            SPEs.print("SPEs");
//
//        }
//    }
//
//
