
#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

#include "iS3D.h"
#include "deltafReader.h"
#include "ParameterReader.h"
#include "readindata.h"

using namespace std;


Deltaf_Data::Deltaf_Data(ParameterReader * paraRdr_in)
{
  paraRdr = paraRdr_in;

  hrg_eos = paraRdr->getVal("hrg_eos");
  mode = paraRdr->getVal("mode");
  df_mode = paraRdr->getVal("df_mode");
  include_baryon = paraRdr->getVal("include_baryon");

  if(hrg_eos == 1)
  {
    hrg_eos_path = urqmd;
  }
  else if(hrg_eos == 2)
  {
    hrg_eos_path = smash;
  }
  else if(hrg_eos == 3)
  {
    hrg_eos_path = smash_box;
  }
  else
  {
    printf("Error: please choose hrg_eos = (1,2,3)\n");
    exit(-1);
  }
}

Deltaf_Data::~Deltaf_Data()
{
  // is there any harm in deallocating memory, while it's being used?
  //gsl_spline_free(c0_spline);
  //gsl_spline_free(c2_spline);
  //gsl_spline_free(c3_spline);
//
  //gsl_spline_free(F_spline);
  //gsl_spline_free(betabulk_spline);
  //gsl_spline_free(betaV_spline);
  //gsl_spline_free(betapi_spline);
//
  //gsl_spline_free(lambda_squared_spline);
  //gsl_spline_free(z_spline);
}

void Deltaf_Data::load_df_coefficient_data()
{
    std::cout << "Reading in 14 moment and Chapman-Enskog coefficient tables..." << std::endl;

    // Define file paths
    std::vector<std::string> file_names = {"c0.dat", "c1.dat", "c2.dat", "c3.dat", "c4.dat",
                                           "F.dat", "G.dat", "betabulk.dat", "betaV.dat", "betapi.dat"};
    std::vector<std::ifstream> files;

    // Open files and check for errors
    for (const auto& name : file_names) {
        std::string full_path = hrg_eos_path + name;
        files.emplace_back(full_path);
        if (!files.back().is_open()) {
            throw std::runtime_error("Couldn't open coefficient file: " + full_path);
        }
    }

    // Read dimensions from the first line of each file
    files[0] >> points_T >> points_muB;
    for (size_t i = 1; i < files.size(); ++i) {
        int temp_T, temp_muB;
        files[i] >> temp_T >> temp_muB;
        if (temp_T != points_T || temp_muB != points_muB) {
            throw std::runtime_error("Mismatch in T and muB dimensions across files.");
        }
    }

    if (!include_baryon) {
        points_muB = 1;
    }

    // Skip headers in all files
    for (auto& file : files) {
        std::string header;
        std::getline(file, header);
        std::getline(file, header);
    }

    // Initialize data storage
    T_array.resize(points_T);
    muB_array.resize(points_muB);

    std::vector<std::vector<std::vector<double>>> data(10, std::vector<std::vector<double>>(points_muB, std::vector<double>(points_T)));

    // Read data into arrays
    for (int iB = 0; iB < points_muB; ++iB) {
        for (int iT = 0; iT < points_T; ++iT) {
            for (size_t i = 0; i < files.size(); ++i) {
                files[i] >> T_array[iT] >> muB_array[iB] >> data[i][iB][iT];
            }
        }
    }

    // Assign data to respective members
    c0_data = std::move(data[0]);
    c1_data = std::move(data[1]);
    c2_data = std::move(data[2]);
    c3_data = std::move(data[3]);
    c4_data = std::move(data[4]);
    F_data = std::move(data[5]);
    G_data = std::move(data[6]);
    betabulk_data = std::move(data[7]);
    betaV_data = std::move(data[8]);
    betapi_data = std::move(data[9]);

    // Compute grid properties
    T_min = T_array.front();
    muB_min = muB_array.front();
    dT = std::fabs(T_array[1] - T_array[0]);
    dmuB = std::fabs(muB_array[1] - muB_array[0]);

    std::cout << "done" << std::endl;
}


void Deltaf_Data::compute_jonah_coefficients(particle_info * particle_data, int Nparticle)
{
  // allocate memory for the arrays
  lambda_squared_array = (double *)calloc(jonah_points, sizeof(double));
  z_array = (double *)calloc(jonah_points, sizeof(double));
  bulkPi_over_Peq_array = (double *)calloc(jonah_points, sizeof(double));

  bulkPi_over_Peq_max = -1.0;      // default to lowest value

  // get the average temperature, energy density, pressure
  Plasma QGP;
  QGP.load_thermodynamic_averages();

  const double T = QGP.temperature;    // GeV (assumes freezeout surface of constant temperature)

  // gauss laguerre roots and weights
  Gauss_Laguerre gla;
  gla.load_roots_and_weights("tables/gla_roots_weights_32_points.txt");

  const int pbar_pts = gla.points;

  double * pbar_root1 = gla.root[1];
  double * pbar_root2 = gla.root[2];

  double * pbar_weight1 = gla.weight[1];
  double * pbar_weight2 = gla.weight[2];

  // calculate the interpolation points of z(bulkPi/P), lambda(bulkPi/P)
  for(int i = 0; i < jonah_points; i++)
  {
    double lambda = lambda_min + (double)i * delta_lambda;

    double E = 0.0;                       // energy density (computed with kinetic theory)
    double P = 0.0;                       // pressure
    double E_mod = 0.0;                   // modified energy density
    double P_mod = 0.0;                   // modified pressure

    // calculate modified energy density (sum over hadron resonance contributions)
    for(int n = 0; n < Nparticle; n++)
    {
      double degeneracy = (double)particle_data[n].gspin;
      double mass = particle_data[n].mass;
      double sign = (double)particle_data[n].sign;

      double mbar = mass / T;

      if(mass == 0.0) continue;   // I skip the photon (Gamma) because the calculation breaks down for lambda = -1.0

      // ignore common prefactor = pow(T,4) / two_pi2_hbarC3 (since they will cancel out)
      E += degeneracy * Gauss1D_mod(E_mod_int, pbar_root2, pbar_weight2, pbar_pts, mbar, 0.0, sign);
      P += (1.0 / 3.0) * degeneracy * Gauss1D_mod(P_mod_int, pbar_root2, pbar_weight2, pbar_pts, mbar, 0.0, sign);

      E_mod += degeneracy * Gauss1D_mod(E_mod_int, pbar_root2, pbar_weight2, pbar_pts, mbar, lambda, sign);
      P_mod += (1.0 / 3.0) * degeneracy * Gauss1D_mod(P_mod_int, pbar_root2, pbar_weight2, pbar_pts, mbar, lambda, sign);
    }

    // jonah's formula (ignoring detLambda factor i.e. n_mod / n -> 1)
    double z = E / E_mod;
    double bulkPi_over_Peq = (P_mod / P) * z  -  1.0;

    // set the arrays and update the max bulk pressure
    lambda_squared_array[i] = lambda * lambda;
    z_array[i] = z;
    bulkPi_over_Peq_array[i] = bulkPi_over_Peq;
    bulkPi_over_Peq_max = max(bulkPi_over_Peq_max, bulkPi_over_Peq);

    //cout << lambda_squared_array[i] << "\t" << z_array[i] << "\t" << bulkPi_over_Peq_array[i] << endl;
  }

  // now construct cubic splines for lambda(bulkPi/Peq) and z(bulkPi/Peq)
  //lambda_squared_spline = gsl_spline_alloc(gsl_interp_cspline, jonah_points);
  //z_spline = gsl_spline_alloc(gsl_interp_cspline, jonah_points);
//
  //gsl_spline_init(lambda_squared_spline, bulkPi_over_Peq_array, lambda_squared_array, jonah_points);
  //gsl_spline_init(z_spline, bulkPi_over_Peq_array, z_array, jonah_points);
}


void Deltaf_Data::construct_cubic_splines()
{
  // Allocate memory for cubic splines
  c0_spline = gsl_spline_alloc(gsl_interp_cspline, points_T);
  c2_spline = gsl_spline_alloc(gsl_interp_cspline, points_T);
  c3_spline = gsl_spline_alloc(gsl_interp_cspline, points_T);

  F_spline = gsl_spline_alloc(gsl_interp_cspline, points_T);
  betabulk_spline = gsl_spline_alloc(gsl_interp_cspline, points_T);
  betapi_spline = gsl_spline_alloc(gsl_interp_cspline, points_T);
  betaV_spline = gsl_spline_alloc(gsl_interp_cspline, points_T);

  // Initialize the cubic splines
  //gsl_spline_init(c0_spline, T_array, c0_data[0], points_T);
  //gsl_spline_init(c2_spline, T_array, c2_data[0], points_T);
  //gsl_spline_init(c3_spline, T_array, c3_data[0], points_T);
//
  //gsl_spline_init(F_spline, T_array, F_data[0], points_T);
  //gsl_spline_init(betabulk_spline, T_array, betabulk_data[0], points_T);
  //gsl_spline_init(betaV_spline, T_array, betaV_data[0], points_T);
  //gsl_spline_init(betapi_spline, T_array, betapi_data[0], points_T);

}


deltaf_coefficients Deltaf_Data::cubic_spline(double T, double E, double P, double bulkPi)
{
  deltaf_coefficients df;

  gsl_interp_accel * accel_T = gsl_interp_accel_alloc();    // for temperature dependent functions
  gsl_interp_accel * accel_bulk = gsl_interp_accel_alloc(); // for bulkPi/Peq dependent functions

  switch(df_mode)
  {
    case 1: // 14 moment
    {
      // undo the temperature power scaling of coefficients
      double T4 = T * T * T * T;

      df.c0 = gsl_spline_eval(c0_spline, T, accel_T) / T4;
      df.c1 = 0.0;
      df.c2 = gsl_spline_eval(c2_spline, T, accel_T) / T4;
      df.c3 = 0.0;
      df.c4 = 0.0;
      df.shear14_coeff = 2.0 * T * T * (E + P);

      break;
    }
    case 2: // Chapman Enskog
    case 3: // Modified (Mike)
    {
      // undo the temperature power scaling of coefficients
      double T4 = T * T * T * T;

      df.F = gsl_spline_eval(F_spline, T, accel_T) * T;
      df.G = 0.0;
      df.betabulk = gsl_spline_eval(betabulk_spline, T, accel_T) * T4;
      df.betaV = 1.0;
      df.betapi = gsl_spline_eval(betapi_spline, T, accel_T) * T4;

      break;
    }
    case 4: // Modified (Jonah)
    {
      // undo the temperature power scaling of betapi
      double T4 = T * T * T * T;

      double lambda_squared = gsl_spline_eval(lambda_squared_spline, (bulkPi / P), accel_bulk);
      if(bulkPi < 0.0)
      {
        df.lambda = - sqrt(lambda_squared);
      }
      else if(bulkPi > 0.0)
      {
        df.lambda = sqrt(lambda_squared);
      }
      df.z = gsl_spline_eval(z_spline, (bulkPi / P), accel_bulk);
      df.betapi = gsl_spline_eval(betapi_spline, T, accel_T) * T4;
      // linearized correction to lambda, z
      df.delta_lambda = bulkPi / (5.0 * df.betapi -  3.0 * P * (E + P) / E);
      df.delta_z = - 3.0 * df.delta_lambda * P / E;

      break;
    }
    case 8:
    {
      df.shear14_coeff = 2.0 * T * T * (E + P);
      std::cout << "Using CCAKE  df correction" << std::endl;
    }
    default:
    {
      printf("Error: choose df_mode = (1,2,3,4)\n");
      exit(-1);
    }
  }

  gsl_interp_accel_free(accel_T);
  gsl_interp_accel_free(accel_bulk);

  return df;
}

double Deltaf_Data::calculate_bilinear(std::vector<std::vector<double>>& f_data, double T, double muB, double TL, double TR, double muBL, double muBR, int iTL, int iTR, int imuBL, int imuBR)
{
    // Bilinear interpolation for f(muB, T)
    // Diagram of points:
    //  f_LL    f_RL (T-axis: Left to Right)
    //
    //  f_LR    f_RR
    // (muB-axis: Bottom to Top)

    //std::cout << "f starting" << std::endl;
    //std::cout << "iTL = " << iTL << ", iTR = " << iTR << ", imuBL = " << imuBL << ", imuBR = " << imuBR << std::endl;
//
    // Ensure indices are within bounds
    if (iTL >= points_T || iTR >= points_T || imuBL >= points_muB || imuBR >= points_muB || iTL < 0 || iTR < 0 || imuBL < 0 || imuBR < 0) {
        throw std::out_of_range("Index out of range in calculate_bilinear.");
    }

    // Access the required points
    double f_LL = f_data[imuBL][iTL]; // Bottom-left
    double f_RL = f_data[imuBL][iTR]; // Bottom-right
    double f_LR = f_data[imuBR][iTL]; // Top-left
    double f_RR = f_data[imuBR][iTR]; // Top-right

    // Perform bilinear interpolation
    double result = ((f_LL * (TR - T) + f_RL * (T - TL)) * (muBR - muB) +
                     (f_LR * (TR - T) + f_RR * (T - TL)) * (muB - muBL)) /
                    (dT * dmuB);

    return result;
}


deltaf_coefficients Deltaf_Data::bilinear_interpolation(double T, double muB, double E, double P, double bulkPi, int icell)
{
  // left and right T, muB indices
  int iTL = (int)floor((T - T_min) / dT);
  int iTR = iTL + 1;
  //std::cout << muB << " " << muB_min << " " << dmuB << std::endl;
  int imuBL = (int)floor((muB - muB_min) / dmuB);
  //std::cout << muB << " " << muB_min << " " << dmuB << " " << imuBL << std::endl;
  int imuBR = imuBL + 1;
  double TL, TR, muBL, muBR;

  if(!(iTL >= 0 && iTR < points_T) || !(imuBL >= 0 && imuBR < points_muB))
  {
    printf("Error: (T,muB) outside df coefficient table. Exiting...\n");
    //exit(-1);
  }
  else
  {
    TL = T_array[iTL];
    TR = T_array[iTR];
    muBL = muB_array[imuBL];
    muBR = muB_array[imuBR];
  }

  deltaf_coefficients df;
  switch(df_mode)
  {
    case 1:
    {
      double T3 = T * T * T;
      double T4 = T3 * T;
      double T5 = T4 * T;

      // bilinear interpolated values & undo temperature power scaling
      df.c0 = calculate_bilinear(c0_data, T, muB, TL, TR, muBL, muBR, iTL, iTR, imuBL, imuBR) / T4;
      df.c1 = calculate_bilinear(c1_data, T, muB, TL, TR, muBL, muBR, iTL, iTR, imuBL, imuBR) / T3;
      df.c2 = calculate_bilinear(c2_data, T, muB, TL, TR, muBL, muBR, iTL, iTR, imuBL, imuBR) / T4;
      df.c3 = calculate_bilinear(c3_data, T, muB, TL, TR, muBL, muBR, iTL, iTR, imuBL, imuBR) / T4;
      df.c4 = calculate_bilinear(c4_data, T, muB, TL, TR, muBL, muBR, iTL, iTR, imuBL, imuBR) / T5;
      df.shear14_coeff = 2.0 * T * T * (E + P);

      break;
    }
    case 2:
    case 3:
    {
      double T3 = T * T * T;
      double T4 = T3 * T;

      // bilinear interpolated values & undo temperature power scaling
      df.F = calculate_bilinear(F_data, T, muB, TL, TR, muBL, muBR, iTL, iTR, imuBL, imuBR) * T;
      df.G = calculate_bilinear(G_data, T, muB, TL, TR, muBL, muBR, iTL, iTR, imuBL, imuBR);
      df.betabulk = calculate_bilinear(betabulk_data, T, muB, TL, TR, muBL, muBR, iTL, iTR, imuBL, imuBR) * T4;
      df.betaV = calculate_bilinear(betaV_data, T, muB, TL, TR, muBL, muBR, iTL, iTR, imuBL, imuBR) * T3;
      df.betapi = calculate_bilinear(betapi_data, T, muB, TL, TR, muBL, muBR, iTL, iTR, imuBL, imuBR) * T4;

      break;
    }
    case 4:
    {
      printf("Bilinear interpolation error: Jonah df doesn't work for nonzero muB. Exiting..\n");
      exit(-1);
    }
    case 8:
    {
      df.shear14_coeff = 2.0 * T * T * (E + P);
      break;
    }
    default:
    {
      printf("Bilinear interpolation error: choose df_mode = (1,2,3). Exiting..\n");
      exit(-1);
    }
  }

  return df;
}

deltaf_coefficients Deltaf_Data::evaluate_df_coefficients(double T, double muB, double E, double P, double bulkPi, int icell)
{
  // evaluate the df coefficients by interpolating the data

  deltaf_coefficients df;

  if(!include_baryon)
  {
    df = cubic_spline(T, E, P, bulkPi);   // cubic spline interpolation wrt T (at muB = 0)
  }
  else
  {
    // muB on freezeout surface should be nonzero in general
    // otherwise should set include_baryon = 0
    if(T== 0.187394)
    {
      std::cout << "T = " << T << ", muB = " << muB << ", E = " << E << ", P = " << P << ", bulkPi = " << bulkPi << std::endl;
      int icell = 1;
    } 
    df = bilinear_interpolation(T, muB, E, P, bulkPi, icell);  // bilinear wrt (T, muB)
  }

  return df;
}


void Deltaf_Data::test_df_coefficients(double bulkPi_over_P)
{
  // test print the output of the df coefficients at average temperature, etc (and a fixed value of bulkPi)

  Plasma QGP;
  QGP.load_thermodynamic_averages();

  double E = QGP.energy_density;
  double T = QGP.temperature;
  double P = QGP.pressure;
  double muB = QGP.baryon_chemical_potential;
  double bulkPi = bulkPi_over_P * P;

  deltaf_coefficients df = evaluate_df_coefficients(T, muB, E, P, bulkPi, 0);

  if(df_mode == 1)
  {
    printf("\n(c0, c1, c2, c3, c4, shear14) = (%lf, %lf, %lf, %lf, %lf, %lf)\n\n", df.c0, df.c1, df.c2, df.c3, df.c4, df.shear14_coeff);
  }
  else if(df_mode == 2 || df_mode == 3)
  {
    printf("\n(F, G, betabulk, betaV, betapi) = (%lf, %lf, %lf, %lf, %lf)\n\n", df.F, df.G, df.betabulk, df.betaV, df.betapi);
  }
  else if(df_mode == 4)
  {
    printf("\n(lambda, z, dlambda, dz, betapi) = (%lf, %lf, %lf, %lf, %lf)\n\n", df.lambda, df.z, df.delta_lambda, df.delta_z, df.betapi);
  }
}

void Deltaf_Data::compute_particle_densities(particle_info * particle_data, int Nparticle)
{
  // get the average temperature, energy density, pressure, etc.
  Plasma QGP;
  QGP.load_thermodynamic_averages();

  const double T = QGP.temperature;    // GeV
  const double E = QGP.energy_density; // GeV / fm^3
  const double P = QGP.pressure;       // GeV / fm^3
  const double muB = QGP.baryon_chemical_potential;
  const double nB = QGP.net_baryon_density;
  const double muS = QGP.strange_chemical_potential;
  const double nS = QGP.net_strange_density;
  const double muQ = QGP.charge_chemical_potential;
  const double nQ = QGP.net_charge_density;
  

  deltaf_coefficients df = evaluate_df_coefficients(T, muB, E, P, 0.0, 0);

  double alphaB = muB / T;
  double alphaQ = muQ / T;
  double alphaS = muS / T;
  double baryon_enthalpy_ratio = nB / (E + P);

  // gauss laguerre roots and weights
  Gauss_Laguerre gla;
  gla.load_roots_and_weights("tables/gla_roots_weights_32_points.txt");

  const int pbar_pts = gla.points;

  double * pbar_root1 = gla.root[1];
  double * pbar_root2 = gla.root[2];
  double * pbar_root3 = gla.root[3];

  double * pbar_weight1 = gla.weight[1];
  double * pbar_weight2 = gla.weight[2];
  double * pbar_weight3 = gla.weight[3];


  // calculate the equilibrium densities and the
  // bulk / diffusion corrections of each particle
  for(int i = 0; i < Nparticle; i++)
  {
    double mass = particle_data[i].mass;
    double degeneracy = (double)particle_data[i].gspin;
    double baryon = (double)particle_data[i].baryon;
    double charge = (double)particle_data[i].charge;
    double strange = (double)particle_data[i].strange;
    double sign = (double)particle_data[i].sign;
    double mbar = mass / T;

    // equilibrium density
    double neq_fact = degeneracy * pow(T,3) / two_pi2_hbarC3;
    double neq = neq_fact * GaussThermal(neq_int, pbar_root1, pbar_weight1, pbar_pts, mbar, alphaB, baryon,
                                         alphaQ, charge, alphaS, strange, sign);

    // bulk and diffusion density corrections
    double dn_bulk = 0.0;
    double dn_diff = 0.0;

    switch(df_mode)
    {
      case 1: // 14 moment (not sure what the status is)
      {
        double c0 = df.c0;
        double c1 = df.c1;
        double c2 = df.c2;
        double c3 = df.c3;
        double c4 = df.c4;

        double J10_fact = degeneracy * pow(T,3) / two_pi2_hbarC3;
        double J20_fact = degeneracy * pow(T,4) / two_pi2_hbarC3;
        double J30_fact = degeneracy * pow(T,5) / two_pi2_hbarC3;
        double J31_fact = degeneracy * pow(T,5) / two_pi2_hbarC3 / 3.0;

        double J10 =  J10_fact *GaussThermal(J10_int, pbar_root1, pbar_weight1, pbar_pts, mbar, alphaB, baryon,
                                         alphaQ, charge, alphaS, strange, sign);
        double J20 = J20_fact * GaussThermal(J20_int, pbar_root2, pbar_weight2, pbar_pts, mbar, alphaB, baryon,
                                         alphaQ, charge, alphaS, strange, sign);
        double J30 = J30_fact * GaussThermal(J30_int, pbar_root3, pbar_weight3, pbar_pts, mbar, alphaB, baryon,
                                         alphaQ, charge, alphaS, strange, sign);
        double J31 = J31_fact * GaussThermal(J31_int, pbar_root3, pbar_weight3, pbar_pts, mbar, alphaB, baryon,
                                         alphaQ, charge, alphaS, strange, sign);

        dn_bulk = ((c0 - c2) * mass * mass * J10 +  c1 * baryon * J20  +  (4.0 * c2 - c0) * J30);
        // these coefficients need to be loaded.
        // c3 ~ cV / V
        // c4 ~ 2cW / V
        dn_diff = baryon * c3 * neq * T  +  c4 * J31;   // not sure if this is right...
        break;
      }
      case 2: // Chapman-Enskog
      case 3: // Modified (Mike)
      {
        double F = df.F;
        double G = df.G;
        double betabulk = df.betabulk;
        double betaV = df.betaV;

        double J10_fact = degeneracy * pow(T,3) / two_pi2_hbarC3;
        double J11_fact = degeneracy * pow(T,3) / two_pi2_hbarC3 / 3.0;
        double J20_fact = degeneracy * pow(T,4) / two_pi2_hbarC3;

        double J10 = J10_fact * GaussThermal(J10_int, pbar_root1, pbar_weight1, pbar_pts, mbar, alphaB, baryon,
                                         alphaQ, charge, alphaS, strange, sign);
        double J11 = J11_fact * GaussThermal(J11_int, pbar_root1, pbar_weight1, pbar_pts, mbar, alphaB, baryon,
                                         alphaQ, charge, alphaS, strange, sign);
        double J20 = J20_fact * GaussThermal(J20_int, pbar_root2, pbar_weight2, pbar_pts, mbar, alphaB, baryon,
                                         alphaQ, charge, alphaS, strange, sign);

        dn_bulk = (neq + (baryon * J10 * G) + (J20 * F / pow(T,2))) / betabulk;
        dn_diff = (neq * T * baryon_enthalpy_ratio  -  baryon * J11) / betaV;

        break;
      }
      case 4:
      {
        // bulk/diffusion densities not needed for jonah
        break;
      }
      case 8:
      {
        break;
      }
      default:
      {
        cout << "Please choose df_mode = (1,2,3,4) in parameters.dat" << endl;
        exit(-1);
      }
    } // df_mode

    particle_data[i].equilibrium_density = neq;
    particle_data[i].bulk_density = dn_bulk;
    particle_data[i].diff_density = dn_diff;
  }
}








