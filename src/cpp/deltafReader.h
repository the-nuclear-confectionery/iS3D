#ifndef DELTAFREADER_H
#define DELTAFREADER_H

#include "ParameterReader.h"
#include "readindata.h"
#include "gaussThermal.h"
#include <fstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <vector> // Add this for std::vector

using namespace std;

class Deltaf_Data
{
    private:
        ParameterReader * paraRdr;

        int hrg_eos; // type of pdg file for hadron resonance gas EoS
        int mode; // type of freezeout surface, VH or VAH
        int df_mode; // type of delta-f correction (e.g., 14-moment, CE, or modified distribution)
        int include_baryon;

        string hrg_eos_path;
        string urqmd = "deltaf_coefficients/vh/urqmd/"; // directories of df coefficient tables
        string smash = "deltaf_coefficients/vh/smash/";
        string smash_box = "deltaf_coefficients/vh/smash_box/";

    public:
        int points_T;
        int points_muB;

        double T_min;
        double muB_min;

        double dT;
        double dmuB;

        std::vector<double> T_array; // Changed from raw pointer to std::vector
        std::vector<double> muB_array; // Changed from raw pointer to std::vector

        // Coefficients of 14 moment approximation (vhydro)
        std::vector<std::vector<double>> c0_data; // Changed from double**
        std::vector<std::vector<double>> c1_data;
        std::vector<std::vector<double>> c2_data;
        std::vector<std::vector<double>> c3_data;
        std::vector<std::vector<double>> c4_data;

        // Coefficients of Chapman-Enskog expansion (vhydro)
        std::vector<std::vector<double>> F_data; // Changed from double**
        std::vector<std::vector<double>> G_data;
        std::vector<std::vector<double>> betabulk_data;
        std::vector<std::vector<double>> betaV_data;
        std::vector<std::vector<double>> betapi_data;

        // Cubic splines of the coefficients
        gsl_spline * c0_spline;
        gsl_spline * c2_spline;
        gsl_spline * c3_spline;

        gsl_spline * F_spline;
        gsl_spline * betabulk_spline;
        gsl_spline * betaV_spline;
        gsl_spline * betapi_spline;

        // Jonah coefficients
        const int jonah_points = 301;       // # lambda interpolation points
        const double lambda_min = -1.0;     // lambda min / max values
        const double lambda_max = 2.0;
        const double delta_lambda = (lambda_max - lambda_min) / ((double)jonah_points - 1.0);

        double * lambda_squared_array;      // squared isotropic momentum scale
        double * z_array;                   // renormalization factor (apart from detLambda)
        double * bulkPi_over_Peq_array;     // bulk pressure output
        double bulkPi_over_Peq_max;         // the maximum bulk pressure in the array

        gsl_spline * lambda_squared_spline; // cubic splines for lambda^2(bulkPi/Peq) and z(bulkPi/Peq)
        gsl_spline * z_spline;

        Deltaf_Data(ParameterReader * paraRdr_in);
        ~Deltaf_Data();

        void load_df_coefficient_data();    // Modified function
        void construct_cubic_splines();

        void compute_jonah_coefficients(particle_info * particle_data, int Nparticle);

        deltaf_coefficients evaluate_df_coefficients(double T, double muB, double E, double P, double bulkPi, int icell);

        deltaf_coefficients cubic_spline(double T, double E, double P, double bulkPi);

        double calculate_bilinear(std::vector<std::vector<double>>& f_data, double T, double muB, double TL, double TR, double muBL, double muBR, int iTL, int iTR, int imuBL, int imuBR);

        deltaf_coefficients bilinear_interpolation(double T, double muB, double E, double P, double bulkPi, int icell);

        void test_df_coefficients(double bulkPi_over_P);

        void compute_particle_densities(particle_info * particle_data, int Nparticle);
};

#endif
