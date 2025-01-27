
#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
#include<cmath>
#include<vector>
#include<sys/time.h>

#include "iS3D.h"
#include "Table.h"
#include "readindata.h"
#include "emissionfunction.h"
#include "arsenal.h"
#include "ParameterReader.h"
#include "deltafReader.h"

IS3D::IS3D()
{
}

IS3D::~IS3D()
{
}

void IS3D::read_fo_surf_from_memory(
                              std::vector<double> tau_in,
                              std::vector<double> x_in,
                              std::vector<double> y_in,
                              std::vector<double> eta_in,
                              std::vector<double> dsigma_tau_in,
                              std::vector<double> dsigma_x_in,
                              std::vector<double> dsigma_y_in,
                              std::vector<double> dsigma_eta_in,
                              std::vector<double> E_in,
                              std::vector<double> T_in,
                              std::vector<double> P_in,
                              std::vector<double> ux_in,
                              std::vector<double> uy_in,
                              std::vector<double> un_in,
                              std::vector<double> pixx_in,
                              std::vector<double> pixy_in,
                              std::vector<double> pixn_in,
                              std::vector<double> piyy_in,
                              std::vector<double> piyn_in,
                              std::vector<double> pinn_in,
                              std::vector<double> Pi_in
                              )
{
  tau = tau_in;
  x = x_in;
  y = y_in;
  eta = eta_in;
  dsigma_tau = dsigma_tau_in;
  dsigma_x = dsigma_x_in;
  dsigma_y = dsigma_y_in;
  dsigma_eta = dsigma_eta_in;
  E = E_in;
  T = T_in;
  P = P_in;
  ux = ux_in;
  uy = uy_in;
  un = un_in;
  pixx = pixx_in;
  pixy = pixy_in;
  pixn = pixn_in;
  piyy = piyy_in;
  piyn = piyn_in;
  pinn = pinn_in;
  Pi = Pi_in;
}

void IS3D::run_particlization(int fo_from_file)
{
  cout << "Welcome to iS3D, a program to accelerate particle spectra computation from 3+1D Hydro Freezeout Surfaces!" << endl;
  cout << "Derek Everett, Mike McNelis, Sameed Pervaiz and Lipei Du (2018)" << endl;
  cout << "Based on iSpectra v1.2 : Chun Shen and Zhi Qiu" << endl;

  // Read-in parameters
  cout << "Reading in parameters:\n" << endl;
  ParameterReader *paraRdr = new ParameterReader;
  paraRdr->readFromFile("iS3D_parameters.dat");
  paraRdr->echo();

  string pathToInput = "input";

  // load freeze out information
  long FO_length = 0;
  FO_data_reader freeze_out_data(paraRdr, pathToInput);
  std::cout << "Reading in freezeout surface from file" << std::endl;
  if (fo_from_file == 1) FO_length = freeze_out_data.get_number_cells(); //load length of surface from file
  else FO_length = tau.size(); //load length of surface from memory
  std::cout << "Number of freezeout cells: " << FO_length << std::endl;
  FO_surf* surf_ptr = new FO_surf[FO_length];

  //load freezeout info from 'input/surface.dat'
  if (fo_from_file) freeze_out_data.read_surf_switch(FO_length, surf_ptr);
  //copy freezeout info from vectors to FO_surf
  else
  {
    printf("Reading in freezeout surface from memory \n");
    for (long is = 0; is < FO_length; is++)
    {
      //CHECK UNITS BETWEEN MODULES!

      //contravariant position
      surf_ptr[is].tau = tau[is];
      surf_ptr[is].x = x[is];
      surf_ptr[is].y = y[is];
      surf_ptr[is].eta = eta[is];
      //covariant surface normal vector
      surf_ptr[is].dat = dsigma_tau[is];
      surf_ptr[is].dax = dsigma_x[is];
      surf_ptr[is].day = dsigma_y[is];
      surf_ptr[is].dan = dsigma_eta[is];
      //contravariant flow velocity
      surf_ptr[is].ux = ux[is];
      surf_ptr[is].uy = uy[is];
      surf_ptr[is].un = un[is];

      surf_ptr[is].E = E[is]; //energy density in units [?]
      surf_ptr[is].T = T[is]; //Temperature in units [?]
      surf_ptr[is].P = P[is]; //Pressure in units [?]
      //contravariant shear stress pi^{\mu\nu} in units [?]
      surf_ptr[is].pixx = pixx[is];
      surf_ptr[is].pixy = pixy[is];
      surf_ptr[is].pixn = pixn[is];
      surf_ptr[is].piyy = piyy[is];
      surf_ptr[is].piyn = piyn[is];
      //bulk pressure in units [?]
      surf_ptr[is].bulkPi = Pi[is];
    }
  }


  // load particle info from one of the pdg files (depends on hrg_eos parameter)
  particle_info *particle_data = new particle_info[Maxparticle];
  PDG_Data pdg(paraRdr);
  int Nparticle = pdg.read_resonances(particle_data);


  // df coefficient data
  Deltaf_Data * df_data = new Deltaf_Data(paraRdr);
  df_data->load_df_coefficient_data();
  df_data->construct_cubic_splines();
  df_data->compute_jonah_coefficients(particle_data, Nparticle);
  df_data->compute_particle_densities(particle_data, Nparticle);
  df_data->test_df_coefficients(-0.1);
    


  //FOR THIS READ IN TO WORK PROPERLY, chosen_particles.dat MUST HAVE AN EMPTY ROW AT THE END!
  //switch to different method of reading chosen_particles.dat file that doesn't
  //have this undesirable feature
  Table chosen_particles("PDG/chosen_particles.dat"); // skip others except for these particles
  
  cout << "Total number of freezeout cells: " <<  FO_length << endl;
  cout << "Number of chosen particles: " << chosen_particles.getNumberOfRows() << endl;

  Table pT_tab("tables/pT_gauss_legendre_table.dat");   // pT value and weight table
  Table phi_tab("tables/phi_gauss_legendre_table.dat"); // phi value and weight table
  Table y_tab("tables/y_trapezoid_table_21pt.dat"); // y values and weights
  string etaTableFile = "tables/eta/eta_trapezoid_table_241pt.dat"; // for smooth C.F.
  int operation = paraRdr->getVal("operation");
  if (operation == 2) etaTableFile = "tables/eta/eta_trapezoid_table_41pt.dat";
  Table eta_tab(etaTableFile); //eta_s values and weights
  std::cout << "Emitting particles..." << std::endl;
  EmissionFunctionArray efa(paraRdr, &chosen_particles, &pT_tab, &phi_tab, &y_tab, &eta_tab, particle_data, Nparticle, surf_ptr, FO_length, df_data);
  std::cout << "Finished emitting particles." << std::endl;
  // old version for jetscape 
  //std::vector<Sampled_Particle> particle_event_list_in;

  //for oversampling in jetscape
  std::vector< std::vector< Sampled_Particle > >  particle_event_list_in;
  std::cout << "Calculating spectra..." << std::endl;
  efa.calculate_spectra(particle_event_list_in);
  std::cout << "Finished calculating spectra." << std::endl;

  //copy final particle list to memory to pass to JETSCAPE module
  if (operation == 2)
  {
    cout << "Copying final event particle lists to memory" << endl;
    cout << "Event particle list contains " << particle_event_list_in.size() << " events" << endl;
    final_particles_ = particle_event_list_in;
  }
  std::cout << "Finished copying final event particle lists to memory." << std::endl;
  delete [] surf_ptr;
  std::cout << "Finished deleting surf_ptr." << std::endl;
  delete paraRdr;
  std::cout << "Finished deleting paraRdr." << std::endl;
  delete df_data;
  std::cout << "Finished deleting df_data." << std::endl;

  cout << "Done calculating particle spectra. Output stored in results folder. Goodbye!" << endl;

}
