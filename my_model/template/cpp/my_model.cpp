// Copyright (C) 2022 National Center for Atmospheric Research,
// National Technology & Engineering Solutions of Sandia, LLC (NTESS),
// and the U.S. Environmental Protection Agency (USEPA)
//
// SPDX-License-Identifier: Apache-2.0
//
#include "my_model.hpp"

#include <aero/aero.hpp>
#include <cstring>

#ifdef AERO_USE_NETCDF
#include <netcdf.h>

// NetCDF error code handler
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#endif

// Aerosol state specific to this model
class MyState : public aero::State {
public:
  // (All data in the state is public for simplicity.)
  std::vector<aero::Real> od,      // aerosol optical depth [m]
                          od_ssa,  // aerosol single scattering albedo [-]
                          od_asym, // aerosol asymmetric scattering optical
                                   // depth [m]
                          od_work; // working array for optical depths

  // Constructor
  MyState()
    : aero::State() {

    // We initialize vectors with data pulled from
    // https://acp.copernicus.org/articles/18/7815/2018/acp-18-7815-2018-f03.pdf
    // For this state, we use the "mixed type" measures, indicated by triangles.
    // Note the units in the fields of the my_model_data_t struct at the top of
    // this file. The data here corresponds to the grid interfaces shown in the
    // figure.
    // The data arrays are ordered so they correspond to the following
    // wavelengths, expressed in descending order:
    // {1020.0, 870.0, 675.0, 440.0} [nm]
    // This corresponds to a grid with interfaces expressed in ascending
    // wave numbers [m-1].
    od      = {0.27, 0.35, 0.5, 0.75};       // top left
    od_ssa  = {0.88, 0.895, 0.905, 0.88};    // middle left
    od_asym = {0.3, 0.035, 0.045, 0.09};     // top right
    od_work = std::vector<aero::Real>(4, 0.0);
  }
};

MyModel::MyModel(const char* description_file)
  : aero::Model(),
    grid_(create_grid_()) {

  // Initialize the aerosol grid with wavelength data pulled from
  // https://acp.copernicus.org/articles/18/7815/2018/acp-18-7815-2018-f03.pdf
  // We specify wavelengths in descending order so their wave numbers appear in
  // ascending order in the grid interfaces array.
  aero::Real wavelengths[] = {1020.0, 870.0, 675.0, 440.0}; // [nm]

#ifdef AERO_USE_NETCDF
  // read some NetCDF data
  if (strlen(description_file)>0) {
    int ncid, retval;;
    if(retval = nc_open(description_file, NC_NOWRITE, &ncid)) ERR(retval);
    if(retval = nc_close(ncid)) ERR(retval);
  }
#endif

  // Convert to wave numbers for the grid's interfaces.
  std::vector<aero::Real> wave_numbers;
  for (size_t i = 0; i < 4; ++i) {
    wave_numbers.push_back(1e-9 / wavelengths[i]); // [m-1]
  }
}

MyModel::~MyModel() {
  delete grid_;
}

// This helper method creates the grid for optical properties.
aero::Grid* MyModel::create_grid_() {
  // Initialize the aerosol grid with wavelength data pulled from
  // https://acp.copernicus.org/articles/18/7815/2018/acp-18-7815-2018-f03.pdf
  aero::Real wavelengths[] = {440.0, 675.0, 870.0, 1020.0}; // [nm]

  // Convert to wave numbers for the grid's interfaces.
  std::vector<aero::Real> wave_numbers;
  for (size_t i = 0; i < 4; ++i) {
    wave_numbers.push_back(1e-9 / wavelengths[i]); // [m-1]
  }

  // Create an interfaces array and, from it, a grid.
  aero::Array *interfaces = new aero::Array(wave_numbers);
  return new aero::Grid(interfaces);
}

std::string MyModel::name() const {
  return "my model";
}

aero::State* MyModel::create_state() const {
  return new MyState();
}

aero::Grid* MyModel::optics_grid() const {
  return new aero::Grid(grid_->interfaces().clone());
}

void MyModel::compute_optics(aero::State& state,
                             aero::Array& od,
                             aero::Array& od_ssa,
                             aero::Array& od_asym) const {
  // We simply copy the state's optics data into place.
  MyState& my_state = dynamic_cast<MyState&>(state);
  for (size_t i=0; i<my_state.od_work.size(); ++i)
    my_state.od_work[i] = my_state.od[i];
  od.copy_in(my_state.od_work);
  for (size_t i=0; i<my_state.od_work.size(); ++i)
    my_state.od_work[i] *= my_state.od_ssa[i];
  od_ssa.copy_in(my_state.od_work);
  for (size_t i=0; i<my_state.od_work.size(); ++i)
    my_state.od_work[i] *= my_state.od_asym[i];
  od_asym.copy_in(my_state.od_work);
}
