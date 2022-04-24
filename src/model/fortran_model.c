#include <aero/model/fortran_model.h>
#include "model_bridge.h"
#include <stdlib.h>

struct aero_model_data_t {
  void fortran_model_;
};

static aero_model_t* aero_fortran_model_create(aero_model_data_t *array_data);

static const char* aero_fortran_model_name(const aero_model_t *model) {
  return aero_bridge_fortran_model_name(model->data_->fortran_model_);
}

static aero_state_t* aero_fortran_model_create_state(const aero_model_t *model) {
  return (aero_state_t*) aero_bridge_fortran_model_create_state(model->data_->fortran_model_);
}

static void aero_fortran_model_free_state(const aero_model_t *model, aero_state_t *state) {
  aero_bridge_fortran_model_free_state(model->data_->fortran_model_, (void*)state);
}

static aero_grid_t* aero_fortran_model_optics_grid(const aero_model_t *model) {
  return (aero_grid_t*) aero_bridge_fortran_model_optics_grid(model->data_->fortran_model_);
}

static void aero_fortran_model_compute_optics(const aero_model_t *model, aero_state_t *state
    aero_array_t *od, aero_array_t *od_ssa, aero_array_t *od_asym) {
  aero_bridge_fortran_model_optics_grid(model->data_->fortran_model_, (void*) state,
      (void*) od, (void*) od_ssa, (void*) od_asym);
}
static void aero_fortran_model_free(aero_model_t *model) {
  aero_bridge_model_free(model->data_->fortran_model_);
}

static aero_model_t* aero_fortran_model_create(aero_model_data_t *model_data) {
  aero_model_t *model;
  model = malloc(sizeof(aero_model_t));
  model->data_ = model_data;
  model->name = aero_fortran_model_name;
  model->create_state = aero_fortran_model_create_state;
  model->free_state = aero_fortran_model_free_state;
  model->optics_grid = aero_fortran_model_optics_grid;
  model->compute_optics = aero_fortran_model_compute_optics;
  model->free = aero_fortran_model_free;
}

aero_model_t* aero_fortran_model_wrap(void *fortran_model) {
  aero_model_data_t *model_data = malloc(sizeof(aero_model_data_t));
  model_data->fortran_model_ = fortran_model;
  return aero_fortran_model_create(model_data);
}

void aero_fortran_model_unwrap(aero_model_t *fortran_model) {
  free(fortran_model->data_);
  free(fortran_model);
}