#include <dlfcn.h>

#include "FWCore/Utilities/interface/Exception.h"
#include "L1Trigger/MLUtilities/interface/HLS4MLModelWrapper.h"

l1t::HLS4MLModelWrapper::HLS4MLModelWrapper() : model_lib_{nullptr}, model_{nullptr}, model_name_{""} {}

l1t::HLS4MLModelWrapper::HLS4MLModelWrapper(std::string const& model_name)
    : model_lib_{nullptr}, model_{nullptr}, model_name_{model_name} {
  load();
}

l1t::HLS4MLModelWrapper::~HLS4MLModelWrapper() { reset(); }

void l1t::HLS4MLModelWrapper::reset() {
  if (model_loaded()) {
    model_.reset();
  }

  if (model_lib_ != nullptr) {
    dlclose(model_lib_);
    model_lib_ = nullptr;
  }

  model_name_ = "";
}

void l1t::HLS4MLModelWrapper::reset(std::string const& model_name) {
  reset();
  model_name_ = model_name;
  load();
}

void l1t::HLS4MLModelWrapper::prepare_input(std::any input) const {
  if (not model_loaded()) {
    throw cms::Exception("InvalidInput") << "failed call to prepare_input: hls4ml model not loaded yet !";
  }

  model_->prepare_input(input);
}

void l1t::HLS4MLModelWrapper::predict() const {
  if (not model_loaded()) {
    throw cms::Exception("InvalidInput") << "failed call to predict: hls4ml model not loaded yet !";
  }

  model_->predict();
}

void l1t::HLS4MLModelWrapper::read_result(std::any result) const {
  if (not model_loaded()) {
    throw cms::Exception("InvalidInput") << "failed call to read_result: hls4ml model not loaded yet !";
  }

  model_->read_result(result);
}

void l1t::HLS4MLModelWrapper::load() {
  if (model_loaded()) {
    throw cms::Exception("LogicError") << "failed to load hls4ml model: a model is already loaded !";
  }

  // open the shared library containing the implementation of the model
  std::string const model_lib_name = model_name_ + ".so";
  model_lib_ = dlopen(model_lib_name.c_str(), RTLD_LAZY | RTLD_LOCAL);
  if (model_lib_ == nullptr) {
    throw cms::Exception("InvalidInput") << "hls4ml model library dlopen failure: cannot load library \""
                                         << model_lib_name << "\" !";
  }

  // "create_model" function: it creates an instance of the model and returns a pointer to the model
  create_model_cls* create_model = (create_model_cls*)dlsym(model_lib_, "create_model");
  if (dlerror()) {
    throw cms::Exception("InvalidInput") << "hls4ml emulator failed to load 'create_model' symbol from library \""
                                         << model_lib_name << "\" !";
  }

  // "destroy_model" function: deletes the model
  destroy_model_cls* destroy_model = (destroy_model_cls*)dlsym(model_lib_, "destroy_model");
  if (dlerror()) {
    throw cms::Exception("InvalidInput") << "hls4ml emulator failed to load 'destroy_model' symbol from library \""
                                         << model_lib_name << "\" !";
  }

  // smart pointer to the model with its own custom deleter
  model_ = model_ptr(create_model(), destroy_model);
  if (model_ == nullptr) {
    throw cms::Exception("InvalidInput") << "hls4ml emulator failed to load model (nullptr) from library \""
                                         << model_lib_name << "\" !";
  }
}
