#ifndef L1Trigger_MLUtilities_HLS4MLModelWrapper_h
#define L1Trigger_MLUtilities_HLS4MLModelWrapper_h

#include <any>
#include <functional>
#include <memory>
#include <string>

namespace l1t {

  class HLS4MLModel {
  public:
    virtual void prepare_input(std::any input) = 0;
    virtual void predict() = 0;
    virtual void read_result(std::any result) = 0;
    virtual ~HLS4MLModel() = default;
  };

  class HLS4MLModelWrapper {
  public:
    HLS4MLModelWrapper();
    HLS4MLModelWrapper(std::string const& model_name);
    ~HLS4MLModelWrapper();

    HLS4MLModelWrapper(HLS4MLModelWrapper const&) = delete;
    HLS4MLModelWrapper& operator=(HLS4MLModelWrapper const&) = delete;
    HLS4MLModelWrapper(HLS4MLModelWrapper&&) = delete;
    HLS4MLModelWrapper& operator=(HLS4MLModelWrapper&&) = delete;

    void reset();
    void reset(std::string const& model_name);

    void prepare_input(std::any input) const;
    void predict() const;
    void read_result(std::any result) const;

    std::string const& model_name() const { return model_name_; }

    bool model_loaded() const { return (model_ != nullptr); }

  private:
    using create_model_cls = HLS4MLModel*();
    using destroy_model_cls = void(HLS4MLModel*);
    using model_ptr = std::unique_ptr<HLS4MLModel, std::function<destroy_model_cls>>;

    void load();

    std::shared_ptr<void> model_lib_;
    model_ptr model_;
    std::string model_name_;
  };

}  // namespace l1t

#endif
