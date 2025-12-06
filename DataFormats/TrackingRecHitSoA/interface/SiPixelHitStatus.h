#ifndef DataFormats_TrackingRecHitSoA_interface_SiPixelHitStatus_h
#define DataFormats_TrackingRecHitSoA_interface_SiPixelHitStatus_h

#include <cstdint>

class SiPixelHitStatus {
public:
    SiPixelHitStatus() : data_{0} {}

    bool isBigX() const { return (data_ >> 0) & 0b1; }
    bool isOneX() const { return (data_ >> 1) & 0b1; }
    bool isBigY() const { return (data_ >> 2) & 0b1; }
    bool isOneY() const { return (data_ >> 3) & 0b1; }
    uint8_t qBin() const { return (data_ >> 4) & 0b111; }

    void set_isBigX(bool const val) { set(val, 0b1, 0); }
    void set_isOneX(bool const val) { set(val, 0b1, 1); }
    void set_isBigY(bool const val) { set(val, 0b1, 2); }
    void set_isOneY(bool const val) { set(val, 0b1, 3); }
    void set_qBin(uint8_t const val) { set(val, 0b111, 4); }

private:
    void set(uint8_t const val, uint8_t const mask, uint8_t const pos) {
      data_ &= ~(mask << pos);
      data_ |= (val & mask) << pos;
    }

    uint8_t data_;
};

class SiPixelHitStatusAndCharge {
public:
    SiPixelHitStatusAndCharge() : charge_{0, 0, 0}, status_{} {}

    SiPixelHitStatus const& status() const { return status_; }
    SiPixelHitStatus& status() { return status_; }

    void set_status(SiPixelHitStatus const val) { status_ = val; }

    uint32_t charge() const {
      uint32_t ret{0};
      ret |= charge_[0];
      ret |= charge_[1] << 8;
      ret |= charge_[2] << 16;
      return ret;
    }

    void set_charge(uint32_t const val) {
      charge_[0] = kChargeMask & val;
      charge_[1] = kChargeMask & (val >> 8);
      charge_[2] = kChargeMask & (val >> 16);
    }

private:
    static constexpr uint8_t kChargeMask = 0b11111111;

    uint8_t charge_[3];
    SiPixelHitStatus status_;
};

#endif  // DataFormats_TrackingRecHitSoA_interface_SiPixelHitStatus_h
