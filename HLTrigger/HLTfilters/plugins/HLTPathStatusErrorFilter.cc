#include <memory>
#include <string>
#include <vector>
#include <regex>
#include <unordered_map>
#include <iomanip>
#include <algorithm>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/RegexMatch.h"
#include "FWCore/Utilities/interface/transform.h"
#include "DataFormats/Common/interface/HLTPathStatus.h"

class HLTPathStatusErrorFilter : public edm::stream::EDFilter<> {
public:
  explicit HLTPathStatusErrorFilter(edm::ParameterSet const&);
  ~HLTPathStatusErrorFilter() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  bool filter(edm::Event&, edm::EventSetup const&) override;

private:
  void beginStream(edm::StreamID) override;

  auto logTrace() const {
    auto const& moduleType = moduleDescription().moduleName();
    auto const& moduleLabel = moduleDescription().moduleLabel();
    return LogTrace(moduleType) << "[" << moduleType << "] (" << moduleLabel << ") ";
  }

  struct PatternData {
    PatternData(std::string const& aStr, std::regex const& aRegex, bool const hasMatch = false)
        : str(aStr), regex(aRegex), matched(hasMatch) {}
    std::string str;
    std::regex regex;
    bool matched;
  };

  bool const ignoreInvalidPathNames_;
  std::vector<PatternData> hltPathStatusKeepPatterns_;
  std::vector<PatternData> hltPathStatusDropPatterns_;
  std::unordered_map<std::string, edm::EDGetTokenT<edm::HLTPathStatus>> hltPathStatusTokensMap_;
};

HLTPathStatusErrorFilter::HLTPathStatusErrorFilter(edm::ParameterSet const& config)
    : ignoreInvalidPathNames_{config.getParameter<bool>("ignoreInvalidPathNames")} {
  hltPathStatusKeepPatterns_ =
      edm::vector_transform(config.getParameter<std::vector<std::string>>("pathNames"), [](std::string const& pattern) {
        return PatternData(pattern, std::regex(edm::glob2reg(pattern), std::regex::extended));
      });

  hltPathStatusDropPatterns_ = edm::vector_transform(
      config.getParameter<std::vector<std::string>>("pathNamesToSkip"), [](std::string const& pattern) {
        return PatternData(pattern, std::regex(edm::glob2reg(pattern), std::regex::extended));
      });

  if (hltPathStatusKeepPatterns_.empty()) {
    return;
  }

  // consume all matching Paths
  callWhenNewProductsRegistered([this](edm::BranchDescription const& branch) {
    if (branch.branchType() == edm::InEvent and branch.className() == "edm::HLTPathStatus") {
      auto const& pathName = branch.moduleLabel();
      for (auto& patternKeep : hltPathStatusKeepPatterns_) {
        if (std::regex_match(pathName, patternKeep.regex)) {
          logTrace() << "Path named \"" << pathName << "\" matches input pattern \"" << patternKeep.str << "\"";
          patternKeep.matched = true;

          bool skipPath = false;
          for (auto& patternDrop : hltPathStatusDropPatterns_) {
            if (not std::regex_match(pathName, patternDrop.regex)) {
              continue;
            }

            logTrace() << "  but Path named \"" << pathName << "\" is vetoed by pattern \"" << patternDrop.str << "\"";
            patternDrop.matched = true;

            skipPath = true;
#ifndef EDM_ML_DEBUG
            break;
#endif
          }

          if (not skipPath and hltPathStatusTokensMap_.find(pathName) == hltPathStatusTokensMap_.end()) {
            hltPathStatusTokensMap_[pathName] = consumes<edm::HLTPathStatus>(
                edm::InputTag(pathName, branch.productInstanceName(), branch.processName()));
          }
        }
      }
    }
  });
}

void HLTPathStatusErrorFilter::beginStream(edm::StreamID) {
  // if throw=True, check if any of the input patterns had zero matches (and if so, throw an exception)
  if (not(hltPathStatusKeepPatterns_.empty() or ignoreInvalidPathNames_)) {
    auto const unmatchedPatternsExist = std::any_of(hltPathStatusKeepPatterns_.cbegin(),
                                                    hltPathStatusKeepPatterns_.cend(),
                                                    [](auto foo) { return (not foo.matched); });
    if (unmatchedPatternsExist) {
      cms::Exception excpt("HLTPathStatusErrorFilterInvalidPathName");
      excpt << "the parameter \"pathNames\" contains patterns with zero matches"
            << " for the available edm::HLTPathStatus collections - invalid patterns are:";
      for (auto const& pattern : hltPathStatusKeepPatterns_)
        if (not pattern.matched)
          excpt << "\n\t" << pattern.str;
      throw excpt;
    }
  }
}

void HLTPathStatusErrorFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::vector<std::string>>("pathNames", {});
  desc.add<std::vector<std::string>>("pathNamesToSkip", {});
  desc.add<bool>("ignoreInvalidPathNames", false);
  descriptions.addWithDefaultLabel(desc);
}

bool HLTPathStatusErrorFilter::filter(edm::Event& event, edm::EventSetup const& setup) {
  bool ret = false;
  for (auto const& [pathName, hltPathStatusToken] : hltPathStatusTokensMap_) {
    auto const& hltPathStatusHandle = event.getHandle(hltPathStatusToken);
    if (not hltPathStatusHandle.isValid()) {
      continue;
    }

    logTrace() << "Path=" << pathName << " state=" << hltPathStatusHandle->state()
               << " wasrun=" << hltPathStatusHandle->wasrun() << " accept=" << hltPathStatusHandle->accept()
               << " error=" << hltPathStatusHandle->error();

    if (hltPathStatusHandle->error()) {
      ret = true;
#ifndef EDM_ML_DEBUG
      break;
#endif
    }
  }

  logTrace() << "filter decision = " << ret;

  return ret;
}

// register as framework plugin
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HLTPathStatusErrorFilter);
