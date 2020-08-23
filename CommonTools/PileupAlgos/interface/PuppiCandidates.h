#ifndef CommonTools_PileupAlgos_PuppiCandidates_h
#define CommonTools_PileupAlgos_PuppiCandidates_h

#include <vector>
#include <cassert>

#include "FWCore/SOA/interface/Column.h"
#include "FWCore/SOA/interface/Table.h"
#include "CommonTools/Utils/interface/KinematicColumns.h"

namespace edm {

  namespace soa {

    namespace col {

      SOA_DECLARE_COLUMN(Rapidity, float, "rapidity");
      SOA_DECLARE_COLUMN(Id, int, "id");

    }

    using PtEtaRapPhiIdTable = edm::soa::Table<col::Pt, col::Eta, col::Rapidity, col::Phi, col::Id>;

    template <class Object>
      PtEtaRapPhiIdTable makePtEtaRapPhiIdTable(std::vector<Object> const& objects) {
      return {objects,
	  edm::soa::column_fillers(
				   col::Pt::filler([](Object const& x) { return x.pt; }),
				   col::Eta::filler([](Object const& x) { return x.eta; }),
				   col::Rapidity::filler([](Object const& x) { return x.rapidity; }),
				   col::Phi::filler([](Object const& x) { return x.phi; }),
				   col::Id::filler([](Object const& x) { return x.id; })
				   )};
    }

  }  // namespace soa

}  // namespace edm

#endif
