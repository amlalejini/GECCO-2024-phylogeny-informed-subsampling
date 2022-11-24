#define CATCH_CONFIG_MAIN

#include "Catch2/single_include/catch2/catch.hpp"
#include "Phylogeny.hpp"

TEST_CASE("Test phenotype_info::RecordPhenotype", "[Phylogeny]") {
  emp::vector<double> phenotype{0.0, 2.0, 0.0};
  emp::vector<bool> traits_evaluated{false, true, false};

  phenotype_info info;
  info.RecordPhenotype(
    phenotype,
    traits_evaluated
  );

  // phenotype info should match what we recorded
  REQUIRE(phenotype == info.GetPhenotype());
  REQUIRE(traits_evaluated == info.GetTraitsEvaluated());

  // alter local variables, should not affect info
  phenotype[0] = 1.0;
  traits_evaluated[0] = true;
  REQUIRE(phenotype != info.GetPhenotype());
  REQUIRE(traits_evaluated != info.GetTraitsEvaluated());

};

struct SimpleOrg {
  size_t org_info=0;
  size_t pos=0;
  emp::vector<bool> traits_evaluated;
};

TEST_CASE("Test Phylogeny::NearestRelativeWithTraitEval", "[Phylogeny]") {

  SECTION("Simple tree") {
    using phylo_t = Phylogeny<SimpleOrg, size_t>;
    using taxon_t = typename phylo_t::taxon_t;
    // const size_t num_traits = 4;

    // Create empty phylogeny
    phylo_t phylo(
      [](SimpleOrg& org) { return org.org_info; }
    );

    SimpleOrg org0;
    SimpleOrg org1;
    SimpleOrg org2;
    org0.org_info = 0;
    org0.pos = 0;
    org0.traits_evaluated = {true, false, false, false};

    org1.org_info = 1;
    org1.pos = 1;
    org1.traits_evaluated = {false, true, true, false};

    org2.org_info = 2;
    org2.pos = 2;
    org2.traits_evaluated = {false, false, true, false};

    // TODO - bug in systematics? Providing no worldpos crashes?
    emp::Ptr<taxon_t> tax0 = phylo.AddOrg(org0, emp::WorldPosition(0), nullptr);
    emp::Ptr<taxon_t> tax1 = phylo.AddOrg(org1, emp::WorldPosition(1), tax0);
    emp::Ptr<taxon_t> tax2 = phylo.AddOrg(org2, emp::WorldPosition(2), tax1);
    tax0->GetData().RecordPhenotype({0,0,0}, org0.traits_evaluated);
    tax1->GetData().RecordPhenotype({0,0,0}, org1.traits_evaluated);
    tax2->GetData().RecordPhenotype({0,0,0}, org2.traits_evaluated);

    // Search for trait 0, starting from tax2.
    auto found_tax = phylo.NearestAncestorWithTraitEval(tax2, 0);
    REQUIRE(found_tax);
    REQUIRE(found_tax.value()->GetID() == 1);
    // Search for trait 0, set max dist such that search should succeed
    found_tax = phylo.NearestAncestorWithTraitEval(tax2, 0, 2);
    REQUIRE(found_tax);
    REQUIRE(found_tax.value()->GetID() == 1);
    // Search for trait 0, set max dist such that search should fail
    found_tax = phylo.NearestAncestorWithTraitEval(tax2, 0, 1);
    REQUIRE(!found_tax);
    // Search for trait 0, set max dist such that search should fail
    found_tax = phylo.NearestAncestorWithTraitEval(tax2, 0, 0);
    REQUIRE(!found_tax);

    // Search for trait 0, starting from tax1
    found_tax = phylo.NearestAncestorWithTraitEval(tax1, 0);
    REQUIRE(found_tax);
    REQUIRE(found_tax.value()->GetID() == 1);
    // Search for trait 0, starting from tax1
    found_tax = phylo.NearestAncestorWithTraitEval(tax1, 0, 1);
    REQUIRE(found_tax);
    REQUIRE(found_tax.value()->GetID() == 1);

    // Search for trait 0, starting from tax0
    found_tax = phylo.NearestAncestorWithTraitEval(tax0, 0);
    REQUIRE(found_tax);
    REQUIRE(found_tax.value()->GetID() == 1);
    // Search for trait 0, starting from tax0
    found_tax = phylo.NearestAncestorWithTraitEval(tax0, 0, 0);
    REQUIRE(found_tax);
    REQUIRE(found_tax.value()->GetID() == 1);

    // Search for trait 1, starting from tax2
    found_tax = phylo.NearestAncestorWithTraitEval(tax2, 1, 0);
    REQUIRE(!found_tax);
    // Search for trait 1, starting from tax2
    found_tax = phylo.NearestAncestorWithTraitEval(tax2, 1, 1);
    REQUIRE(found_tax);
    REQUIRE(found_tax.value()->GetID() == 2);

    // Search for trait 3, starting from tax2
    found_tax = phylo.NearestAncestorWithTraitEval(tax2, 3);
    REQUIRE(!found_tax);
    // Search for trait 3
    found_tax = phylo.NearestAncestorWithTraitEval(tax1, 3);
    REQUIRE(!found_tax);
    // Search for trait 3
    found_tax = phylo.NearestAncestorWithTraitEval(tax0, 3);
    REQUIRE(!found_tax);

    // if (found_tax) {
    //   std::cout << "Search results: " << found_tax.value()->GetID() << std::endl;
    // } else {
    //   std::cout << "Search failed" << std::endl;
    // }

  }

};