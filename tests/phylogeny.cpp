#define CATCH_CONFIG_MAIN

#include "Catch2/single_include/catch2/catch.hpp"

#include "emp/Evolve/Systematics.hpp"

#include "phylogeny/phylogeny_utils.hpp"

TEST_CASE("Test phenotype_info::RecordPhenotype", "[Phylogeny]") {
  emp::vector<double> phenotype{0.0, 2.0, 0.0};
  emp::vector<bool> traits_evaluated{false, true, false};

  phylo::phenotype_info info;
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
  emp::vector<double> scores;
};

TEST_CASE("Test Phylogeny::NearestAncestorWithTraitEval", "[Phylogeny]") {

  SECTION("Simple tree") {
    using phylo_t = emp::Systematics<SimpleOrg, size_t, phylo::phenotype_info>;
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
    auto found_tax = phylo::NearestAncestorWithTraitEval(tax2, 0);
    REQUIRE(found_tax);
    REQUIRE(found_tax.value()->GetID() == 1);
    // Search for trait 0, set max dist such that search should succeed
    found_tax = phylo::NearestAncestorWithTraitEval(tax2, 0, 2);
    REQUIRE(found_tax);
    REQUIRE(found_tax.value()->GetID() == 1);
    // Search for trait 0, set max dist such that search should fail
    found_tax = phylo::NearestAncestorWithTraitEval(tax2, 0, 1);
    REQUIRE(!found_tax);
    // Search for trait 0, set max dist such that search should fail
    found_tax = phylo::NearestAncestorWithTraitEval(tax2, 0, 0);
    REQUIRE(!found_tax);

    // Search for trait 0, starting from tax1
    found_tax = phylo::NearestAncestorWithTraitEval(tax1, 0);
    REQUIRE(found_tax);
    REQUIRE(found_tax.value()->GetID() == 1);
    // Search for trait 0, starting from tax1
    found_tax = phylo::NearestAncestorWithTraitEval(tax1, 0, 1);
    REQUIRE(found_tax);
    REQUIRE(found_tax.value()->GetID() == 1);

    // Search for trait 0, starting from tax0
    found_tax = phylo::NearestAncestorWithTraitEval(tax0, 0);
    REQUIRE(found_tax);
    REQUIRE(found_tax.value()->GetID() == 1);
    // Search for trait 0, starting from tax0
    found_tax = phylo::NearestAncestorWithTraitEval(tax0, 0, 0);
    REQUIRE(found_tax);
    REQUIRE(found_tax.value()->GetID() == 1);

    // Search for trait 1, starting from tax2
    found_tax = phylo::NearestAncestorWithTraitEval(tax2, 1, 0);
    REQUIRE(!found_tax);
    // Search for trait 1, starting from tax2
    found_tax = phylo::NearestAncestorWithTraitEval(tax2, 1, 1);
    REQUIRE(found_tax);
    REQUIRE(found_tax.value()->GetID() == 2);

    // Search for trait 3, starting from tax2
    found_tax = phylo::NearestAncestorWithTraitEval(tax2, 3);
    REQUIRE(!found_tax);
    // Search for trait 3
    found_tax = phylo::NearestAncestorWithTraitEval(tax1, 3);
    REQUIRE(!found_tax);
    // Search for trait 3
    found_tax = phylo::NearestAncestorWithTraitEval(tax0, 3);
    REQUIRE(!found_tax);

    // if (found_tax) {
    //   std::cout << "Search results: " << found_tax.value()->GetID() << std::endl;
    // } else {
    //   std::cout << "Search failed" << std::endl;
    // }
  }

};

TEST_CASE("Test Phylogeny::NearestRelativeWithTraitEval", "[Phylogeny]") {

  /***
   * Tree:
   * 0 <= 1, 2, 3
   * 1 <= 4, 5, 6
   * 2 <= 7
   * 3 <= 8, 9
  */

  using phylo_t = emp::Systematics<SimpleOrg, size_t, phylo::phenotype_info>;
  using taxon_t = typename phylo_t::taxon_t;

  // Create empty phylogeny
  phylo_t phylo(
    [](SimpleOrg& org) { return org.org_info; }
  );

  // Initialize organisms, each organism starts with a single evaluated trait (id)
  // No organisms have been evaluted on the final trait
  emp::vector<SimpleOrg> orgs(10);
  emp::vector<size_t> taxon_ids(orgs.size(), 0);
  for (size_t i = 0; i < orgs.size(); ++i) {
    orgs[i].org_info = i;
    orgs[i].pos = i;
    orgs[i].scores = emp::vector<double>(orgs.size()+1, 0.0);
    orgs[i].traits_evaluated = emp::vector<bool>(orgs.size()+1, false);
    orgs[i].traits_evaluated[i] = true;
  }
  // Build the tree
  emp::vector<emp::Ptr<taxon_t>> taxa(orgs.size(), nullptr);
  taxa[0] = phylo.AddOrg(orgs[0], emp::WorldPosition(orgs[0].pos), nullptr);
  taxa[1] = phylo.AddOrg(orgs[1], emp::WorldPosition(orgs[1].pos), taxa[0]);
  taxa[2] = phylo.AddOrg(orgs[2], emp::WorldPosition(orgs[2].pos), taxa[0]);
  taxa[3] = phylo.AddOrg(orgs[3], emp::WorldPosition(orgs[3].pos), taxa[0]);
  taxa[4] = phylo.AddOrg(orgs[4], emp::WorldPosition(orgs[4].pos), taxa[1]);
  taxa[5] = phylo.AddOrg(orgs[5], emp::WorldPosition(orgs[5].pos), taxa[1]);
  taxa[6] = phylo.AddOrg(orgs[6], emp::WorldPosition(orgs[6].pos), taxa[1]);
  taxa[7] = phylo.AddOrg(orgs[7], emp::WorldPosition(orgs[7].pos), taxa[2]);
  taxa[8] = phylo.AddOrg(orgs[8], emp::WorldPosition(orgs[8].pos), taxa[3]);
  taxa[9] = phylo.AddOrg(orgs[9], emp::WorldPosition(orgs[9].pos), taxa[3]);
  // Record each taxon's phenotype to match associated organism
  std::unordered_map<size_t, size_t> taxonID_to_orgID;
  for (size_t i = 0; i < orgs.size(); ++i) {
    taxa[i]->GetData().RecordPhenotype(orgs[i].scores, orgs[i].traits_evaluated);
    taxon_ids[i] = taxa[i]->GetID();
    taxonID_to_orgID[taxa[i]->GetID()] = orgs[i].org_info;
  }

  // Search from organism 0 for trait 9
  auto found_tax = phylo::NearestRelativeWithTraitEval(taxa[0], 9);
  REQUIRE(found_tax);
  REQUIRE(taxonID_to_orgID[found_tax.value()->GetID()] == 9);
  // Search from organism 0 for trait 9, max dist 3
  found_tax = phylo::NearestRelativeWithTraitEval(taxa[0], 9, 3);
  REQUIRE(found_tax);
  REQUIRE(taxonID_to_orgID[found_tax.value()->GetID()] == 9);
  // Search from organism 0 for trait 9, max dist 2
  found_tax = phylo::NearestRelativeWithTraitEval(taxa[0], 9, 1);
  REQUIRE(!found_tax);

  // Search from organism 7 for trait 0
  found_tax = phylo::NearestRelativeWithTraitEval(taxa[7], 0);
  REQUIRE(found_tax);
  REQUIRE(taxonID_to_orgID[found_tax.value()->GetID()] == 0);
  // Search from organism 7 for trait 0, max dist 2
  found_tax = phylo::NearestRelativeWithTraitEval(taxa[7], 0, 1);
  REQUIRE(!found_tax);

  // Search from organism 5 to trait 7
  found_tax = phylo::NearestRelativeWithTraitEval(taxa[5], 7);
  REQUIRE(found_tax);
  REQUIRE(taxonID_to_orgID[found_tax.value()->GetID()] == 7);
  // Search from organism 5 to trait 7, max dist 4
  found_tax = phylo::NearestRelativeWithTraitEval(taxa[5], 7, 4);
  REQUIRE(found_tax);
  REQUIRE(taxonID_to_orgID[found_tax.value()->GetID()] == 7);
  // Search from organism 5 to trait 7, max dist 3
  found_tax = phylo::NearestRelativeWithTraitEval(taxa[5], 7, 3);
  REQUIRE(!found_tax);

  // Search from organism 1 to trait 3
  found_tax = phylo::NearestRelativeWithTraitEval(taxa[1], 3);
  REQUIRE(found_tax);
  REQUIRE(taxonID_to_orgID[found_tax.value()->GetID()] == 3);
  // Search from organism 5 to trait 7, max dist 4
  found_tax = phylo::NearestRelativeWithTraitEval(taxa[1], 3, 2);
  REQUIRE(found_tax);
  REQUIRE(taxonID_to_orgID[found_tax.value()->GetID()] == 3);
  // Search from organism 5 to trait 7, max dist 3
  found_tax = phylo::NearestRelativeWithTraitEval(taxa[1], 3, 1);
  REQUIRE(!found_tax);

  // Search from organism 0 to trait 10
  found_tax = phylo::NearestRelativeWithTraitEval(taxa[0], 10);
  REQUIRE(!found_tax);

  // Now org 2 and org 3 both evaluated on trait 3
  orgs[2].traits_evaluated[3] = true;
  taxa[2]->GetData().RecordPhenotype(orgs[2].scores, orgs[2].traits_evaluated);

  // Search from organism 1 to trait 3
  found_tax = phylo::NearestRelativeWithTraitEval(taxa[1], 3);
  REQUIRE(found_tax);
  REQUIRE( emp::Has({2, 3}, taxonID_to_orgID[found_tax.value()->GetID()]) );

  // Now only org 3 and 5 evaluted on trait 10
  orgs[3].traits_evaluated[10] = true;
  taxa[3]->GetData().RecordPhenotype(orgs[3].scores, orgs[3].traits_evaluated);
  orgs[5].traits_evaluated[10] = true;
  taxa[5]->GetData().RecordPhenotype(orgs[5].scores, orgs[5].traits_evaluated);
  found_tax = phylo::NearestRelativeWithTraitEval(taxa[0], 10);
  REQUIRE(found_tax);
  REQUIRE(taxonID_to_orgID[found_tax.value()->GetID()] == 3);

  // if (found_tax) {
  //   std::cout << "Search results (org id): " << taxonID_to_orgID[found_tax.value()->GetID()] << std::endl;
  // } else {
  //   std::cout << "Search failed" << std::endl;
  // }

};