#ifndef ZP2MUANALYSIS_H
#define ZP2MUANALYSIS_H

#include <iosfwd>
#include <vector>
#include <string>

#include "TLorentzVector.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"

#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/OneToOne.h"
#include "DataFormats/Common/interface/OwnVector.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/RecLevelHelper.h"

namespace reco {
  // JMTBAD this is included in 170 and above
  typedef edm::RefToBase<reco::Track> TrackBaseRef;
}

// details about the number of quality cuts available
// JMTBAD which are not actually available yet
const int NUM_Q_SETS       = 8;
const int NUM_L2_CUTS      = 4;
const int NUM_L3_CUTS      = 9;
const int NUM_TRACKER_CUTS = 4;

// A simple class for passing around the generator-level particles
// from the resonant interaction.
struct InteractionParticles {
  InteractionParticles() : genQuark(0), genResonance(0),
			   genLepPlus(0), genLepMinus(0),
			   genLepPlusNoIB(0), genLepMinusNoIB(0) {}
  const reco::Candidate* genQuark;
  const reco::Candidate* genResonance;
  const reco::Candidate* genLepPlus;
  const reco::Candidate* genLepMinus;
  const reco::Candidate* genLepPlusNoIB;
  const reco::Candidate* genLepMinusNoIB;
};

class Zprime2muAnalysis : public edm::EDAnalyzer {
 public:
  explicit Zprime2muAnalysis(const edm::ParameterSet&);
  virtual ~Zprime2muAnalysis() {}

  virtual void beginJob(const edm::EventSetup&) {}
  virtual void endJob();

  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  static const double PTMIN;

 protected:
  // Verbosity levels.
  enum VERBOSITY { VERBOSITY_NONE, VERBOSITY_SIMPLE,
		   VERBOSITY_LOTS, VERBOSITY_TOOMUCH };

  // Type of hits to count in nHits().
  enum HITSTYPE { HITS_OTH=1, HITS_MU=1, HITS_EL=1,
		  HITS_TRK, HITS_PIX, HITS_SIL, HITS_ALL };

  // Lepton location codes, used by whereIs(Di)Lepton() methods.
  enum WhereLepton { W_BARREL=0, W_OVERLAP, W_ENDCAP, W_OUTSIDE };
  enum WhereDilepton { W_BARRELBARREL=0, W_BARRELOVERLAP,  W_BARRELENDCAP,
		       W_BARRELOUTSIDE,  W_OVERLAPOVERLAP, W_OVERLAPENDCAP, 
		       W_OVERLAPOUTSIDE, W_ENDCAPENDCAP,   W_ENDCAPOUTSIDE,
		       W_OUTSIDEOUTSIDE };
  
  // Some convenient typedefs.
  typedef reco::Particle::LorentzVector LorentzVector;
  typedef std::vector<reco::CandidateBaseRef> LeptonRefVector;

  // Hard-coded parameters defined in the .cc file.
  static const unsigned int MAX_DILEPTONS;
  static const double ETA_CUT;
  static const double ENDCAP_BARREL_CUT;
  static const double TRIGGER_ETA_CUT[NUM_REC_LEVELS];
  static const bool DO_QCUTS;
  static const int QSEL;
  static const double L2QCUT[NUM_Q_SETS][NUM_L2_CUTS];
  static const double L3QCUT[NUM_Q_SETS][NUM_L3_CUTS];
  static const double TRACKERQCUT[NUM_Q_SETS][NUM_TRACKER_CUTS];
  static bool cutTrig[NUM_REC_LEVELS];

  ////////////////////////////////////////////////////////////////////
  // Parameters read or determined from the config file:
  ////////////////////////////////////////////////////////////////////
  
  // whether to date the histogram output postscript pages
  bool dateHistograms;

  // whether we are looking at electrons instead of muons;
  bool doingElectrons;

  // determines whether one or two dileptons are kept (useful for
  // H->ZZ studies);
  bool doingHiggs;

  // whether to allow construction of generator-level dileptons from
  // generated leptons -- default behavior is no, and to take the
  // actual resonance from the PYTHIA event record, but if there is
  // none as is the case in some COMPHEP-generated samples, then the
  // code will not find any generated dilepton unless this parameter
  // is true;
  bool constructGenDil;

  // whether to look at only generator-level muons (i.e. don't bother
  // trying to store globalMuons, standAloneMuons, etc);
  bool generatedOnly;

  // whether to look at Geant4 particles in addition to Pythia particles;
  bool doingGeant4;

  // whether to look at only reconstructed muons (as in the real data sets);
  bool reconstructedOnly;

  // whether to include the extra muon reconstructors (FMS, PMR, etc);
  bool useOtherMuonRecos;

  // whether to use the TMR cocktail or Piotr's
  bool useTMRforBest;

  // if the input file is only AOD, then don't use EDProducts that are
  // included only in the full RECO tier;
  bool usingAODOnly;

  // whether trigger information is supposed to be present in the
  // input file;
  bool useTriggerInfo;

  // basic quantities for the chosen lepton (muon or electron);
  unsigned int leptonFlavor; // PDG ID
  double leptonMass;         // in GeV/c^2

  // and trigger path information.
  std::vector<l1extra::L1ParticleMap::L1TriggerType> l1paths;
  std::vector<std::string> hltModules[2]; // in order: L2, L3
  std::vector<std::string> hltPaths;

  ////////////////////////////////////////////////////////////////////
  // Event data
  ////////////////////////////////////////////////////////////////////
  
  // Keep track of the event number for printing out.
  int eventNum;

  // Trigger results.
  bool passTrig[NUM_REC_LEVELS];
  unsigned int trigWord[NUM_REC_LEVELS];

  // Keep the collection of generator-level particles around; it is
  // used in many places.
  const reco::CandidateCollection* genParticles;

  // Stored information from the MC record about the hard interaction,
  // which is assumed to be of the form q qbar -> resonance ->
  // l+l-. Store the quark, the resonance, the final-state leptons,
  // and the leptons before bremsstrahlung (the last two that are
  // suffixed with NoIB). We explicitly assume that there is only one
  // resonance per event.
  InteractionParticles intParticles;

  // The lepton collections, stored as CandidateBaseRefs.
  // bestLeptons, at least, cannot be a CandidateBaseRefVector since we
  // want to mix CandidateBaseRefs from different collections with
  // different product ids. So be consistent and make the rest of the
  // lepton collections not CandidateBaseRefVectors.
  LeptonRefVector allLeptons[MAX_LEVELS];
  LeptonRefVector bestLeptons;

  // The dilepton collections, stored as CompositeCandidates in an
  // OwnVector (for compatibility with CandCombiner, if we switch to
  // using it.)
  reco::CandidateCollection allDileptons[MAX_LEVELS];
  reco::CandidateCollection bestDileptons;

  // Store separately the "resonance" four-vectors, i.e. leptons plus
  // photons, since the four-vectors in allDileptons are from just
  // leptons. dileptonResonances[MAX_LEVELS] are bestDileptons'
  // resonances.
  std::vector<LorentzVector> dileptonResonances[MAX_LEVELS+1];
  
  // Map product IDs from CandidateBaseRef to whatever rec level we
  // choose.
  std::map<int,int> recLevelMap;
  
  // Map lepton CandidateBaseRefs to photon ones.
  reco::CandViewMatchMap photonMatchMap[MAX_LEVELS];

  // Store seed indices for global fits (indexed in the vector by
  // id(cand)).
  std::vector<int> seedIndices[MAX_LEVELS];

  // Store an AssociationMap for each rec level which maps leptons to
  // leptons at other rec levels, matching by closest in delta R. (We
  // are wasteful, not using entries for which the rec levels are
  // equal, but storage is easier this way.)
  reco::CandViewMatchMap closestMatchMap[MAX_LEVELS][MAX_LEVELS];

  // Store an AssociationMap for each rec level of global fits
  // (i.e. L3 and above), matching by "seed" lepton. (Currently only
  // meaningful for muons, and we are even more wasteful than before,
  // in addition not using entries for which rec levels are less than
  // l3.)
  reco::CandViewMatchMap seedMatchMap[MAX_LEVELS][MAX_LEVELS];

  // An invalid CandidateBaseRef so that some of the methods above can
  // return by reference and still be able to return a null reference.
  reco::CandidateBaseRef invalidRef;

  // A handle to the TFileService for convenience.
  edm::Service<TFileService> fs;

 private:
  // Object to help with per-rec-level information.
  RecLevelHelper recLevelHelper;

  // verbosity controls the amount of debugging information printed;
  // levels are defined using the VERBOSITY_* codes above
  VERBOSITY verbosity;
  edm::InputTag l1ParticleMap;
  edm::InputTag hltResults;
  edm::InputTag standAloneMuons;
  edm::InputTag photons;

  edm::Handle<reco::TrackCollection> seedTracks;
  edm::Handle<reco::PhotonCollection> photonCollection;

  ////////////////////////////////////////////////////////////////////
  // Initialization
  ////////////////////////////////////////////////////////////////////

  // Initialize ROOT settings to how we like them.
  void InitROOT();

  // Make sure all the various lepton vectors and trigger bitmaps,
  // etc., are cleared.
  void clearValues();

  ////////////////////////////////////////////////////////////////////
  // Storing leptons/dileptons
  ////////////////////////////////////////////////////////////////////

  // Store a sorted vector of CandidateBaseRefs to leptons at level rec,
  // returning a success flag.
  bool storeLeptons(const edm::Event& event, const int rec);

  // Construct a new dilepton from the two specified lepton
  // candidates, using AddFourMomenta to set up the dilepton.
  reco::CompositeCandidate* newDilepton(const reco::CandidateBaseRef& dau1,
					const reco::CandidateBaseRef& dau2);

  // Prune a dilepton collection: if more than one dilepton was
  // formed, we accept only those containing distinct leptons. The
  // preference is given to higher mass dilepton over lower mass
  // dilepton (i.e. to the earlier dilepton in the current sort
  // order).
  void removeDileptonOverlap(reco::CandidateCollection& dileptons);

  // Consider all possible combinations of l+ and l-, and form
  // dilepton candidates. 
  void makeDileptons(const int rec);

  // If there is a photon candidate at the phi-eta distance dRmax from 
  // either of the leptons, combine its 4-momentum with that of the dilepton
  // and store it in the collection of resonance vectors.
  void addBremCandidates(const int rec);

  // Instead of adding brem candidates for the generator-level
  // dilepton, take the true resonance value from e.g. PYTHIA.
  void addTrueResonance(const edm::Event& event);

  ////////////////////////////////////////////////////////////////////
  // Using trigger info
  ////////////////////////////////////////////////////////////////////

  // Get Level-1 decisions for trigger paths we are interested in,
  // storing them in a bitmap.
  void storeL1Decision(const edm::Event& event);
  // Same idea, but for levels 2 and 3 of the HLT.
  void storeHLTDecision(const edm::Event& event);

  // Function to translate the algorithm algo defined for trigger
  // level lvl, requiring nmu muons (e.g. single or dimuon trigger),
  // and return whether we judge the event to pass the trigger at this
  // level for this algorithm.
  bool TriggerTranslator(const std::string& algo, const unsigned int lvl, 
			 const unsigned int nmu) const;

  // Use TriggerTranslator to compare the "official" decisions stored
  // by store*Decision with what we calculate.
  void compareTrigDecision(const edm::Event& event, bool old=false) const;
  
  ////////////////////////////////////////////////////////////////////
  // Picking "best" leptons
  ////////////////////////////////////////////////////////////////////

  // Our implementation of cocktail methods for picking muons from the
  // various TeV muon reconstructors (either TMR, picking between
  // tracker-only and tracker+first muon station, or Piotr's, picking
  // between those two and also PMR); returns a reference to the one
  // picked.
  const reco::CandidateBaseRef&
    cocktailMuon(const reco::CandidateBaseRef& trk,
		 const bool doTMR, const bool debug) const;

  // Keep some statistics on what the cocktail picked.
  int best_ntrk, best_ngmr, best_nfms, best_npmr, best_ngpr, best_ntot;

  // A driver routine which uses the cocktail methods above to pick
  // "best" leptons (only implemented for muons).
  void findBestLeptons();

 protected:
  ////////////////////////////////////////////////////////////////////
  // General utility functions
  ////////////////////////////////////////////////////////////////////

  // Build a string containing a name for a histogram -- useful inside
  // loops when booking histos.
  std::string nameHist(const char* s, int i, int j = -1) {
    std::string x = std::string(s) + char(i + 48);
    if (j >= 0) x += char(j + 48);
    return x;
  }

  ////////////////////////////////////////////////////////////////////
  // Lepton/dilepton rec level utility functions
  ////////////////////////////////////////////////////////////////////

  // Return a unique id for the lepton cand (currently implemented by
  // the reference's index into the collection).
  int id(const reco::CandidateRef& cand) const;
  int id(const reco::CandidateBaseRef& cand) const;

  // Translate the Ref's product id to one of our rec levels, using the
  // cached map.
  int recLevel(const reco::CandidateBaseRef& cand) const
  { return recLevelHelper.recLevel(cand); }

  // Get the rec level for a dilepton, making sure that either all the
  // daughter leptons have the same rec level, or else returning lbest,
  // (since a "best" dilepton can be made up of leptons at different rec
  // levels).
  int recLevel(const reco::Candidate& cand) const;

  // Provide uniform access to allLeptons and bestLeptons by returning a
  // reference to the appropriate collection.
  const LeptonRefVector& getLeptons(const int rec) const;

  // Same as getLeptons() but for dileptons.
  const reco::CandidateCollection& getDileptons(const int rec) const;

  // Get the four-vector of the closest photon found for cand.
  LorentzVector closestPhoton(const reco::CandidateBaseRef& cand) const
  { return recLevelHelper.closestPhoton(cand); }

  // Get the seed index (i.e. the index into the stand-alone muon
  // collection) of the candidate.
  int seedIndex(const reco::CandidateBaseRef& cand) const 
  { return recLevelHelper.seedIndex(cand); }

  ////////////////////////////////////////////////////////////////////
  // Lepton/dilepton matching
  ////////////////////////////////////////////////////////////////////

  // Const access to lep's closest match in delta R at another rec level.
  const reco::CandidateBaseRef&
    closestLepton(const reco::CandidateBaseRef& lep,
		  const int level) const
    { return recLevelHelper.closestLepton(lep, level); }

  // Const access to lep's same-seed match at another rec level.
  const reco::CandidateBaseRef&
    sameSeedLepton(const reco::CandidateBaseRef& lep,
		   const int level) const
    { return recLevelHelper.sameSeedLepton(lep, level); }

  // Const access to lep's "best" match, i.e. same-seed if available,
  // closest in delta R if not, at another rec level.
  const reco::CandidateBaseRef&
    matchedLepton(const reco::CandidateBaseRef& lep,
		  const int level) const
    { return recLevelHelper.matchedLepton(lep, level); }

  // Try to find the dilepton at the new rec level that has the same
  // two leptons as dil. Return success, and return the other dilepton in
  // newdil.
  bool matchDilepton(const reco::Candidate& dil,
		     const int level,
		     const reco::Candidate* newdil) const;

  ////////////////////////////////////////////////////////////////////
  // Generator-level utility functions
  ////////////////////////////////////////////////////////////////////

  // Work around Cartesian v. polar coordinates for LorentzVectors for now.
  void SetP4M(reco::Particle::LorentzVector& v,
	      double pt, double phi, double p,
	      double theta, double mass) const;

  // Return the lepton cand's mother. If cand's mother is the same
  // lepton but before brem, go up the decay chain (what we call the
  // "non-brem mother"). This is mainly to avoid stopping at leptons
  // in documentation lines, which are declared to be ancestors of
  // muons produced in hard interaction (i.e., ancestors of
  // themselves).
  const reco::Candidate* mother(const reco::CandidateBaseRef& cand) const;

  // Return the PDG id of the non-brem mother of cand. If the mother
  // pointer isn't valid, return pdgId = 0.
  int motherId(const reco::CandidateBaseRef& cand) const;
  
  // Return the PDG id of the grandmother of cand (i.e. a quark or
  // antiquark, or a gluon). If the mother or grandmother pointer isn't
  // valid, return pdgId = 0.
  int grandmotherId(const reco::CandidateBaseRef& cand) const;

  // If cand1 and cand2 have the same non-brem mothers, return a
  // pointer to the mother candidate, else return null.
  const reco::Candidate* sameMother(const reco::CandidateBaseRef& cand1,
				    const reco::CandidateBaseRef& cand2) const;

  // Test to see if the pdg ID is that of a resonance we want to
  // analyze, currently one of Z0 (i.e. Drell-Yan), Z', or G*.
  bool isResonance(int pid) const;

  // Store the generator-level momenta for the particles in the resonant
  // interaction from the MC record (assuming the event was of the form
  // q qbar -> resonance -> l+l-); return success as a bool. The returned
  // momenta include the true resonance, the final-state (after-brem)
  // l-l+, the muons before bremsstrahlung, and the quark that entered
  // the hard interaction.
  bool storeInteractionParticles(const reco::CandidateCollection& genParticles,
				 int eventNum, InteractionParticles& ip) const;

  ////////////////////////////////////////////////////////////////////
  // Track utility functions
  ////////////////////////////////////////////////////////////////////

  // Get the "main" track associated with the lepton. The definition
  // of "main" depends on the lepton type and its rec level. For electrons,
  // only the GsfTrack is appropriate for now. For muons at each rec level:
  //   L2: return the muon system track
  //   L3 or higher: return the combined track (tracker+muon system)
  // If something doesn't make sense (e.g. requesting an L1 track for
  // which a Track is not reconstructed), return an invalid reference.
  reco::TrackBaseRef getMainTrack(const reco::CandidateBaseRef& cand) const;

  // There is no pError() in Track/TrackBase; calculate error on p
  // ourselves from error on qoverp().
  template <typename TrackType> double pError(const TrackType& cand) const;
  // There is a ptError() in Track/TrackBase but for symmetry with the
  // above let's have another.
  template <typename TrackType> double ptError(const TrackType& cand) const;
  
  // Propagate inverse errors...
  double invError(double val, double err) const { return 1/val/val*err; } 
  // ... to 1/pT
  template <typename TrackType> double invPtError(const TrackType& tk) const;
  // ... and 1/p
  template <typename TrackType> double invPError(const TrackType& tk) const;

  // Provide the same track error methods to be called directly on
  // candidates, for convenience.
  double ptError(const reco::CandidateBaseRef& cand) const;
  double pError(const reco::CandidateBaseRef& cand) const;
  double invPtError(const reco::CandidateBaseRef& cand) const;
  double invPError(const reco::CandidateBaseRef& cand) const;
  
  // Return by reference the number of rec hits counted by pixel,
  // silicon, and "other-system", which could be the muon system or
  // the ecal for electrons.
  template <typename TrackType>
    void numSubDetHits(const TrackType& theTrack,
		       int& nPixHits, int& nSilHits, int& nOtherHits,
		       DetId::Detector otherDet) const;

  // Get the number of hits on the appropriate track (specified by
  // type using the codes HITS_* above) depending on the lepton type
  // and rec level.
  int nHits(const reco::CandidateBaseRef& cand, const int type) const;

  // Return whether cand1 and cand2 are "close"; instead of a circle
  // in eta-phi space, we look at a square .5 on a side.
  bool matchTracks(const reco::CandidateBaseRef& cand1,
		   const reco::CandidateBaseRef& cand2) const;

  // Mimic what HLTMuonPrefilter does for pt cut: convert 50%
  // efficiency threshold to 90% efficiency threshold.
  double ptLx(const reco::TrackRef& theTrack, const int rec) const;

  // Simplify usage of chi2 prob for tracks (used in selecting "best"
  // muons).
  double cumulativeChiSquare(const reco::CandidateBaseRef& mu) const;

  ////////////////////////////////////////////////////////////////////
  // Dilepton utility functions
  ////////////////////////////////////////////////////////////////////

  // Get the "resonance" four-vector, i.e. the dilepton four-vector
  // plus any associated photons.
  const LorentzVector& resV(const int rec, const int idil) const;

  // Count the number of daughters the dilepton has in the specified
  // acceptance in eta.
  int numDaughtersInAcc(const reco::Candidate& dil, const double etaCut) const;

  // Return a reference to the ith daughter lepton of the dilepton, or
  // an invalid reference if i is out of bounds.
  const reco::CandidateBaseRef
    dileptonDaughter(const reco::Candidate& dil,
		     const unsigned i) const;

  // Return a reference to the daughter lepton of the dilepton with
  // specified charge (if it is a same-sign dilepton, will return the
  // first one found), or else an invalid reference if not found.
  const reco::CandidateBaseRef
    dileptonDaughterByCharge(const reco::Candidate& dil,
			     const int charge) const;

  ////////////////////////////////////////////////////////////////////
  // Print-outs
  ////////////////////////////////////////////////////////////////////

  // Print out all the relevant information about the lepton; but this
  // method is just as useful as documentation on how to access this
  // information.
  void dumpLepton(std::ostream& output, reco::CandidateBaseRef cand) const;

  // Dump the masses of the dileptons formed at each level of
  // reconstruction.
  void dumpDileptonMasses() const;

  // Dump all quality information about a dilepton's daughter lepton's
  // tracks.
  void dumpDilQuality() const;

  // Dump the event, printing out the specified information at each
  // level of lepton reconstruction.
  void dumpEvent(const bool printGen = false, const bool printL1 = false,
                 const bool printL2 = false, const bool printL3 = false,
                 const bool printBest = false,
		 const bool printSeeds = false) const;

  ////////////////////////////////////////////////////////////////////
  // Quality cuts
  ////////////////////////////////////////////////////////////////////

  // Quality cuts on tracks, not yet implemented.
  bool TrackQCheck(const reco::CandidateBaseRef& lepton, const int qsel,
		   int& ncut) const;

  // Apply the above track quality cuts to each of the daughters of
  // the dilepton, not yet implemented.
  bool dilQCheck(const reco::Candidate& dilepton, const int qsel,
		 int& ncut_mum, int& ncut_mup) const;

  ////////////////////////////////////////////////////////////////////
  // Trigger
  ////////////////////////////////////////////////////////////////////

  // Return whether the event passed the trigger at irec level.
  bool passTrigger(const int irec) const;

  // Return whether the event passed the entire trigger (L1+HLT).
  bool passTrigger() const;

  ////////////////////////////////////////////////////////////////////
  // Lepton/dilepton acceptance
  ////////////////////////////////////////////////////////////////////

  // Helper method to return a code (one of the W_* ones above) based
  // on where the lepton is in the muon system by eta; in the barrel,
  // in the overlap region, in the endcap, or outside acceptance
  // (nominally 2.4).
  // JMTBAD are the codes for electrons useful?
  WhereLepton whereIsLepton(const reco::CandidateBaseRef& lepton);

  // Helper method to return a code based on where the leptons of a
  // dilepton are in the lepton system (using the above whereIsLepton method
  // definitions of location).
  WhereDilepton whereIsDilepton(const reco::Candidate& dil);

  ////////////////////////////////////////////////////////////////////
  // Analysis level function (cuts, etc.)
  ////////////////////////////////////////////////////////////////////

  // Return whether the event is "interesting", used in determining
  // whether to override the verbosity level for dumping the event.
  bool eventIsInteresting();

  // Checks our defined cut methods and returns a bitmap of which
  // passed. This can be directly used as a comparison, since if no
  // cuts are made, then the return value is 0 == false; otherwise it
  // is > 0 == true.
  unsigned leptonIsCut(const reco::CandidateBaseRef& lepton);
};

// Sorting functors.
struct reverse_momentum_sort {
  bool operator()(const reco::Candidate& lhs, const reco::Candidate& rhs) const {
    return lhs.p() > rhs.p();
  }
  bool operator()(const reco::CandidateBaseRef& lhs, const reco::CandidateBaseRef& rhs) const {
    return lhs->p() > rhs->p();
  }
};

struct reverse_mass_sort {
  bool operator()(const reco::Candidate& lhs, const reco::Candidate& rhs) const {
    return lhs.mass() > rhs.mass();
  }
};

// Functions to cast base types of Candidates to concrete derived types
// e.g. to a reco::Muon.
template <typename T>
inline const T& toConcrete(const reco::Candidate& cand) {
  return *dynamic_cast<const T*>(&cand);
}

template <typename T>
inline const T& toConcrete(const reco::CandidateRef& cand) {
  return *dynamic_cast<const T*>(&*cand);
}

template <typename T>
inline const T& toConcrete(const reco::CandidateBaseRef& cand) {
  return *dynamic_cast<const T*>(&*cand);
}

// For pretty-printing.
std::ostream& operator<<(std::ostream& out, const TLorentzVector& vect);
std::ostream& operator<<(std::ostream& out, const reco::Candidate& par);

#endif // ZP2MUANALYSIS_H
