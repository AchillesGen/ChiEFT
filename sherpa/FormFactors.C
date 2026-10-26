// FormFactors.C -- Sherpa form factor plugin for the ChiPT UFO model.
//
// Sherpa-generate-model does not write this file.  Copy it into
// ChiPT_UFO/.sherpa/ and register it in .sherpa/CMakeLists.txt, which
// sherpa/build_model.sh does.
//
// Each UFO form factor calls one of the functions below with
//   arg = Q^2 = -q^2 of the gauge boson      (all except FDeltaLam)
//   arg = p^2 of the Delta leg               (FDeltaLam)
//
// All coefficients are read from block CHIPTFF of the param card, so they can
// be changed without recompiling.  Values are cached on first use: start a new
// run after editing the card.  Do not compile with -ffast-math (the per-call
// cache uses a NaN sentinel).

#include "ATOOLS/Org/Message.H"
#include "METOOLS/Explicit/Form_Factor.H"
#include "METOOLS/Explicit/Vertex_Key.H"
#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Model_Base.H"

#include "ATOOLS/Org/Exception.H"

#include <cmath>
#include <limits>
#include <map>

using namespace ATOOLS;
using namespace METOOLS;

namespace {

  // Masses in GeV (nucleon, neutral pion, Delta pole).
  constexpr double MN = 0.938272, MPi = 0.1349768, MD = 1.209;
  // Charged pion mass, for the charged-current pion poles.
  constexpr double MPiC = 0.1395704;

  // ---- param card access ---------------------------------------------------
  struct Params {
    // Required entry.  Sherpa returns 0 for a missing entry, which would
    // silently zero the cross section, so fail loudly instead.
    static double Get(const std::string &name) {
      const auto &c = MODEL::s_model->ScalarConstants();
      auto it = c.find(name);
      if (it == c.end() || it->second == 0.0)
        THROW(fatal_error, "CHIPTFF parameter '" + name + "' missing or zero "
              "in the param card.");
      return it->second;
    }
    // Optional entry with a default.
    static double Opt(const std::string &name, double dflt) {
      const auto &c = MODEL::s_model->ScalarConstants();
      auto it = c.find(name);
      return it == c.end() ? dflt : it->second;
    }
    // The model's own sin^2(theta_W) = 1 - MW^2/MZ^2 (not Sw2FF).
    static double Sw2Model() {
      static const double sw2 = MODEL::s_model->ComplexConstant("csin2_thetaW").real();
      return sw2;
    }
  };

  // ---- nucleon vector form factors (Kelly) ---------------------------------
  inline double KellyRatio(double a1, double b1, double b2, double b3, double x) {
    return (1 + a1 * x) / (1 + b1 * x + b2 * x * x + b3 * x * x * x);
  }

  // Sachs form factors GE, GM for proton and neutron.
  void Sachs(double Q2, double &Gep, double &Gen, double &Gmp, double &Gmn) {
    static const double ep1 = Params::Get("KEp1"), ep2 = Params::Get("KEp2"),
                        ep3 = Params::Get("KEp3"), ep4 = Params::Get("KEp4"),
                        lam = Params::Get("KLamSq"),
                        en1 = Params::Get("KEn1"), en2 = Params::Get("KEn2"),
                        mup = Params::Get("MuP"),
                        mp1 = Params::Get("KMp1"), mp2 = Params::Get("KMp2"),
                        mp3 = Params::Get("KMp3"), mp4 = Params::Get("KMp4"),
                        mun = Params::Get("MuN"),
                        mn1 = Params::Get("KMn1"), mn2 = Params::Get("KMn2"),
                        mn3 = Params::Get("KMn3"), mn4 = Params::Get("KMn4");
    const double tau = Q2 / (4 * MN * MN);
    const double d = 1 + Q2 / lam;
    Gep = KellyRatio(ep1, ep2, ep3, ep4, tau);
    Gen = 1.0 / (d * d) * en1 * tau / (1 + en2 * tau);
    Gmp = mup * KellyRatio(mp1, mp2, mp3, mp4, tau);
    Gmn = mun * KellyRatio(mn1, mn2, mn3, mn4, tau);
  }

  // Dirac F1 and Pauli F2 from the Sachs form factors.
  void DiracPauli(double Q2, double &F1p, double &F1n, double &F2p, double &F2n) {
    double Gep, Gen, Gmp, Gmn;
    Sachs(Q2, Gep, Gen, Gmp, Gmn);
    const double tau = Q2 / (4 * MN * MN);
    F1p = (Gep + tau * Gmp) / (1 + tau);
    F1n = (Gen + tau * Gmn) / (1 + tau);
    F2p = (Gmp - Gep) / (1 + tau);
    F2n = (Gmn - Gen) / (1 + tau);
  }

  // ---- nucleon axial form factor (z-expansion) -----------------------------
  double FA_raw(double Q2) {
    static const double tcut = Params::Get("AxTcut"), t0 = Params::Get("AxT0");
    static const double s0 = std::sqrt(tcut - t0);
    static const double a[9] = {Params::Get("AxA0"), Params::Get("AxA1"), Params::Get("AxA2"),
                                Params::Get("AxA3"), Params::Get("AxA4"), Params::Get("AxA5"),
                                Params::Get("AxA6"), Params::Get("AxA7"), Params::Get("AxA8")};
    const double s = std::sqrt(tcut + Q2);
    const double z = (s - s0) / (s + s0);
    double r = 0;                       // Horner
    for (int i = 8; i > 0; --i) { r += a[i]; r *= z; }
    return r + a[0];
  }

  // F_A(Q^2)/F_A(0): only the shape is used, the model keeps its own g_A.
  double FA_ratio(double Q2) {
    static const double inv0 = 1.0 / FA_raw(0.0);
    return FA_raw(Q2) * inv0;
  }

  // rho-dominance factor for the axial half of the contact term.
  double FrhoCT(double Q2) {
    static const double mrho = Params::Get("MRhoCT"), m2 = mrho * mrho;
    return m2 / (m2 + Q2);
  }

  // ---- N -> Delta axial: C5A(Q^2)/C5A(0) -----------------------------------
  // The model's own C5A(0) = sqrt(2/3) hD, with hD rebuilt from Gamma(Delta -> N pi).
  double C5A_model() {
    const double hDcard = Params::Opt("hD", 0.0);
    if(hDcard > 0.0) return std::sqrt(2.0 / 3.0) * hDcard;
    const double mn  = Params::Opt("MNucleon", MN);
    const double md  = Params::Opt("MDelta", MD);
    const double mpi = Params::Opt("MPiP", MPiC);
    const double fpi = Params::Get("fPi");
    const double gam = Params::Get("GamNPi");
    const double a = md * md - mn * mn - mpi * mpi;
    const double k = std::sqrt(a * a - 4.0 * mpi * mpi * mn * mn) / (2.0 * md);
    const double en = std::sqrt(k * k + mn * mn);
    const double hD = std::sqrt(12.0 * M_PI * fpi * fpi * md * gam
                                / (k * k * k * (en + mn)));
    return std::sqrt(2.0 / 3.0) * hD;
  }

  // Modified dipole.  CA5Norm > 0 overrides C5A(0) (e.g. 0.867 for the HNV fit,
  // together with MADeltaSq = 0.970225); 0 keeps the model's prediction.
  double C5A_ratio(double Q2) {
    static const double M2 = Params::Get("MADeltaSq");
    static const double norm = Params::Opt("CA5Norm", 0.0);
    static const double scale = norm > 0.0 ? norm / C5A_model() : 1.0;
    const double d = 1 + Q2 / M2;
    return scale / (d * d) / (1 + Q2 / (3 * M2));
  }

  // ---- pion poles (PCAC) ----------------------------------------------------
  // Nucleon induced pseudoscalar, F_P/F_A without the 2 MN (that is in the coupling).
  double FPP(double Q2)  { return FA_ratio(Q2) / (Q2 + MPiC * MPiC); }   // W
  double FPPZ(double Q2) { return FA_ratio(Q2) / (Q2 + MPi * MPi); }     // Z

  // Pion-pole diagram (boson -> pion -> pi N N), dimension 1/mass^2.
  double FPPpi(double Q2)  { return FrhoCT(Q2) / (Q2 + MPiC * MPiC); }   // W
  double FPPZpi(double Q2) { return FrhoCT(Q2) / (Q2 + MPi * MPi); }     // Z

  // Delta induced pseudoscalar: C6A = C5A MN^2/(Q^2 + mpi^2), MN^2 in the operator.
  double FC6A(double Q2) { return C5A_ratio(Q2) / (Q2 + MPiC * MPiC); }

  // ---- N -> Delta vector (CVC) dipole, 1 at Q^2 = 0 -------------------------
  double VND_ratio(double Q2) {
    static const double M2 = Params::Get("MV2ND");
    const double d = 1 + Q2 / M2;
    return 1.0 / (d * d) / (1 + Q2 / (4 * M2));
  }

  // ---- neutral-current combinations, J^Z = J^3 - 2 sw^2 J^EM ----------------
  // Uses sw^2 = Sw2FF, the low-Q value.
  void FZ(double Q2, bool pauli, double &Zp, double &Zn) {
    double F1p, F1n, F2p, F2n;
    DiracPauli(Q2, F1p, F1n, F2p, F2n);
    const double a = pauli ? F2p : F1p, b = pauli ? F2n : F1n;
    static const double s = Params::Get("Sw2FF");
    Zp = 0.5 * (a - b) - 2 * s * a;
    Zn = 0.5 * (b - a) - 2 * s * b;
  }

  // ---- Delta virtuality form factor on the strong pi-N-Delta vertex ---------
  // F(p^2) = Lambda^4 / (Lambda^4 + (p^2 - M_Delta^2)^2), Lambda = LambdaFF.
  // Equals 1 for an on-shell Delta.
  double DeltaLam(double p2) {
    static const double lam = Params::Get("LambdaFF"), L4 = lam * lam * lam * lam;
    static const double md = Params::Opt("MDelta", MD), md2 = md * md;
    const double d = p2 - md2;
    return L4 / (L4 + d * d);
  }

  // ---- form factor names ----------------------------------------------------
  enum class FFId {
    F1p, F1n, F2p, F2n,
    FW1, FW2, FWA,
    FZ1p, FZ1n, FZ2p, FZ2n, FZA, FZ1pi,
    FCA5, FVND, FPP, FPPZ, FPPpi, FPPZpi, FC6A, FrhoCT,
    FDeltaLam
  };

  FFId Resolve(const std::string &id) {
    static const std::map<std::string, FFId> ids = {
      {"F1p", FFId::F1p},       {"F1n", FFId::F1n},
      {"F2p", FFId::F2p},       {"F2n", FFId::F2n},
      {"FW1", FFId::FW1},       {"FW2", FFId::FW2},       {"FWA", FFId::FWA},
      {"FZ1p", FFId::FZ1p},     {"FZ1n", FFId::FZ1n},
      {"FZ2p", FFId::FZ2p},     {"FZ2n", FFId::FZ2n},
      {"FZA", FFId::FZA},       {"FZ1pi", FFId::FZ1pi},
      {"FCA5", FFId::FCA5},     {"FVND", FFId::FVND},
      {"FPP", FFId::FPP},       {"FPPZ", FFId::FPPZ},
      {"FPPpi", FFId::FPPpi},   {"FPPZpi", FFId::FPPZpi},
      {"FC6A", FFId::FC6A},     {"FrhoCT", FFId::FrhoCT},
      {"FDeltaLam", FFId::FDeltaLam},
    };
    const auto it = ids.find(id);
    if(it == ids.end())
      THROW(fatal_error, "Form factor not implemented '" + id + "'");
    return it->second;
  }

}

class FormFactorInterface : public METOOLS::Form_Factor {
  private:
    std::string m_id;
    FFId m_ff;

    // Cache of the last call (Comix repeats the same argument per helicity).
    mutable double m_arg{std::numeric_limits<double>::quiet_NaN()};
    mutable double m_val{0.0};

    double Compute(double arg) const {
        // nucleon electromagnetic
        if(m_ff == FFId::F1p || m_ff == FFId::F1n || m_ff == FFId::F2p || m_ff == FFId::F2n) {
            double F1p, F1n, F2p, F2n;
            DiracPauli(arg, F1p, F1n, F2p, F2n);
            if(m_ff == FFId::F1p) return F1p;
            if(m_ff == FFId::F1n) return F1n;
            if(m_ff == FFId::F2p) return F2p;
            return F2n;
        }

        // charged current: isovector combinations (CVC) and the axial ratio
        if(m_ff == FFId::FW1 || m_ff == FFId::FW2) {
            double F1p, F1n, F2p, F2n;
            DiracPauli(arg, F1p, F1n, F2p, F2n);
            return m_ff == FFId::FW1 ? F1p - F1n : F2p - F2n;
        }
        if(m_ff == FFId::FWA) return FA_ratio(arg);

        // neutral current.  FZ1p/FZ1n divide by the model's own point-like
        // coupling, so they also correct sin^2(theta_W) (FZ1p(0) is not 1).
        if(m_ff == FFId::FZ1p || m_ff == FFId::FZ1n) {
            double Zp, Zn; FZ(arg, false, Zp, Zn);
            return m_ff == FFId::FZ1p ? Zp / (0.5 - 2 * Params::Sw2Model())
                                      : Zn / (-0.5);
        }
        if(m_ff == FFId::FZ2p || m_ff == FFId::FZ2n) {
            double Zp, Zn; FZ(arg, true, Zp, Zn);
            return m_ff == FFId::FZ2p ? Zp : Zn;
        }
        if(m_ff == FFId::FZA) return FA_ratio(arg);
        // Z-pi-pi and Z contact: F1V with the same sin^2(theta_W) correction.
        if(m_ff == FFId::FZ1pi) {
            double F1p, F1n, F2p, F2n;
            DiracPauli(arg, F1p, F1n, F2p, F2n);
            return (F1p - F1n) * (1.0 - 2.0 * Params::Get("Sw2FF"))
                               / (1.0 - 2.0 * Params::Sw2Model());
        }

        // N -> Delta and pion poles
        if(m_ff == FFId::FCA5) return C5A_ratio(arg);
        if(m_ff == FFId::FVND) return VND_ratio(arg);
        if(m_ff == FFId::FPP) return FPP(arg);
        if(m_ff == FFId::FPPZ) return FPPZ(arg);
        if(m_ff == FFId::FPPpi) return FPPpi(arg);
        if(m_ff == FFId::FPPZpi) return FPPZpi(arg);
        if(m_ff == FFId::FC6A) return FC6A(arg);
        if(m_ff == FFId::FrhoCT) return FrhoCT(arg);

        // Delta virtuality (arg = p_Delta^2)
        if(m_ff == FFId::FDeltaLam) return DeltaLam(arg);

        THROW(fatal_error, "Form factor not implemented '" + m_id + "'");
    }

  public:
    FormFactorInterface(const Vertex_Key &key, const std::string &id)
      : Form_Factor(id, key), m_id(id), m_ff(Resolve(id)) {}

    double FF(double arg) const override {
        if(arg == m_arg) return m_val;
        m_arg = arg;
        return m_val = Compute(arg);
    }
};


// Registers one form factor name with Sherpa.
#define FORM_FACTOR_IMPLEMENTATION(NAME)                                                         \
    namespace achilles {                                                                         \
    class NAME : public Getter_Function<Form_Factor, Vertex_Key, std::less<std::string> > {      \
        static NAME s_initializer;                                                               \
                                                                                                 \
      protected:                                                                                 \
        void PrintInfo(std::ostream &str, const size_t) const {                                  \
            str << #NAME;                                                                        \
        }                                                                                        \
        Object_Type *operator()(const Parameter_Type &args) const {                              \
            return new FormFactorInterface(args, #NAME);                                         \
        }                                                                                        \
                                                                                                 \
      public:                                                                                    \
        NAME(const std::string &name)                                                            \
            : ATOOLS::Getter_Function<Form_Factor, Vertex_Key, std::less<std::string> >(name) {} \
    };                                                                                           \
    NAME NAME::s_initializer(#NAME);                                                             \
    }

FORM_FACTOR_IMPLEMENTATION(F1p)
FORM_FACTOR_IMPLEMENTATION(F1n)
FORM_FACTOR_IMPLEMENTATION(F2p)
FORM_FACTOR_IMPLEMENTATION(F2n)
FORM_FACTOR_IMPLEMENTATION(FDeltaLam)
FORM_FACTOR_IMPLEMENTATION(FW1)
FORM_FACTOR_IMPLEMENTATION(FW2)
FORM_FACTOR_IMPLEMENTATION(FWA)
FORM_FACTOR_IMPLEMENTATION(FZ1p)
FORM_FACTOR_IMPLEMENTATION(FZ1n)
FORM_FACTOR_IMPLEMENTATION(FZ2p)
FORM_FACTOR_IMPLEMENTATION(FZ2n)
FORM_FACTOR_IMPLEMENTATION(FZA)
FORM_FACTOR_IMPLEMENTATION(FZ1pi)
FORM_FACTOR_IMPLEMENTATION(FCA5)
FORM_FACTOR_IMPLEMENTATION(FVND)
FORM_FACTOR_IMPLEMENTATION(FPP)
FORM_FACTOR_IMPLEMENTATION(FPPZ)
FORM_FACTOR_IMPLEMENTATION(FPPpi)
FORM_FACTOR_IMPLEMENTATION(FPPZpi)
FORM_FACTOR_IMPLEMENTATION(FC6A)
FORM_FACTOR_IMPLEMENTATION(FrhoCT)

#undef FORM_FACTOR_IMPLEMENTATION
