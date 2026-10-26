# ChiPT: π, N and Δ(1232) model

FeynRules model and UFO for single-pion production on nucleons, from a chiral
Lagrangian with nucleons and the Δ(1232).

- γ, W and Z couple through the chiral covariant derivative, so EM, CC and NC
  processes come from one Lagrangian: `γ N → π N`, `e N → e' N π`, `ν N → ℓ N' π`.
- πNΔ and γNΔ couplings are the gauge-invariant (Pascalutsa) ones of
  Fernández-Ramírez, Moya de Guerra and Udías,
  [nucl-th/0509020](https://arxiv.org/abs/nucl-th/0509020).
- Nucleon and N → Δ form factors are UFO form factors, placed to respect CVC
  and PCAC.
- Tree level, up to O(1/f_π²).

## Contents

| Path | |
|---|---|
| `ChiPT.fr` | FeynRules model file |
| `ChiPT_UFO/`, `param_card.dat` | Exported UFO and default param card |
| `sherpa/FormFactors.C` | Sherpa plugin implementing the form factors |
| `sherpa/build_model.sh` | Builds and installs the Sherpa model library with the plugin |
| `sherpa/Sherpa.yaml` | Example run card |
| `paper/` | Scripts, run cards and data to reproduce the paper's figures ([paper/README.md](paper/README.md)) |
| `ExportUFO.nb` | Notebook that regenerates the UFO and param card |
| `feynrules/*.patch` | Patch to FeynRules 2.3.49 needed for the export |

## Running with Sherpa

The UFO only names the form factors; their values come from
`sherpa/FormFactors.C`. This needs a Sherpa with the `Form_Factor::FF(double)`
interface, UFO 2.0 form factors and spin-3/2 particles. Until these are
released, use the `570-update-achilles-interface` branch of Sherpa
([sherpa-team/sherpa!899](https://gitlab.com/sherpa-team/sherpa/-/merge_requests/899),
commit 1cd6f079f used for the paper). Its fixed-target fixes are also needed
for runs with a nucleon beam and `PDF_SET: None`, as in `paper/free_nucleon/`.

```sh
sherpa/build_model.sh                  # or: sherpa/build_model.sh /path/to/Sherpa-generate-model
```

This runs `Sherpa-generate-model` on `ChiPT_UFO/`, adds `FormFactors.C` to the
generated library in `ChiPT_UFO/.sherpa/`, then builds and installs it into
Sherpa's library directory. `NCORE` sets the number of build jobs (default 4).

- Use this repository's `param_card.dat` (`UFO_PARAM_CARD`), since
  `Sherpa-generate-model` may overwrite it.
- Form factor parameters are read from the param card, so no recompilation is
  needed to change them.
- `Form factor not implemented 'F1p'` means the installed library does not
  contain `FormFactors.C`: rerun `sherpa/build_model.sh`.

### With Achilles

Achilles runs the model through its Sherpa backend. This needs Achilles with
the Sherpa 3 interface (<Achilles version>; the paper used the `sherpa3-0`
branch at commit ff313b48).

In Achilles runs with this model, the nucleon form factors are ChiPT's, with
parameters from the param card. Achilles registers its own nucleon form
factors as `AchF1p`, `AchF1n`, `AchF2p`, `AchF2n`, so they don't clash with this
plugin's `F1p`, `F1n`, `F2p`, `F2n`. Older Achilles builds, including ff313b48,
register `F1p` … `F2n` themselves. Sherpa then warns `Doubled identifier "F1p"`
and uses the model's, which are loaded later, so the result is the same.

## Regenerating the UFO

Stock FeynRules does not work. Patch a fresh FeynRules 2.3.49:

```sh
curl -LO https://feynrules.irmp.ucl.ac.be/downloads/feynrules-current.tar.gz
# sha256 85eddac7b5b61eac772c93df93a0c44bc34f531069985a719fe454f0a54afc81
tar xzf feynrules-current.tar.gz && cd feynrules-current
patch -p1 < /path/to/ChiPT/feynrules/feynrules-2.3.49-chipt.patch
```

The patch fixes Mathematica 14 compatibility (`Commutator`, `MatrixSymbol`),
switches to the Python 3 UFO writer, and writes form factors into the Lorentz
structures (UFO 2.0) instead of the couplings.

Then set the FeynRules path in the first cell of `ExportUFO.nb` and evaluate
it. It writes `ChiPT_UFO/` and `param_card.dat` and runs basic consistency checks.

## Parameters

| Parameter | Block | Default | Meaning |
|---|---|---|---|
| `fPi` | CHIPTPION | 0.093 GeV | f_π |
| `ga` | CHIPTPION | 1.2754 | g_A |
| `MDelta`, `WDelta*` | MASS, DECAY | 1.209, 0.102 GeV | Δ pole mass and fixed width |
| `GamNPi` | CHIPTDELTA | 0.1014 GeV | Γ(Δ → Nπ), sets `hD` |
| `GamNGam` | CHIPTDELTA | 0.00057 GeV | Γ(Δ → Nγ), sets G1, G2 |
| `G2overG1` | CHIPTDELTA | 0.398 | E2/M1 mixing |
| `LambdaFF` | CHIPTFF | 1.121 GeV | Δ off-shell cutoff Λ |
| `MADeltaSq` | CHIPTFF | 1.1025 GeV² | N → Δ axial mass² |
| `CA5Norm` | CHIPTFF | 0 | C5A(0) override (0 = √(2/3) `hD`) |
| `Sw2FF` | CHIPTFF | 0.23863 | sin²θ_W for NC form factors |

Couplings are derived from these, so they follow changes to masses and widths.
Vector form factors use the Kelly fit (block `CHIPTKELLY`); the nucleon axial
form factor uses a z-expansion (block `CHIPTZEXP`, plus `t_cut`, `t_0` in `CHIPTFF`).

## Form factors

| Function | Vertex |
|---|---|
| `F1p`, `F1n`, `F2p`, `F2n` | γNN (Kelly) |
| `FW1`, `FW2`, `FWA` | WNN vector and axial |
| `FZ1p`, `FZ1n`, `FZ2p`, `FZ2n`, `FZA` | ZNN |
| `FPP`, `FPPZ`, `FPPpi`, `FPPZpi` | induced pseudoscalar / pion pole |
| `FZ1pi`, `FrhoCT` | contact and pion-in-flight terms (Ward identities) |
| `FCA5`, `FC6A`, `FVND` | N → Δ axial and vector |
| `FDeltaLam` | πNΔ: Λ⁴ / (Λ⁴ + (p² − M_Δ²)²) |

All depend on Q², except `FDeltaLam` which depends on p_Δ² and suppresses
off-shell Δ contributions.

## Limitations and plans

- Fixed Δ width; tree level only; Cabibbo mixing only (V_ud).
- Planned for v2: multi-pion production, higher resonances, strange baryons.

## License

GPL-3.0, see `LICENSE`.
