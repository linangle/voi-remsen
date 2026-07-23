# voi-remsen

**Value of remote-sensing information for harmful-algal-bloom harvest decisions.**

How much is a satellite signal worth for deciding when to close a shellfish harvest,
given that satellites can (imperfectly) inform *biomass* but the regulatory trigger is
*tissue toxicity* (domoic acid ≥ 20 µg/g)?

Design principle: **remote sensing appears only in the observation channel; the world
model is built entirely from in-situ truth.**

## Status

| layer | state |
|---|---|
| Data pipeline (calHABMAP + CDPH, coordinate-matched) | done |
| Latent ecological world model | done — factorial biomass × toxin-quota, Santa Cruz Wharf |
| Tissue-burden kinetics | done — measured, and shown *not* to need a separate state |
| Seasonality | supported but not established (LFO-corrected) |
| Remote-sensing observation channel | **not built** |
| Actions, rewards, policy, VoI | **not built** |

No current claim about sensor break-even accuracy is supported by this code.

## Headline finding

At identical biomass, tissue burden differs **9.3×** and closure probability ~40×
depending on the latent toxin quota (0.8% vs 33%). Biomass does not determine
toxicity, so a biomass-only sensor faces an irreducible ceiling on its value — this is
the mechanism behind EVPI(biomass) < EVPI(toxicity). Biomass is decisive for ruling
danger *out*, and nearly uninformative for ruling it *in*.

## Read this first

**[`docs/methodology.md`](docs/methodology.md)** — the full handoff: concepts from the
ground up, the mathematics and scientific basis for every modelling choice, validation
evidence, current results, errors found and corrected, limitations, and next steps.

## Reproduce

```bash
pip install -e ".[worldmodel]"
python scripts/00_build_panels.py            # raw -> interim panels
python scripts/03_match_radius_scan.py       # coordinate matching, 3 km vs 5 km
python scripts/01_fit_worldmodel.py 2 3 4 5  # K-scan with convergence gating
python scripts/02_fit_seasonal.py 1          # seasonal Q at the production K
```
