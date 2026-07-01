# Prediction Market Visualization

This repository contains the article and a small interactive simulator for its model.

## Run

```bash
python3 market_sim.py
```

The app uses only `numpy`, `scipy`, `matplotlib`, and Python stdlib `tkinter`. Adjust sliders, then press `Recompute now` to update the plots. Use `Load article example` to restore the paper example and the fixed article posterior `N(0.55, 0.50^2)`. It renders Matplotlib plots through the non-GUI Agg backend, so it does not require `PIL.ImageTk`.

For a non-interactive numerical check against the article example:

```bash
python3 market_sim.py --check
```

## What The Simulator Shows

The simulator starts from synthetic resolved historical markets, fits a posterior over the calibration coefficient `beta`, and then computes recommendations for two currently open CPMM markets. The recommendation stages are:

1. Raw frictionless Kelly using the bot forecast directly.
2. Raw exact-CPMM Kelly with slippage.
3. Isolated posterior Kelly after beta shrinkage.
4. Joint posterior Kelly for the two current markets with shared-beta dependence.
5. Final dynamic recommendation after reserving cash for Poisson-arriving future markets.

## Main Controls

- History controls set the synthetic sample size, random seed, prior, true beta, edge size, and wrong-side rate.
- Current-market controls set wealth, quotes, bot forecasts, and CPMM liquidity for markets A and B.
- Future controls set Poisson arrival rate, truncation, future market quote, future bot forecast, wrong-edge rate, and liquidity.

The implementation uses the finite two-stage Poisson model from the article rather than a full continuous-time Bellman PDE solver.
