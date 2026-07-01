#!/usr/bin/env python3
"""Interactive simulation for the prediction-market allocation article."""

from __future__ import annotations

import argparse
import base64
import io
import math
import tkinter as tk
from dataclasses import dataclass, replace
from tkinter import ttk

import numpy as np
import matplotlib

matplotlib.use("Agg")
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
from scipy.optimize import brentq, minimize, minimize_scalar
from scipy.special import expit, gammaln
from scipy.stats import norm


EPS = 1e-12


def clip_prob(x: float | np.ndarray) -> float | np.ndarray:
    return np.clip(x, 1e-6, 1.0 - 1e-6)


def logit(p: float | np.ndarray) -> float | np.ndarray:
    p = clip_prob(p)
    return np.log(p / (1.0 - p))


@dataclass(frozen=True)
class Market:
    name: str
    y: float
    n: float
    r: float
    p_star: float

    @property
    def q(self) -> float:
        return self.r * self.n / ((1.0 - self.r) * self.y + self.r * self.n)

    def with_quote_and_liquidity(self, q: float, liquidity: float, p_star: float) -> "Market":
        total = max(10.0, liquidity)
        q = float(clip_prob(q))
        # For r=1/2, q = n / (y+n). Keep total pool depth fixed.
        return Market(self.name, total * (1.0 - q), total * q, 0.5, float(clip_prob(p_star)))


@dataclass(frozen=True)
class Params:
    wealth: float = 300.0
    history_n: int = 80
    seed: int = 7
    prior_mu: float = 0.0
    prior_sigma: float = 1.0
    true_beta: float = 0.7
    quote_dispersion: float = 1.0
    edge_scale: float = 1.0
    wrong_side_rate: float = 0.12
    beta_grid_min: float = -2.5
    beta_grid_max: float = 3.0
    beta_grid_n: int = 321
    poisson_lambda: float = 3.0
    poisson_n_max: int = 15
    future_wrong_rate: float = 0.0
    current_a_q: float = 0.35
    current_a_p: float = 0.70
    current_a_liq: float = 1000.0
    current_b_q: float = 0.25
    current_b_p: float = 0.60
    current_b_liq: float = 500.0
    future_q: float = 0.50
    future_p: float = 0.75
    future_liq: float = 800.0


@dataclass
class History:
    q: np.ndarray
    p_star: np.ndarray
    y: np.ndarray
    p_true: np.ndarray


@dataclass
class Posterior:
    beta: np.ndarray
    weights: np.ndarray

    @property
    def mean(self) -> float:
        return float(np.sum(self.beta * self.weights))

    @property
    def std(self) -> float:
        return float(np.sqrt(np.sum(((self.beta - self.mean) ** 2) * self.weights)))

    def prob_lt(self, x: float) -> float:
        return float(np.sum(self.weights[self.beta < x]))

    def prob_between(self, lo: float, hi: float) -> float:
        return float(np.sum(self.weights[(self.beta >= lo) & (self.beta <= hi)]))


@dataclass
class StrategyResult:
    label: str
    shares: np.ndarray
    costs: np.ndarray


@dataclass
class SimulationResult:
    params: Params
    history: History
    posterior: Posterior
    markets: tuple[Market, Market, Market]
    pbar: np.ndarray
    joint_current: np.ndarray
    corr_ab: float
    stages: list[StrategyResult]
    future_policy: list[tuple[int, float, float, float]]


def p_beta(q: float | np.ndarray, p_star: float | np.ndarray, beta: np.ndarray) -> np.ndarray:
    return expit(logit(q) + np.outer(beta, np.ravel(logit(p_star) - logit(q)))).squeeze()


def generate_history(params: Params) -> History:
    rng = np.random.default_rng(params.seed)
    logits_q = rng.normal(0.0, params.quote_dispersion, params.history_n)
    q = clip_prob(expit(logits_q))
    edge = rng.normal(0.0, params.edge_scale, params.history_n)
    wrong = rng.random(params.history_n) < params.wrong_side_rate
    edge[wrong] *= -1.0
    p_star = clip_prob(expit(logit(q) + edge))
    p_true = clip_prob(expit(logit(q) + params.true_beta * (logit(p_star) - logit(q))))
    y = rng.binomial(1, p_true)
    return History(q=q, p_star=p_star, y=y.astype(float), p_true=p_true)


def fit_posterior(params: Params, history: History) -> Posterior:
    beta = np.linspace(params.beta_grid_min, params.beta_grid_max, params.beta_grid_n)
    prior_sigma = max(0.05, params.prior_sigma)
    logw = -0.5 * ((beta - params.prior_mu) / prior_sigma) ** 2
    d = logit(history.p_star) - logit(history.q)
    probs = clip_prob(expit(logit(history.q)[None, :] + beta[:, None] * d[None, :]))
    logw += np.sum(history.y[None, :] * np.log(probs) + (1.0 - history.y[None, :]) * np.log1p(-probs), axis=1)
    logw -= np.max(logw)
    weights = np.exp(logw)
    weights /= np.sum(weights)
    return Posterior(beta=beta, weights=weights)


def normal_posterior(mu: float, sigma: float, n: int = 501) -> Posterior:
    beta = np.linspace(mu - 5.0 * sigma, mu + 5.0 * sigma, n)
    weights = norm.pdf(beta, mu, sigma)
    weights /= np.sum(weights)
    return Posterior(beta=beta, weights=weights)


def predictive_probs(markets: list[Market] | tuple[Market, ...], posterior: Posterior) -> np.ndarray:
    out = []
    for market in markets:
        p = expit(logit(market.q) + posterior.beta * (logit(market.p_star) - logit(market.q)))
        out.append(float(np.sum(p * posterior.weights)))
    return np.array(out)


def joint_probs(markets: list[Market] | tuple[Market, ...], posterior: Posterior) -> np.ndarray:
    m = len(markets)
    probs_by_beta = []
    for market in markets:
        probs_by_beta.append(clip_prob(expit(logit(market.q) + posterior.beta * (logit(market.p_star) - logit(market.q)))))
    rows = []
    for mask in range(2**m):
        prob_beta = np.ones_like(posterior.beta)
        for i in range(m):
            yi = (mask >> (m - i - 1)) & 1
            p = probs_by_beta[i]
            prob_beta *= p if yi else (1.0 - p)
        rows.append(float(np.sum(prob_beta * posterior.weights)))
    return np.array(rows)


def cpmm_cost(market: Market, shares: float) -> float:
    shares = max(0.0, float(shares))
    if shares <= 0.0:
        return 0.0
    y, n, r = market.y, market.n, market.r
    if abs(r - 0.5) < 1e-12:
        bcoef = y + n - shares
        disc = bcoef * bcoef + 4.0 * n * shares
        return max(0.0, (-bcoef + math.sqrt(disc)) / 2.0)

    invariant_log = r * math.log(y) + (1.0 - r) * math.log(n)

    def f(b: float) -> float:
        yp = y + b - shares
        npool = n + b
        if yp <= 0.0 or npool <= 0.0:
            return -1e100
        return r * math.log(yp) + (1.0 - r) * math.log(npool) - invariant_log

    lo = max(0.0, shares - y + 1e-10)
    hi = max(1.0, shares, y + n)
    while f(hi) <= 0.0:
        hi *= 2.0
    return float(brentq(f, lo, hi, maxiter=100))


def cpmm_marginal_after(market: Market, shares: float) -> float:
    b = cpmm_cost(market, shares)
    yp = market.y + b - shares
    npool = market.n + b
    return market.r * npool / ((1.0 - market.r) * yp + market.r * npool)


def max_affordable_shares(market: Market, cash: float) -> float:
    if cash <= 0.0:
        return 0.0
    y, n = market.y, market.n
    return max(0.0, y + cash - (y * n) / (n + cash))


def isolated_exact(market: Market, p: float, wealth: float) -> tuple[float, float]:
    p = float(clip_prob(p))
    if p <= market.q:
        return 0.0, 0.0
    upper = 0.999 * max_affordable_shares(market, wealth * 0.999)

    def objective(s: float) -> float:
        c = cpmm_cost(market, s)
        w0 = wealth - c
        w1 = wealth - c + s
        if w0 <= 0.0 or w1 <= 0.0:
            return 1e30
        return -((1.0 - p) * math.log(w0) + p * math.log(w1))

    res = minimize_scalar(objective, bounds=(0.0, upper), method="bounded", options={"xatol": 1e-7})
    s = max(0.0, float(res.x))
    return s, cpmm_cost(market, s)


def frictionless_kelly(market: Market, p: float, wealth: float) -> tuple[float, float]:
    if p <= market.q:
        return 0.0, 0.0
    fraction = min(0.99, max(0.0, (p - market.q) / (1.0 - market.q)))
    cost = fraction * wealth
    return cost / market.q, cost


def optimize_static(markets: tuple[Market, Market], outcome_probs: np.ndarray, wealth: float) -> tuple[np.ndarray, np.ndarray]:
    def objective(svec: np.ndarray) -> float:
        s0, s1 = np.maximum(svec, 0.0)
        c0, c1 = cpmm_cost(markets[0], s0), cpmm_cost(markets[1], s1)
        wealths = np.array([
            wealth - c0 - c1,
            wealth - c0 - c1 + s1,
            wealth - c0 - c1 + s0,
            wealth - c0 - c1 + s0 + s1,
        ])
        if np.any(wealths <= 0.0):
            return 1e20 + 1e10 * float(np.sum(np.minimum(wealths, 0.0) ** 2))
        return -float(np.sum(outcome_probs * np.log(wealths)))

    bnds = [
        (0.0, 0.999 * max_affordable_shares(markets[0], wealth * 0.999)),
        (0.0, 0.999 * max_affordable_shares(markets[1], wealth * 0.999)),
    ]
    x0 = np.array([0.2 * bnds[0][1], 0.2 * bnds[1][1]])
    res = minimize(objective, x0=x0, method="Nelder-Mead", options={"maxiter": 700, "xatol": 1e-6, "fatol": 1e-8})
    x = np.clip(res.x, [b[0] for b in bnds], [b[1] for b in bnds])
    res2 = minimize(objective, x0=x, method="L-BFGS-B", bounds=bnds, options={"maxiter": 300, "ftol": 1e-10})
    if res2.success or res2.fun < res.fun:
        x = res2.x
    x = np.maximum(x, 0.0)
    return x, np.array([cpmm_cost(markets[0], x[0]), cpmm_cost(markets[1], x[1])])


def second_stage_value(
    current_markets: tuple[Market, Market],
    future_market: Market,
    posterior: Posterior,
    s_ab: np.ndarray,
    n_arrivals: int,
    wealth: float,
) -> tuple[float, float, float, float]:
    ca = cpmm_cost(current_markets[0], s_ab[0])
    cb = cpmm_cost(current_markets[1], s_ab[1])
    base_cash = wealth - ca - cb
    if base_cash <= 0.0:
        return -1e30, 0.0, 0.0, base_cash
    if n_arrivals <= 0:
        jp = joint_probs(current_markets, posterior)
        ws = np.array([base_cash, base_cash + s_ab[1], base_cash + s_ab[0], base_cash + s_ab[0] + s_ab[1]])
        if np.any(ws <= 0.0):
            return -1e30, 0.0, 0.0, base_cash
        return float(np.sum(jp * np.log(ws))), 0.0, 0.0, base_cash

    p_a = clip_prob(expit(logit(current_markets[0].q) + posterior.beta * (logit(current_markets[0].p_star) - logit(current_markets[0].q))))
    p_b = clip_prob(expit(logit(current_markets[1].q) + posterior.beta * (logit(current_markets[1].p_star) - logit(current_markets[1].q))))
    p_f = clip_prob(expit(logit(future_market.q) + posterior.beta * (logit(future_market.p_star) - logit(future_market.q))))
    max_u = 0.999 * max_affordable_shares(future_market, base_cash * 0.999 / n_arrivals)

    binom_coeff = np.array([
        math.exp(gammaln(n_arrivals + 1) - gammaln(m + 1) - gammaln(n_arrivals - m + 1))
        for m in range(n_arrivals + 1)
    ])

    def value_for_u(u: float) -> float:
        cf = cpmm_cost(future_market, u)
        cash = base_cash - n_arrivals * cf
        if cash <= 0.0:
            return -1e30
        total_by_beta = np.zeros_like(posterior.beta)
        for a in (0, 1):
            pra = p_a if a else 1.0 - p_a
            for b in (0, 1):
                prab = pra * (p_b if b else 1.0 - p_b)
                for m in range(n_arrivals + 1):
                    wt = cash + a * s_ab[0] + b * s_ab[1] + m * u
                    if wt <= 0.0:
                        return -1e30
                    pm = binom_coeff[m] * (p_f ** m) * ((1.0 - p_f) ** (n_arrivals - m))
                    total_by_beta += prab * pm * math.log(wt)
        return float(np.sum(posterior.weights * total_by_beta))

    opt = minimize_scalar(lambda u: -value_for_u(u), bounds=(0.0, max_u), method="bounded", options={"xatol": 1e-5})
    u = max(0.0, float(opt.x))
    cf = cpmm_cost(future_market, u)
    cash_left = base_cash - n_arrivals * cf
    return value_for_u(u), u, cf, cash_left


def optimize_dynamic(
    current_markets: tuple[Market, Market],
    future_market: Market,
    posterior: Posterior,
    params: Params,
) -> tuple[np.ndarray, np.ndarray, list[tuple[int, float, float, float]]]:
    lam = max(0.0, params.poisson_lambda)
    n_max = max(0, int(params.poisson_n_max))
    ns = np.arange(n_max + 1)
    logp = -lam + ns * math.log(lam if lam > 0 else 1.0) - gammaln(ns + 1)
    pois = np.exp(logp)
    if lam == 0.0:
        pois[:] = 0.0
        pois[0] = 1.0
    pois /= np.sum(pois)

    cache: dict[tuple[float, float, int], tuple[float, float, float, float]] = {}

    def value(svec: np.ndarray) -> float:
        s = np.maximum(svec, 0.0)
        key_base = (round(float(s[0]), 5), round(float(s[1]), 5))
        total = 0.0
        for n, pr in zip(ns, pois):
            key = (key_base[0], key_base[1], int(n))
            if key not in cache:
                cache[key] = second_stage_value(current_markets, future_market, posterior, s, int(n), params.wealth)
            total += pr * cache[key][0]
        return float(total)

    bnds = [
        (0.0, 0.999 * max_affordable_shares(current_markets[0], params.wealth * 0.999)),
        (0.0, 0.999 * max_affordable_shares(current_markets[1], params.wealth * 0.999)),
    ]
    jp = joint_probs(current_markets, posterior)
    static_s, _ = optimize_static(current_markets, jp, params.wealth)
    res = minimize(lambda x: -value(x), x0=static_s * 0.85, method="Nelder-Mead", options={"maxiter": 420, "xatol": 1e-5, "fatol": 1e-7})
    x = np.clip(res.x, [0.0, 0.0], [bnds[0][1], bnds[1][1]])
    res2 = minimize(lambda x: -value(x), x0=x, method="L-BFGS-B", bounds=bnds, options={"maxiter": 160, "ftol": 1e-9})
    if res2.success or -res2.fun > value(x):
        x = res2.x
    x = np.maximum(x, 0.0)
    costs = np.array([cpmm_cost(current_markets[0], x[0]), cpmm_cost(current_markets[1], x[1])])
    policy = []
    for n in range(1, min(n_max, 8) + 1):
        _, u, cf, cash = second_stage_value(current_markets, future_market, posterior, x, n, params.wealth)
        policy.append((n, u, cf, cash))
    return x, costs, policy



def run_simulation(params: Params, article_posterior: bool = False) -> SimulationResult:
    hist = generate_history(params)
    posterior = normal_posterior(0.55, 0.50) if article_posterior else fit_posterior(params, hist)
    future_d = logit(params.future_p) - logit(params.future_q)
    future_effective_p = float(clip_prob(expit(logit(params.future_q) + (1.0 - 2.0 * params.future_wrong_rate) * future_d)))
    markets = (
        Market("A", 650, 350, 0.5, 0.70).with_quote_and_liquidity(params.current_a_q, params.current_a_liq, params.current_a_p),
        Market("B", 375, 125, 0.5, 0.60).with_quote_and_liquidity(params.current_b_q, params.current_b_liq, params.current_b_p),
        Market("F", 400, 400, 0.5, 0.75).with_quote_and_liquidity(params.future_q, params.future_liq, future_effective_p),
    )
    pbar = predictive_probs(markets, posterior)
    joint = joint_probs(markets[:2], posterior)
    corr = (joint[3] - pbar[0] * pbar[1]) / math.sqrt(pbar[0] * (1.0 - pbar[0]) * pbar[1] * (1.0 - pbar[1]))

    fk = [frictionless_kelly(markets[i], markets[i].p_star, params.wealth) for i in range(2)]
    raw = [isolated_exact(markets[i], markets[i].p_star, params.wealth) for i in range(2)]
    post = [isolated_exact(markets[i], pbar[i], params.wealth) for i in range(2)]
    static_s, static_c = optimize_static(markets[:2], joint, params.wealth)
    dyn_s, dyn_c, future_policy = optimize_dynamic(markets[:2], markets[2], posterior, params)
    stages = [
        StrategyResult("Raw frictionless Kelly", np.array([fk[0][0], fk[1][0]]), np.array([fk[0][1], fk[1][1]])),
        StrategyResult("Raw exact CPMM", np.array([raw[0][0], raw[1][0]]), np.array([raw[0][1], raw[1][1]])),
        StrategyResult("Posterior isolated", np.array([post[0][0], post[1][0]]), np.array([post[0][1], post[1][1]])),
        StrategyResult("Joint posterior", static_s, static_c),
        StrategyResult("Final dynamic", dyn_s, dyn_c),
    ]
    return SimulationResult(params, hist, posterior, markets, pbar, joint, corr, stages, future_policy)


def weighted_quantile(
    values: np.ndarray,
    weights: np.ndarray,
    quantiles: float | list[float] | np.ndarray,
) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    sorter = np.argsort(values)
    sorted_values = values[sorter]
    sorted_weights = weights[sorter]
    total = np.sum(sorted_weights)
    if total <= 0.0:
        raise ValueError("weights must sum to a positive value")
    sorted_weights = sorted_weights / total
    cumulative = np.cumsum(sorted_weights)
    q = np.asarray(quantiles, dtype=float)
    return np.interp(q, cumulative, sorted_values)


def stage_by_label(result: SimulationResult, label: str) -> StrategyResult:
    for stage in result.stages:
        if stage.label == label:
            return stage
    available = ", ".join(stage.label for stage in result.stages)
    raise KeyError(f"No strategy stage labelled {label!r}. Available stages: {available}")


def independent_joint_probs(p_a: float, p_b: float) -> np.ndarray:
    return np.array([
        (1.0 - p_a) * (1.0 - p_b),
        (1.0 - p_a) * p_b,
        p_a * (1.0 - p_b),
        p_a * p_b,
    ])


def _set_axis_grid(axis) -> None:
    axis.grid(True, color="#d9dee7", linewidth=0.7, alpha=0.75)
    axis.set_axisbelow(True)


def _safe_reduction(static_share: float, dynamic_share: float) -> float:
    if abs(static_share) < EPS:
        return 0.0
    return 100.0 * (static_share - dynamic_share) / static_share


def draw_dashboard(fig: Figure, result: SimulationResult) -> None:
    fig.clear()
    market_colors = {"A": "#2364aa", "B": "#b23a48", "F": "#2a9d8f"}
    static_color = "#6d597a"
    dynamic_color = "#e76f51"
    cash_color = "#d8dee9"

    gs = fig.add_gridspec(
        3,
        3,
        height_ratios=[1.0, 1.0, 0.22],
        hspace=0.50,
        wspace=0.34,
    )
    ax_beta = fig.add_subplot(gs[0, 0])
    ax_shrink = fig.add_subplot(gs[0, 1])
    ax_dependence = fig.add_subplot(gs[0, 2])
    ax_cpmm = fig.add_subplot(gs[1, 0])
    ax_budget = fig.add_subplot(gs[1, 1])
    ax_policy = fig.add_subplot(gs[1, 2])
    ax_metrics = fig.add_subplot(gs[2, :])
    ax_metrics.axis("off")

    # 1. Calibration posterior.
    beta = result.posterior.beta
    density = result.posterior.weights / max(EPS, beta[1] - beta[0])
    ymax = max(float(np.max(density)), EPS)
    ax_beta.axvspan(beta[0], 0.0, color="#e76f51", alpha=0.15)
    ax_beta.axvspan(0.0, 1.0, color="#2a9d8f", alpha=0.15)
    ax_beta.axvspan(1.0, beta[-1], color="#f4a261", alpha=0.15)
    ax_beta.plot(beta, density, color="#244c9a", linewidth=2.0)
    ax_beta.axvline(0.0, color="#555555", linewidth=1.0)
    ax_beta.axvline(1.0, color="#555555", linewidth=1.0, linestyle="--")
    ax_beta.axvline(result.posterior.mean, color="#111111", linewidth=1.8)
    ci_lo, ci_hi = weighted_quantile(beta, result.posterior.weights, [0.05, 0.95])
    ax_beta.axvline(ci_lo, color="#111111", linewidth=1.0, linestyle=":")
    ax_beta.axvline(ci_hi, color="#111111", linewidth=1.0, linestyle=":")
    ci_y = 0.055 * ymax
    ax_beta.hlines(ci_y, ci_lo, ci_hi, color="#111111", linewidth=2.0)
    ax_beta.text((ci_lo + ci_hi) / 2.0, ci_y + 0.035 * ymax, "90% CI", ha="center", va="bottom", fontsize=8)
    y_region = 1.0
    # ax_beta.text(max(beta[0] + 0.25, -0.75), y_region, "β < 0", ha="center", fontsize=8.5)
    # ax_beta.text(0.62, y_region, "0 < β < 1", ha="center", fontsize=8.5)
    # ax_beta.text(min(beta[-1] - 0.30, 1.65), y_region, "β > 1", ha="center", fontsize=8.5)
    beta_text = (
        f"mean  {result.posterior.mean:.3f}\n"
        f"sd    {result.posterior.std:.3f}\n"
        f"90% CI [{ci_lo:.3f}, {ci_hi:.3f}]\n"
        f"P(β<0) {result.posterior.prob_lt(0):.3f}\n"
        f"P(β>1) {1.0 - result.posterior.prob_lt(1):.3f}"
    )
    ax_beta.text(
        0.98,
        0.95,
        beta_text,
        transform=ax_beta.transAxes,
        ha="right",
        va="top",
        fontsize=8.2,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "0.7", "alpha": 0.92},
    )
    ax_beta.set_title("Calibration posterior", fontsize=11.5)
    ax_beta.set_xlabel("calibration slope β", fontsize=9.5)
    ax_beta.set_ylabel("posterior density", fontsize=9.5)
    ax_beta.set_ylim(0.0, max(ymax * 1.12, 1.15))
    ax_beta.tick_params(labelsize=8.5)
    _set_axis_grid(ax_beta)

    # 2. Forecast shrinkage.
    ax_shrink.set_title("Quote, posterior belief, and forecast", fontsize=11.5)
    y_positions = np.arange(3)
    q_vals = [m.q for m in result.markets]
    pstar_vals = [m.p_star for m in result.markets]
    marker_specs = [("q", "o", q_vals), ("pbar", "D", result.pbar), ("p*", "^", pstar_vals)]
    for i, market in enumerate(result.markets):
        vals = [market.q, result.pbar[i], market.p_star]
        ax_shrink.hlines(i, min(vals), max(vals), color="#9aa3ad", linewidth=1.2)
    for label, marker, vals in marker_specs:
        ax_shrink.scatter(vals, y_positions, marker=marker, s=60, label=label, zorder=3)
    for i, market in enumerate(result.markets):
        q, pbar, pstar = market.q, float(result.pbar[i]), market.p_star
        q_x_offset = -10 if abs(q - pbar) < 0.08 else 0
        pbar_x_offset = 10 if abs(q - pbar) < 0.08 else 0
        pstar_x_offset = 10 if abs(pstar - pbar) < 0.08 else 0
        ax_shrink.annotate(f"{q:.3f}", (q, i), xytext=(q_x_offset, 10), textcoords="offset points", ha="center", va="bottom", fontsize=8.2)
        ax_shrink.annotate(f"{pbar:.3f}", (pbar, i), xytext=(pbar_x_offset, -12), textcoords="offset points", ha="center", va="top", fontsize=8.2)
        ax_shrink.annotate(f"{pstar:.3f}", (pstar, i), xytext=(pstar_x_offset, 10), textcoords="offset points", ha="center", va="bottom", fontsize=8.2)
    ax_shrink.set_xlim(0.0, 1.0)
    ax_shrink.set_yticks(y_positions, [m.name for m in result.markets])
    ax_shrink.set_ylim(2.65, -0.65)
    ax_shrink.set_xlabel("probability", fontsize=9.5)
    ax_shrink.tick_params(labelsize=8.5)
    ax_shrink.legend(loc="lower left", bbox_to_anchor=(0.02, 0.02), ncol=1, fontsize=8.5, frameon=False)
    _set_axis_grid(ax_shrink)

    # 3. Dependence from common calibration risk.
    ax_dependence.set_title(f"Dependence induced by shared β\nCorr(A,B) = {result.corr_ab:.3f}", fontsize=11.0)
    outcome_labels = ["A=0\nB=0", "A=0\nB=1", "A=1\nB=0", "A=1\nB=1"]
    actual = result.joint_current
    independent = independent_joint_probs(result.pbar[0], result.pbar[1])
    x = np.arange(len(outcome_labels))
    width = 0.36
    ax_dependence.bar(x - width / 2, actual, width, label="shared β", color="#2364aa")
    ax_dependence.bar(x + width / 2, independent, width, label="independent", color="#b8c2cc")
    for i, delta in enumerate(actual - independent):
        y = max(actual[i], independent[i]) + 0.014
        ax_dependence.text(i, y, f"Δ={delta:+.4f}", ha="center", va="bottom", fontsize=8.5)
    ax_dependence.set_xticks(x, outcome_labels)
    ax_dependence.set_ylim(0.0, max(0.35, float(np.max([actual, independent])) * 1.25))
    ax_dependence.set_ylabel("probability", fontsize=9.5)
    ax_dependence.set_xlabel("Outcome (A, B)", fontsize=9.5)
    ax_dependence.tick_params(labelsize=8.5)
    ax_dependence.legend(fontsize=8.5, frameon=False, loc="upper left")
    _set_axis_grid(ax_dependence)

    # 4. CPMM price impact.
    ax_cpmm.set_title("CPMM price impact", fontsize=11.5)
    static = stage_by_label(result, "Joint posterior")
    dynamic = stage_by_label(result, "Final dynamic")
    cpmm_xmax = 10.0
    for idx, market in enumerate(result.markets[:2]):
        color = market_colors[market.name]
        largest_stage_share = max(float(stage.shares[idx]) for stage in result.stages)
        s_max = max(10.0, 1.15 * largest_stage_share)
        cpmm_xmax = max(cpmm_xmax, s_max)
        share_grid = np.linspace(0.0, s_max, 250)
        marginal = np.array([cpmm_marginal_after(market, s) for s in share_grid])
        ax_cpmm.plot(share_grid, marginal, color=color, linewidth=2.0)
        ax_cpmm.axhline(result.pbar[idx], color=color, linestyle=":", linewidth=1.2, alpha=0.75)
        ax_cpmm.text(s_max, marginal[-1], market.name, color=color, ha="left", va="center", fontsize=9, fontweight="bold")
        ax_cpmm.text(s_max, result.pbar[idx], f"p̄{market.name}", color=color, ha="left", va="center", fontsize=8.5)
        for stage, marker, stage_color, strategy_name in ((static, "o", static_color, "static"), (dynamic, "s", dynamic_color, "dynamic")):
            share = float(stage.shares[idx])
            price = cpmm_marginal_after(market, share)
            ax_cpmm.scatter(share, price, marker=marker, s=62, color=stage_color, edgecolor="white", linewidth=0.7, zorder=4)
            offsets = {
                ("A", "static"): (-16, 12),
                ("A", "dynamic"): (14, 12),
                ("B", "static"): (-16, -15),
                ("B", "dynamic"): (14, -15),
            }
            ax_cpmm.annotate(
                f"{share:.1f}",
                xy=(share, price),
                xytext=offsets[(market.name, strategy_name)],
                textcoords="offset points",
                ha="center",
                va="center",
                fontsize=8.2,
                arrowprops={"arrowstyle": "-", "lw": 0.6},
            )
    marker_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=static_color, label="static", markersize=7),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=dynamic_color, label="dynamic", markersize=7),
    ]
    ax_cpmm.legend(handles=marker_handles, fontsize=8.5, frameon=False, loc="upper left")
    ax_cpmm.set_xlim(0.0, cpmm_xmax * 1.10)
    ax_cpmm.set_xlabel("YES shares purchased", fontsize=9.5)
    ax_cpmm.set_ylabel("post-trade marginal probability", fontsize=9.5)
    ax_cpmm.set_ylim(0.0, 1.0)
    ax_cpmm.tick_params(labelsize=8.5)
    _set_axis_grid(ax_cpmm)

    # 5. Static versus dynamic current allocation.
    ax_budget.set_title("Current allocation and retained cash", fontsize=11.5)
    budget_rows = [("Static", static), ("Dynamic", dynamic)]
    for row_idx, (_row_label, stage) in enumerate(budget_rows):
        y = row_idx
        a_cost, b_cost = float(stage.costs[0]), float(stage.costs[1])
        cash_remaining = result.params.wealth - a_cost - b_cost
        ax_budget.barh(y, a_cost, color=market_colors["A"], height=0.42)
        ax_budget.barh(y, b_cost, left=a_cost, color=market_colors["B"], height=0.42)
        ax_budget.barh(y, cash_remaining, left=a_cost + b_cost, color=cash_color, height=0.42)
        segments = [("A", a_cost, 0.0, stage.shares[0]), ("B", b_cost, a_cost, stage.shares[1])]
        for label, cost, left, shares in segments:
            center = left + cost / 2.0
            ax_budget.text(center, y, f"{label}\nM{cost:.1f}", ha="center", va="center", fontsize=8.2, color="white")
        if cash_remaining >= 25.0:
            ax_budget.text(a_cost + b_cost + cash_remaining / 2.0, y, f"cash\nM {cash_remaining:.1f}", ha="center", va="center", fontsize=8.2, color="#222222")
        else:
            ax_budget.text(min(result.params.wealth, a_cost + b_cost + cash_remaining + 3.0), y, f"cash M {cash_remaining:.1f}", ha="left", va="center", fontsize=8)
    static_cash = result.params.wealth - float(np.sum(static.costs))
    dynamic_cash = result.params.wealth - float(np.sum(dynamic.costs))
    ax_budget.text(
        0.5,
        -0.23,
        f"Dynamic policy retains M {dynamic_cash - static_cash:.2f} more",
        transform=ax_budget.transAxes,
        ha="center",
        va="top",
        fontsize=9,
        fontweight="bold",
    )
    ax_budget.set_yticks([0, 1], ["Static", "Dynamic"])
    ax_budget.set_ylim(1.75, -0.75)
    ax_budget.set_xlim(0.0, result.params.wealth)
    ax_budget.set_xlabel("cash allocation (M)", fontsize=9.5)
    ax_budget.tick_params(labelsize=8.5)
    budget_handles = [Patch(color=market_colors["A"], label="A cost"), Patch(color=market_colors["B"], label="B cost"), Patch(color=cash_color, label="retained cash")]
    ax_budget.legend(handles=budget_handles, fontsize=8.5, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 0.02), ncol=3)
    _set_axis_grid(ax_budget)

    # 6. Future allocation policy.
    ax_policy.set_title("Future allocation policy", fontsize=11.5)
    ax_policy.set_xlabel("number of future markets n", fontsize=9.5)
    ax_policy.set_ylabel("shares per F market", fontsize=9.5, color=market_colors["F"])
    ax_policy.tick_params(axis="y", labelcolor=market_colors["F"], labelsize=8.5)
    ax_policy.tick_params(axis="x", labelsize=8.5)
    ax_policy_right = ax_policy.twinx()
    ax_policy_right.set_ylabel("cash (M)", fontsize=9.5, color="#333333")
    ax_policy_right.tick_params(axis="y", labelcolor="#333333", labelsize=8.5)
    if result.future_policy:
        ns = np.array([row[0] for row in result.future_policy])
        shares = np.array([row[1] for row in result.future_policy])
        costs = np.array([row[2] for row in result.future_policy])
        cash = np.array([row[3] for row in result.future_policy])
        left_line = ax_policy.plot(ns, shares, marker="o", color=market_colors["F"], linewidth=2.0, label="shares/F")[0]
        cost_line = ax_policy_right.plot(ns, costs, marker="s", linestyle="--", color="#6d597a", linewidth=1.8, label="cost/F")[0]
        cash_line = ax_policy_right.plot(ns, cash, marker="^", linestyle=":", color="#264653", alpha=0.75, linewidth=1.4, label="cash remaining")[0]
        ax_policy.xaxis.set_major_locator(MaxNLocator(integer=True))
        # ax_policy.annotate("shares/F", xy=(ns[-1], shares[-1]), xytext=(0, 0), textcoords="offset points", va="center",
        #                    fontsize=8.5, color=market_colors["F"])
        # ax_policy_right.annotate("cost/F", xy=(ns[-1], costs[-1]), xytext=(6, 0), textcoords="offset points", va="center",
        #                          fontsize=8.5, color="#6d597a")
        # ax_policy_right.annotate("cash", xy=(ns[-1], cash[-1]), xytext=(6, 0), textcoords="offset points", va="center",
        #                          fontsize=8.5, color="#264653")
        ax_policy.legend([left_line, cost_line, cash_line], ["shares/F", "cost/F", "cash remaining"], fontsize=8.5,
                         frameon=False, loc="lower center", bbox_to_anchor=(0.5, 0.85), ncol=3)

        left_ymin, left_ymax = ax_policy.get_ylim()
        right_ymin, right_ymax = ax_policy_right.get_ylim()

        ax_policy.set_ylim(left_ymin, left_ymax + 0.15 * (left_ymax - left_ymin))
        ax_policy_right.set_ylim(right_ymin, right_ymax + 0.15 * (right_ymax - right_ymin))
    else:
        ax_policy.text(0.5, 0.5, "No future-arrival policy to display", transform=ax_policy.transAxes, ha="center", va="center", fontsize=10)
    _set_axis_grid(ax_policy)
    ax_policy_right.grid(False)

    # Bottom metrics strip.
    static_commit = float(np.sum(static.costs))
    dynamic_commit = float(np.sum(dynamic.costs))
    static_cash = result.params.wealth - static_commit
    dynamic_cash = result.params.wealth - dynamic_commit
    reduction_a = _safe_reduction(float(static.shares[0]), float(dynamic.shares[0]))
    reduction_b = _safe_reduction(float(static.shares[1]), float(dynamic.shares[1]))
    ax_metrics.text(0.01, 0.78, "POSTERIOR", transform=ax_metrics.transAxes, fontsize=9.5, fontweight="bold", va="center")
    ax_metrics.text(0.01, 0.45, f"β = {result.posterior.mean:.3f} ± {result.posterior.std:.3f}", transform=ax_metrics.transAxes, fontsize=9.5, va="center")
    ax_metrics.text(0.01, 0.16, f"p̄A = {result.pbar[0]:.3f}   p̄B = {result.pbar[1]:.3f}   p̄F = {result.pbar[2]:.3f}", transform=ax_metrics.transAxes, fontsize=9.5, va="center")
    ax_metrics.text(0.35, 0.78, "CURRENT ALLOCATION", transform=ax_metrics.transAxes, fontsize=9.5, fontweight="bold", va="center")
    ax_metrics.text(0.35, 0.45, f"Static: M {static_commit:.2f} committed, M {static_cash:.2f} retained", transform=ax_metrics.transAxes, fontsize=9.5, va="center")
    ax_metrics.text(0.35, 0.16, f"Dynamic: M {dynamic_commit:.2f} committed, M {dynamic_cash:.2f} retained", transform=ax_metrics.transAxes, fontsize=9.5, va="center")
    ax_metrics.text(0.70, 0.78, "OPTION VALUE", transform=ax_metrics.transAxes, fontsize=9.5, fontweight="bold", va="center")
    ax_metrics.text(0.70, 0.45, f"Extra retained cash: M {dynamic_cash - static_cash:.2f}", transform=ax_metrics.transAxes, fontsize=9.5, va="center")
    ax_metrics.text(0.70, 0.16, f"A reduction: {reduction_a:.2f}%   B reduction: {reduction_b:.2f}%", transform=ax_metrics.transAxes, fontsize=9.5, va="center")
    fig.subplots_adjust(left=0.055, right=0.970, top=0.930, bottom=0.080)


class MarketApp:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("Prediction Market Kelly / Calibration / Poisson Simulator")
        self.root.minsize(1250, 820)
        self.params = Params()
        self.vars: dict[str, tk.DoubleVar] = {}
        self.value_labels: dict[str, ttk.Label] = {}
        self.use_article_posterior = False
        self.status = tk.StringVar(value="Ready")
        self._build()
        self.refresh()

    def _add_slider(self, parent: ttk.Frame, label: str, name: str, lo: float, hi: float, resolution: float) -> None:
        row = ttk.Frame(parent)
        row.pack(fill="x", padx=4, pady=0)
        ttk.Label(row, text=label, width=14).pack(side="left")
        var = tk.DoubleVar(value=float(getattr(self.params, name)))
        self.vars[name] = var
        value_label = ttk.Label(row, text=self._format_slider_value(name, var.get()), width=7, anchor="e")
        self.value_labels[name] = value_label
        scale = tk.Scale(
            row,
            from_=lo,
            to=hi,
            resolution=resolution,
            orient="horizontal",
            variable=var,
            length=130,
            showvalue=False,
            width=9,
            sliderlength=18,
            highlightthickness=0,
            command=lambda _value, slider_name=name: self._slider_changed(slider_name),
        )
        scale.pack(side="left", fill="x", expand=True)
        value_label.pack(side="left", padx=(3, 0))

    def _format_slider_value(self, name: str, value: float) -> str:
        if name in {"history_n", "seed", "poisson_n_max"}:
            return str(int(round(value)))
        if abs(value) >= 100:
            return f"{value:.0f}"
        return f"{value:.2f}"

    def _slider_changed(self, name: str) -> None:
        if name in self.value_labels:
            self.value_labels[name].configure(text=self._format_slider_value(name, self.vars[name].get()))

    def _build(self) -> None:
        outer = ttk.Frame(self.root)
        outer.pack(fill="both", expand=True)
        controls = ttk.Frame(outer)
        controls.pack(side="left", fill="y")
        plots = ttk.Frame(outer)
        plots.pack(side="right", fill="both", expand=True)

        sections = [
            ("History / beta", [
                ("History N", "history_n", 5, 250, 1),
                ("Seed", "seed", 1, 100, 1),
                ("Prior mu", "prior_mu", -1.5, 1.5, 0.05),
                ("Prior sigma", "prior_sigma", 0.1, 2.0, 0.05),
                ("True beta", "true_beta", -0.5, 2.0, 0.05),
                ("Edge scale", "edge_scale", 0.1, 2.5, 0.05),
                ("Wrong-side", "wrong_side_rate", 0.0, 0.5, 0.01),
            ]),
            ("Current markets", [
                ("Wealth", "wealth", 50, 1000, 10),
                ("A quote", "current_a_q", 0.05, 0.90, 0.01),
                ("A bot p*", "current_a_p", 0.05, 0.95, 0.01),
                ("A liquidity", "current_a_liq", 100, 2500, 25),
                ("B quote", "current_b_q", 0.05, 0.90, 0.01),
                ("B bot p*", "current_b_p", 0.05, 0.95, 0.01),
                ("B liquidity", "current_b_liq", 100, 2500, 25),
            ]),
            ("Future opportunities", [
                ("Poisson lam", "poisson_lambda", 0, 8, 0.1),
                ("Poisson max", "poisson_n_max", 3, 20, 1),
                ("Future quote", "future_q", 0.05, 0.90, 0.01),
                ("Future p*", "future_p", 0.05, 0.95, 0.01),
                ("Future wrong", "future_wrong_rate", 0.0, 0.5, 0.01),
                ("Future liq", "future_liq", 100, 2500, 25),
            ]),
        ]
        for title, fields in sections:
            box = ttk.LabelFrame(controls, text=title)
            box.pack(fill="x", padx=4, pady=2)
            for field in fields:
                self._add_slider(box, *field)
        ttk.Button(controls, text="Recompute now", command=self.refresh).pack(fill="x", padx=6, pady=(4, 2))
        ttk.Button(controls, text="Load article example", command=self.load_article_example).pack(fill="x", padx=6, pady=(0, 4))
        ttk.Label(controls, textvariable=self.status, wraplength=260).pack(fill="x", padx=6, pady=2)

        self.fig = Figure(figsize=(13.5, 8.8), dpi=100)
        self.plot_image: tk.PhotoImage | None = None
        self.image_label = ttk.Label(plots)
        self.image_label.pack(fill="both", expand=True)

    def _params_from_ui(self) -> Params:
        values = {name: var.get() for name, var in self.vars.items()}
        for key in ("history_n", "seed", "poisson_n_max"):
            if key in values:
                values[key] = int(round(values[key]))
        return replace(self.params, **values)

    def _set_slider(self, name: str, value: float) -> None:
        if name not in self.vars:
            return
        self.vars[name].set(value)
        self._slider_changed(name)

    def load_article_example(self) -> None:
        article = Params(
            wealth=300.0,
            current_a_q=0.35,
            current_a_p=0.70,
            current_a_liq=1000.0,
            current_b_q=0.25,
            current_b_p=0.60,
            current_b_liq=500.0,
            future_q=0.50,
            future_p=0.75,
            future_wrong_rate=0.0,
            future_liq=800.0,
            poisson_lambda=3.0,
            poisson_n_max=15,
        )
        for name in self.vars:
            self._set_slider(name, float(getattr(article, name)))
        self.use_article_posterior = True
        self.refresh()

    def refresh(self) -> None:
        article_mode = self.use_article_posterior
        self.use_article_posterior = False
        self.status.set("Computing...")
        self.root.update_idletasks()
        try:
            self.params = self._params_from_ui()
            result = run_simulation(self.params, article_posterior=article_mode)
            self._draw(result)
            prefix = "article posterior; " if article_mode else ""
            dynamic = stage_by_label(result, "Final dynamic")
            dynamic_commitment = float(np.sum(dynamic.costs))
            retained_cash = result.params.wealth - dynamic_commitment
            self.status.set(
                f"{prefix}beta mean={result.posterior.mean:.3f}, sd={result.posterior.std:.3f}; "
                f"corr(A,B)={result.corr_ab:.3f}; dynamic commitment="
                f"{dynamic_commitment:.1f}; retained cash={retained_cash:.1f}"
            )
        except Exception as exc:  # Keep GUI alive while sliders explore awkward corners.
            self.status.set(f"Computation failed: {exc}")

    def _draw(self, result: SimulationResult) -> None:
        draw_dashboard(self.fig, result)
        buffer = io.BytesIO()
        self.fig.savefig(buffer, format="png", dpi=100)
        data = base64.b64encode(buffer.getvalue())
        self.plot_image = tk.PhotoImage(data=data)
        self.image_label.configure(image=self.plot_image)



def article_check_params() -> Params:
    return Params(
        wealth=300,
        current_a_q=0.35,
        current_a_p=0.70,
        current_a_liq=1000,
        current_b_q=0.25,
        current_b_p=0.60,
        current_b_liq=500,
        future_q=0.50,
        future_p=0.75,
        future_liq=800,
        poisson_lambda=3.0,
        poisson_n_max=15,
        beta_grid_n=301,
    )


def run_check() -> None:
    params = article_check_params()
    result = run_simulation(params, article_posterior=True)
    static = result.stages[3]
    dynamic = result.stages[4]
    print("Article-posterior check")
    print(f"pbar: A={result.pbar[0]:.6f}, B={result.pbar[1]:.6f}, F={result.pbar[2]:.6f}")
    print(f"joint: {', '.join(f'{x:.6f}' for x in result.joint_current)}")
    print(f"corr(A,B): {result.corr_ab:.6f}")
    print(f"static shares: A={static.shares[0]:.3f}, B={static.shares[1]:.3f}; costs A={static.costs[0]:.3f}, B={static.costs[1]:.3f}")
    print(f"dynamic shares: A={dynamic.shares[0]:.3f}, B={dynamic.shares[1]:.3f}; costs A={dynamic.costs[0]:.3f}, B={dynamic.costs[1]:.3f}")
    if not (abs(static.shares[0] - 142.5) < 2.0 and abs(static.shares[1] - 109.3) < 2.0):
        raise SystemExit("static benchmark outside tolerance")
    if not (abs(dynamic.shares[0] - 127.6) < 6.0 and abs(dynamic.shares[1] - 99.9) < 6.0):
        raise SystemExit("dynamic benchmark outside tolerance")


def render_check(output_png: str) -> None:
    params = article_check_params()
    result = run_simulation(params, article_posterior=True)
    fig = Figure(figsize=(13.5, 8.8), dpi=100)
    draw_dashboard(fig, result)
    fig.savefig(output_png, dpi=120)
    print(output_png)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="run numerical checks without opening the GUI")
    parser.add_argument(
        "--render-check",
        metavar="OUTPUT_PNG",
        help="render the article-posterior dashboard to a PNG without opening the GUI",
    )
    args = parser.parse_args()
    if args.check:
        run_check()
        return
    if args.render_check:
        render_check(args.render_check)
        return
    root = tk.Tk()
    MarketApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
