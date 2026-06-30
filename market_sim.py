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
    mc_paths: int = 4000
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
    wealth_samples: dict[str, np.ndarray]


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


def simulate_wealth(
    params: Params,
    markets: tuple[Market, Market, Market],
    posterior: Posterior,
    strategy: StrategyResult,
    rng: np.random.Generator,
) -> np.ndarray:
    beta = rng.choice(posterior.beta, size=params.mc_paths, p=posterior.weights)
    probs = []
    for m in markets[:2]:
        probs.append(expit(logit(m.q) + beta * (logit(m.p_star) - logit(m.q))))
    ya = rng.binomial(1, probs[0])
    yb = rng.binomial(1, probs[1])
    return params.wealth - np.sum(strategy.costs) + strategy.shares[0] * ya + strategy.shares[1] * yb


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
    rng = np.random.default_rng(params.seed + 123)
    wealth_samples = {
        "Raw frictionless Kelly": simulate_wealth(params, markets, posterior, stages[0], rng),
        "Joint posterior": simulate_wealth(params, markets, posterior, stages[3], rng),
        "Final dynamic": simulate_wealth(params, markets, posterior, stages[4], rng),
    }
    return SimulationResult(params, hist, posterior, markets, pbar, joint, corr, stages, future_policy, wealth_samples)


class MarketApp:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("Prediction Market Kelly / Calibration / Poisson Simulator")
        self.params = Params()
        self.vars: dict[str, tk.DoubleVar] = {}
        self.value_labels: dict[str, ttk.Label] = {}
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
        ttk.Button(controls, text="Recompute now", command=self.refresh).pack(fill="x", padx=6, pady=4)
        ttk.Label(controls, textvariable=self.status, wraplength=260).pack(fill="x", padx=6, pady=2)

        self.fig = Figure(figsize=(12, 8), dpi=100)
        self.plot_image: tk.PhotoImage | None = None
        self.image_label = ttk.Label(plots)
        self.image_label.pack(fill="both", expand=True)

    def _params_from_ui(self) -> Params:
        values = {name: var.get() for name, var in self.vars.items()}
        for key in ("history_n", "seed", "poisson_n_max"):
            if key in values:
                values[key] = int(round(values[key]))
        return replace(self.params, **values)

    def refresh(self) -> None:
        self.status.set("Computing...")
        self.root.update_idletasks()
        try:
            self.params = self._params_from_ui()
            result = run_simulation(self.params)
            self._draw(result)
            self.status.set(
                f"beta mean={result.posterior.mean:.3f}, sd={result.posterior.std:.3f}; "
                f"corr(A,B)={result.corr_ab:.3f}; final costs="
                f"{result.stages[-1].costs[0]:.1f}, {result.stages[-1].costs[1]:.1f}"
            )
        except Exception as exc:  # Keep GUI alive while sliders explore awkward corners.
            self.status.set(f"Computation failed: {exc}")

    def _draw(self, result: SimulationResult) -> None:
        self.fig.clear()
        gs = self.fig.add_gridspec(3, 3)
        ax_hist = self.fig.add_subplot(gs[0, 0])
        ax_beta = self.fig.add_subplot(gs[0, 1])
        ax_stage = self.fig.add_subplot(gs[0, 2])
        ax_policy = self.fig.add_subplot(gs[1, 0])
        ax_wealth = self.fig.add_subplot(gs[1, 1])
        ax_table = self.fig.add_subplot(gs[1:, 2])
        ax_joint = self.fig.add_subplot(gs[2, 0])
        ax_cost = self.fig.add_subplot(gs[2, 1])

        h = result.history
        ax_hist.scatter(h.q[h.y == 0], h.p_star[h.y == 0], s=18, alpha=0.7, label="resolved NO")
        ax_hist.scatter(h.q[h.y == 1], h.p_star[h.y == 1], s=18, alpha=0.7, label="resolved YES")
        ax_hist.plot([0, 1], [0, 1], color="black", lw=1)
        ax_hist.set_title("Resolved history")
        ax_hist.set_xlabel("market quote q")
        ax_hist.set_ylabel("bot forecast p*")
        ax_hist.legend(fontsize=8)

        beta = result.posterior.beta
        density = result.posterior.weights / max(EPS, beta[1] - beta[0])
        ax_beta.plot(beta, density, color="#244c9a")
        ax_beta.axvline(0, color="gray", lw=1)
        ax_beta.axvline(1, color="gray", lw=1, ls="--")
        ax_beta.set_title("Posterior over beta")
        ax_beta.set_xlabel("beta")

        labels = [s.label.replace(" ", "\n") for s in result.stages]
        x = np.arange(len(result.stages))
        width = 0.38
        ax_stage.bar(x - width / 2, [s.shares[0] for s in result.stages], width, label="A")
        ax_stage.bar(x + width / 2, [s.shares[1] for s in result.stages], width, label="B")
        ax_stage.set_xticks(x)
        ax_stage.set_xticklabels(labels, fontsize=7)
        ax_stage.set_title("Recommended shares by stage")
        ax_stage.legend(fontsize=8)

        if result.future_policy:
            ns = [p[0] for p in result.future_policy]
            us = [p[1] for p in result.future_policy]
            ax_policy.plot(ns, us, marker="o")
        ax_policy.set_title("Second-stage F shares per market")
        ax_policy.set_xlabel("future arrivals n")
        ax_policy.set_ylabel("u_n")

        for name, samples in result.wealth_samples.items():
            ax_wealth.hist(samples, bins=35, alpha=0.45, density=True, label=name)
        ax_wealth.set_title("Current A/B terminal wealth")
        ax_wealth.legend(fontsize=7)

        joint_matrix = np.array([[result.joint_current[0], result.joint_current[1]], [result.joint_current[2], result.joint_current[3]]])
        im = ax_joint.imshow(joint_matrix, cmap="Blues")
        for i in range(2):
            for j in range(2):
                ax_joint.text(j, i, f"{joint_matrix[i, j]:.3f}", ha="center", va="center")
        ax_joint.set_xticks([0, 1], labels=["B=0", "B=1"])
        ax_joint.set_yticks([0, 1], labels=["A=0", "A=1"])
        ax_joint.set_title(f"Joint predictive; corr={result.corr_ab:.3f}")
        self.fig.colorbar(im, ax=ax_joint, fraction=0.046, pad=0.04)

        ax_cost.bar(x - width / 2, [s.costs[0] for s in result.stages], width, label="A")
        ax_cost.bar(x + width / 2, [s.costs[1] for s in result.stages], width, label="B")
        ax_cost.set_xticks(x)
        ax_cost.set_xticklabels(labels, fontsize=7)
        ax_cost.set_title("Cash committed by stage")
        ax_cost.legend(fontsize=8)

        ax_table.axis("off")
        rows = [
            ["beta mean", f"{result.posterior.mean:.4f}"],
            ["beta std", f"{result.posterior.std:.4f}"],
            ["Pr(beta<0)", f"{result.posterior.prob_lt(0):.4f}"],
            ["Pr(0<=beta<=1)", f"{result.posterior.prob_between(0, 1):.4f}"],
            ["Pr(beta>1)", f"{1.0 - result.posterior.prob_lt(1):.4f}"],
            ["pbar A", f"{result.pbar[0]:.4f}"],
            ["pbar B", f"{result.pbar[1]:.4f}"],
            ["pbar F", f"{result.pbar[2]:.4f}"],
        ]
        rows += [[s.label, f"A {s.shares[0]:.1f} / B {s.shares[1]:.1f}"] for s in result.stages]
        table = ax_table.table(cellText=rows, colLabels=["quantity", "value"], loc="center", cellLoc="left")
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1.0, 1.25)

        self.fig.tight_layout()
        buffer = io.BytesIO()
        self.fig.savefig(buffer, format="png", dpi=100)
        data = base64.b64encode(buffer.getvalue())
        self.plot_image = tk.PhotoImage(data=data)
        self.image_label.configure(image=self.plot_image)


def run_check() -> None:
    params = Params(
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
        mc_paths=500,
    )
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="run numerical checks without opening the GUI")
    args = parser.parse_args()
    if args.check:
        run_check()
        return
    root = tk.Tk()
    MarketApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
