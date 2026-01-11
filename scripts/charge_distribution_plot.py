import numpy as np
import matplotlib.pyplot as plt

paths = {
    "bimodal": "results/charge_distribution_study/phase1_filling/bimodal/filled_particles.dat",
    "normal":  "results/charge_distribution_study/phase1_filling/normal/filled_particles.dat",
    "uniform": "results/charge_distribution_study/phase1_filling/uniform/filled_particles.dat",
}

# 読み込み
qs = {k: np.loadtxt(p, comments="#")[:, 3] for k, p in paths.items()}

# 3分布で同じ範囲にする（比較しやすい）
q_all = np.concatenate(list(qs.values()))
qmin, qmax = np.min(q_all), np.max(q_all)
bins = np.linspace(qmin, qmax, 100)

fig, axes = plt.subplots(1, 3, figsize=(12, 3), sharey=True)
for ax, (name, q) in zip(axes, qs.items()):
    ax.hist(q, bins=bins, alpha=0.9)
    ax.set_title(name)
    ax.set_xlabel("charge [C]")
    ax.grid(True, alpha=0.2)

axes[0].set_ylabel("count")
fig.suptitle("Charge distributions at filled state (Phase 1)")
fig.tight_layout()
fig.savefig("charge_distribution_comparison.png", dpi=200)
print("Saved: charge_distribution_comparison.png")