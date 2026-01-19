#!/usr/bin/env python3
"""
壁引き抜き法による安息角シミュレーション結果のプロット
摩擦係数と転がり摩擦係数が安息角に与える影響を可視化
"""

import matplotlib.pyplot as plt
import numpy as np

# 日本語フォント設定（macOS）
plt.rcParams['font.family'] = 'Hiragino Sans'
plt.rcParams['axes.unicode_minus'] = False

# === 摩擦係数 vs 安息角 ===
friction_coeffs = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
angles_friction = np.array([0.93, 9.64, 11.47, 12.08, 11.80, 12.57])  # 度

# プロット1: 摩擦係数 vs 安息角（シミュレーション結果のみ）
fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(friction_coeffs, angles_friction, s=100, c='tab:blue', 
           marker='o', zorder=5, label='Simulation Results')
ax.plot(friction_coeffs, angles_friction, 'tab:blue', alpha=0.7, linewidth=1.5)

ax.set_xlabel(r'Sliding Friction Coefficient $\mu$ [-]', fontsize=12)
ax.set_ylabel(r'Angle of Repose $\theta_r$ [deg]', fontsize=12)
ax.set_xlim(-0.05, 1.05)
ax.set_ylim(0, 15)
ax.grid(True, alpha=0.3)
ax.legend(loc='upper left', fontsize=11)
ax.set_title('Sliding Friction Coefficient vs Angle of Repose', fontsize=14)

plt.tight_layout()
plt.savefig('friction_vs_repose_angle.png', dpi=150, bbox_inches='tight')
plt.close()

# === 転がり摩擦係数 vs 安息角 ===
rolling_friction_coeffs = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
angles_rolling = np.array([1.51, 8.24, 8.90, 9.43, 9.51, 9.69])  # 度

# プロット2: 転がり摩擦係数 vs 安息角
fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(rolling_friction_coeffs, angles_rolling, s=100, c='tab:orange', 
           marker='s', zorder=5, label='Simulation Results')
ax.plot(rolling_friction_coeffs, angles_rolling, 'tab:orange', alpha=0.7)

ax.set_xlabel(r'Rolling Friction Coefficient $\mu_r$ [-]', fontsize=12)
ax.set_ylabel(r'Angle of Repose $\theta_r$ [deg]', fontsize=12)
ax.set_xlim(-0.05, 1.05)
ax.set_ylim(0, 15)
ax.grid(True, alpha=0.3)
ax.legend(loc='upper left', fontsize=11)
ax.set_title('Rolling Friction Coefficient vs Angle of Repose', fontsize=14)

plt.tight_layout()
plt.savefig('rolling_friction_vs_repose_angle.png', dpi=150, bbox_inches='tight')
plt.close()

# === 両方を並べたプロット ===
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# 左: 滑り摩擦
ax1 = axes[0]
ax1.scatter(friction_coeffs, angles_friction, s=100, c='tab:blue', 
            marker='o', zorder=5, label='Simulation Results')
ax1.plot(friction_coeffs, angles_friction, 'tab:blue', alpha=0.7, linewidth=1.5)
ax1.set_xlabel(r'Sliding Friction Coefficient $\mu$ [-]', fontsize=12)
ax1.set_ylabel(r'Angle of Repose $\theta_r$ [deg]', fontsize=12)
ax1.set_xlim(-0.05, 1.05)
ax1.set_ylim(0, 15)
ax1.grid(True, alpha=0.3)
ax1.legend(loc='upper left', fontsize=10)
ax1.set_title('(a) Effect of Sliding Friction', fontsize=12)

# 右: 転がり摩擦
ax2 = axes[1]
ax2.scatter(rolling_friction_coeffs, angles_rolling, s=100, c='tab:orange', 
            marker='s', zorder=5, label='Simulation Results')
ax2.plot(rolling_friction_coeffs, angles_rolling, 'tab:orange', alpha=0.7, linewidth=1.5)
ax2.set_xlabel(r'Rolling Friction Coefficient $\mu_r$ [-]', fontsize=12)
ax2.set_ylabel(r'Angle of Repose $\theta_r$ [deg]', fontsize=12)
ax2.set_xlim(-0.05, 1.05)
ax2.set_ylim(0, 15)
ax2.grid(True, alpha=0.3)
ax2.legend(loc='upper left', fontsize=10)
ax2.set_title('(b) Effect of Rolling Friction', fontsize=12)

plt.tight_layout()
plt.savefig('repose_angle_combined.png', dpi=150, bbox_inches='tight')
plt.close()

print("Plots saved successfully:")
print("  - friction_vs_repose_angle.png")
print("  - rolling_friction_vs_repose_angle.png")
print("  - repose_angle_combined.png")
