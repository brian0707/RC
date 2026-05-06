# -*- coding: utf-8 -*-
"""
將 Xtract 匯出之 tab 檔（case1 test.txt）與 OpenSees CSV 疊繪於一圖。
輸出：與本資料夾 results/MomentCurvature_Op_vs_Xtract_case1.png 相同樣式，檔名 test.png。
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_SCRIPT_DIR = Path(__file__).resolve().parent
_HW2_DIR = _SCRIPT_DIR.parent

# 與 plot_hw2_mphi_xtract.py 一致
COLOR_OPENSEES = "#E94F37"
COLOR_XTRACT = "#2E86AB"
NMM_TO_KNM = 1e-6

# 輸入／輸出
XTRACT_TXT = _SCRIPT_DIR / "case1 test.txt"
OPENSEES_CSV = (
    _HW2_DIR
    / "MainScript"
    / "results_test"
    / "HW2_case1_Data"
    / "MomentCurvature_Case1_28MPa_SD420.csv"
)
OUT_PNG = _SCRIPT_DIR / "test.png"

# Case 1 標題用（與 MomentCurvature_Op_vs_Xtract_case1.png）
FC_MPA = 28
FY_MPA = 420


def load_xtract_case1_txt(path: Path):
    """
    case1 test.txt：表頭後為兩欄 tab，Myy 為 N·m、Kxx 為 1/m（見第 15 列單位列）。
    回傳：曲率 κ (1/m)、彎矩 (kN·m)。
    """
    if not path.is_file():
        raise FileNotFoundError(f"找不到 Xtract 檔：{path}")
    df = pd.read_csv(path, sep="\t", skiprows=16, header=None, engine="python")
    df = df.iloc[:, :2].apply(pd.to_numeric, errors="coerce").dropna()
    moment_nm = df.iloc[:, 0].values.astype(float)
    kappa = df.iloc[:, 1].values.astype(float)
    moment_knm = moment_nm * 1e-3
    kappa = np.maximum(kappa, 0.0)
    return kappa, moment_knm


def load_opensees_mphi_csv(csv_path: Path):
    """OpenSees 匯出 CSV：curvature_1_per_m、moment_kNm（kN·m）。"""
    if not csv_path.is_file():
        raise FileNotFoundError(f"找不到 OpenSees CSV：{csv_path}")
    df = pd.read_csv(csv_path)
    df = df.apply(pd.to_numeric, errors="coerce").dropna()
    kappa = df["curvature_1_per_m"].values.astype(float)
    if "moment_kNm" in df.columns:
        moment_kn_m = df["moment_kNm"].values.astype(float)
    elif "moment_Nmm" in df.columns:
        moment_kn_m = df["moment_Nmm"].values.astype(float) * NMM_TO_KNM
    else:
        raise KeyError(f"CSV 需含 moment_kNm 或 moment_Nmm：{csv_path}")
    return kappa, moment_kn_m


def main():
    kappa_x, m_x = load_xtract_case1_txt(XTRACT_TXT)
    kappa_o, m_o = load_opensees_mphi_csv(OPENSEES_CSV)

    y_max = max(float(np.max(m_x)), float(np.max(m_o))) * 1.02
    x_max = max(float(np.max(np.abs(kappa_x))), float(np.max(np.abs(kappa_o)))) * 1.02

    plt.figure(figsize=(10, 6))
    plt.plot(
        kappa_o,
        m_o,
        color=COLOR_OPENSEES,
        linewidth=2,
        linestyle="-",
        label="OpenSees",
    )
    plt.plot(
        kappa_x,
        m_x,
        color=COLOR_XTRACT,
        linewidth=2,
        linestyle="--",
        label="Xtract",
    )
    plt.xlabel(r"Curvature $\phi$ (m$^{-1}$)", fontsize=12)
    plt.ylabel("Moment (kN·m)", fontsize=12)
    plt.title(
        rf"Moment–curvature Case 1: $f_c'={FC_MPA}$ MPa, $f_y={FY_MPA}$ MPa (OpenSees vs Xtract)",
        fontsize=14,
    )
    plt.legend(loc="best", fontsize=11)
    plt.grid(True)
    plt.xlim(0, x_max)
    plt.ylim(0, y_max)
    plt.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"已儲存：{OUT_PNG}")


if __name__ == "__main__":
    main()
