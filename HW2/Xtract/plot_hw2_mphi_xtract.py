# -*- coding: utf-8 -*-
"""
從 hw2 data_knm.xlsx 工作表「hw2 data」讀取 6 組 moment–curvature（Xtract 匯出），
每組兩欄：Myy、Kyy（**1/m**）。Myy 在檔案中應為 **kN·m**（第 14 列單位列標示 kN-m）；
若舊檔仍為 **N·m**，程式會依單位列或數值量級自動換算為 kN·m 再繪圖，與 OpenSees（CSV 欄位 moment_kNm，或舊檔 moment_Nmm）一致。

另自 HW2/MainScript/results_test/HW2_case{n}_Data/MomentCurvature_*.csv 讀取 OpenSees
之 M–φ，與同案 Xtract 曲線疊畫於一圖（實線／虛線樣式同 data dealer.py 慣例）。

樣式參考 HW1/Xtract/data dealer.py（figsize、線寬、字級、格線、dpi）。

輸出：
  - results/MomentCurvature_Xtract_case1.png … case6.png（僅 Xtract）
  - results/MomentCurvature_Xtract_all_cases.png（Xtract 六案合圖，格式同 MainScript MomentCurvature_all_cases.png）
  - compa xtract/Mphi_*.png（六張雙曲線比較圖，對應 MainScript/compa 之檔名與樣式）
  - compa xtract/Mphi_fy420_fc28_fc35_fc42.png、Mphi_fy550_fc28_fc35_fc42.png（三曲線混凝土強度比較，同 HW2_test）
  - results/MomentCurvature_Op_vs_Xtract_case1.png … case6.png（OpenSees + Xtract）

選項（命令列）：
  --save-excel-knm  若 Excel 仍為 N·m，轉成 kN·m 後寫回 hw2 data_knm.xlsx（請先關閉 Excel 以免檔案鎖定）。

註：若 Excel 欄位順序與 HW2 分析 Case1～6 一致，可對應
    (28,420)、(28,550)、(35,420)、(35,550)、(42,420)、(42,550) MPa（僅用於圖題說明）。
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

_SCRIPT_DIR = Path(__file__).resolve().parent
_HW2_DIR = _SCRIPT_DIR.parent
EXCEL_PATH = _SCRIPT_DIR / "hw2 data_knm.xlsx"
SHEET_NAME = "hw2 data"
RESULTS_DIR = _SCRIPT_DIR / "results"
# 與 MainScript/compa 相同邏輯之雙曲線比較圖（Xtract 資料）
COMPA_XTRACT_DIR = _SCRIPT_DIR / "compa xtract"
# OpenSees 匯出（與 HW2/MainScript/HW2_test.py 之 CSV 欄位一致）
MAINSCRIPT_RESULTS = _HW2_DIR / "MainScript" / "results_test"

# 與 data dealer.py 對比圖之色與線型（OpenSees 實線、Xtract 虛線）
COLOR_OPENSEES = "#E94F37"
COLOR_XTRACT = "#2E86AB"

# 彎矩換算（所有繪圖之縱軸均為 kN·m，與 OpenSees CSV 一致）
# - Xtract Excel：Myy 欄應為 **kN·m**（舊檔 N·m 時由 normalize_excel_moment_to_knm 轉換）
# - OpenSees CSV：moment_kNm 為 **kN·m**；舊檔 moment_Nmm 為 **N·mm** → kN·m 須 ÷1e6
NMM_TO_KNM = 1e-6

# 第 17 列（0-based index 16）起為數值（首列可為 0,0）
DATA_START_ROW = 16

# 六組：每組 (Myy 欄, Kyy 欄)；Myy 單位 **kN·m**（檔案第 14 列），Kyy 單位 **1/m**
COL_PAIRS = [(0, 1), (2, 3), (4, 5), (6, 7), (8, 9), (10, 11)]

# 與 HW2 MainScript 前六案對應之 f'c、fy（MPa），供圖題使用
CASE_FC_FY_MP = (
    (28, 420),
    (28, 550),
    (35, 420),
    (35, 550),
    (42, 420),
    (42, 550),
)

# 與 HW2/MainScript/HW2_test.py 合圖圖例一致
MPHI_LEGEND_LABELS = (
    r"Case 1: $f_c'=28$, $f_y=420$ MPa",
    r"Case 2: $f_c'=28$, $f_y=550$ MPa",
    r"Case 3: $f_c'=35$, $f_y=420$ MPa",
    r"Case 4: $f_c'=35$, $f_y=550$ MPa",
    r"Case 5: $f_c'=42$, $f_y=420$ MPa",
    r"Case 6: $f_c'=42$, $f_y=550$ MPa",
)

# 各案 OpenSees CSV（相對於 MAINSCRIPT_RESULTS）
OPENSEES_CASE_CSV = (
    ("HW2_case1_Data", "MomentCurvature_Case1_28MPa_SD420.csv"),
    ("HW2_case2_Data", "MomentCurvature_Case2_28MPa_SD550.csv"),
    ("HW2_case3_Data", "MomentCurvature_Case3_35MPa_SD420.csv"),
    ("HW2_case4_Data", "MomentCurvature_Case4_35MPa_SD550.csv"),
    ("HW2_case5_Data", "MomentCurvature_Case5_42MPa_SD420.csv"),
    ("HW2_case6_Data", "MomentCurvature_Case6_42MPa_SD550.csv"),
)

# 與 HW2_test.export_mphi_comparison_charts 之 dict 鍵一致
CASE_KEYS = (
    "Case1_28MPa_SD420",
    "Case2_28MPa_SD550",
    "Case3_35MPa_SD420",
    "Case4_35MPa_SD550",
    "Case5_42MPa_SD420",
    "Case6_42MPa_SD550",
)


def normalize_excel_moment_to_knm(df: pd.DataFrame):
    """
    若工作表第 14 列標為 N·m，或數值量級顯示仍為 N·m，將 Myy 欄（0,2,4,6,8,10）
    自第 16 列起除以 1000，並將第 14 列改為 kN-m。
    已為 kN·m 的檔案不變更。
    回傳：(df, 是否曾轉換)
    """
    df = df.copy()
    moment_cols = [0, 2, 4, 6, 8, 10]
    unit_cell = df.iloc[14, 0]
    need_convert = False
    if isinstance(unit_cell, str):
        su = unit_cell.strip().replace("·", "-").replace(" ", "").upper()
        if su.startswith("KN") or "KN" in su:
            return df, False
        if su in ("N-M", "N.M"):
            need_convert = True
    if not need_convert:
        s0 = pd.to_numeric(df.iloc[16:, 0], errors="coerce").dropna()
        if len(s0) and float(np.nanmax(s0)) > 5000.0:
            need_convert = True
    if not need_convert:
        return df, False
    for c in moment_cols:
        for r in range(16, len(df)):
            val = df.iloc[r, c]
            if pd.isna(val):
                continue
            try:
                x = float(val)
            except (TypeError, ValueError):
                continue
            df.iloc[r, c] = x * 1e-3
        df.iloc[14, c] = "kN-m"
    return df, True


def load_mphi_series_xtract(df: pd.DataFrame, col_myy: int, col_kyy: int):
    """
    自 Xtract Excel 讀取一組 M–φ（Myy 須已為 **kN·m**）。
    Kyy：**1/m**。
    回傳：曲率 κ (1/m)、彎矩 **(kN·m)**。
    """
    blk = df.iloc[DATA_START_ROW:, [col_myy, col_kyy]].copy()
    blk = blk.apply(pd.to_numeric, errors="coerce").dropna()
    moment_kn_m = blk.iloc[:, 0].values.astype(float)
    kappa = blk.iloc[:, 1].values.astype(float)
    return kappa, moment_kn_m


def load_opensees_mphi_csv(csv_path: Path):
    """
    OpenSees 匯出 CSV：curvature_1_per_m、moment_kNm（**kN·m**）。
    若為舊版 moment_Nmm（**N·mm**），仍會自動換算為 kN·m。
    """
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
        raise KeyError(
            f"CSV 需含 moment_kNm 或 moment_Nmm 欄位：{csv_path}"
        )
    return kappa, moment_kn_m


def build_xtract_curves_by_case(df_excel: pd.DataFrame):
    """回傳與 HW2_test.export_mphi_comparison_charts 相同鍵名之 (N,2) 曲線（κ 1/m、M kN·m）。"""
    out = {}
    for i, key in enumerate(CASE_KEYS):
        col_m, col_k = COL_PAIRS[i]
        kappa, m_knm = load_mphi_series_xtract(df_excel, col_m, col_k)
        out[key] = np.column_stack([kappa, m_knm])
    return out


def _curve_moment_knm(curve):
    """curve 第二欄為彎矩（kN·m）。"""
    return np.asarray(curve[:, 1], dtype=float)


def plot_moment_curvature_compare_two(curve_a, label_a, curve_b, label_b, out_png: Path, title_line: str):
    """
    同一圖內兩條 M–φ；格式與 HW2_test.plot_moment_curvature_compare_two 一致。
    curve_*: (N,2)，曲率 1/m、彎矩 kN·m
    """
    colors = ["#1f77b4", "#ff7f0e"]
    fig, ax = plt.subplots(figsize=(7, 4.8))
    for curve, lab, c in ((curve_a, label_a, colors[0]), (curve_b, label_b, colors[1])):
        if curve is None or len(curve) == 0:
            continue
        k = curve[:, 0]
        m_knm = _curve_moment_knm(curve)
        ax.plot(k, m_knm, "-", linewidth=1.6, color=c, label=lab)
    ax.set_xlabel(r"Curvature $\phi$ (m$^{-1}$)")
    ax.set_ylabel("Moment (kN·m)")
    ax.set_title(title_line)
    ax.grid(True, linestyle=":", alpha=0.7)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=8, prune="both"))
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.3f"))
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=240, bbox_inches="tight")
    plt.close(fig)
    print(f"已匯出比較圖（Xtract）：{out_png}")


def plot_moment_curvature_compare_three(
    curve_a, label_a, curve_b, label_b, curve_c, label_c, out_png: Path, title_line: str
):
    """
    同一圖內三條 M–φ；與 HW2_test.plot_moment_curvature_compare_three 一致。
    """
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    fig, ax = plt.subplots(figsize=(7, 4.8))
    for curve, lab, c in (
        (curve_a, label_a, colors[0]),
        (curve_b, label_b, colors[1]),
        (curve_c, label_c, colors[2]),
    ):
        if curve is None or len(curve) == 0:
            continue
        k = curve[:, 0]
        m_knm = _curve_moment_knm(curve)
        ax.plot(k, m_knm, "-", linewidth=1.6, color=c, label=lab)
    ax.set_xlabel(r"Curvature $\phi$ (m$^{-1}$)")
    ax.set_ylabel("Moment (kN·m)")
    ax.set_title(title_line)
    ax.grid(True, linestyle=":", alpha=0.7)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=8, prune="both"))
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.3f"))
    ax.legend(loc="lower right", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=240, bbox_inches="tight")
    plt.close(fig)
    print(f"已匯出比較圖（Xtract，三曲線）：{out_png}")


def export_xtract_mphi_comparison_charts(compa_dir: Path, df_excel: pd.DataFrame):
    """
    產生六張雙曲線比較圖至 compa_dir；配對與檔名同 HW2_test.export_mphi_comparison_charts。
    """
    compa_dir.mkdir(parents=True, exist_ok=True)
    curves_by_case = build_xtract_curves_by_case(df_excel)
    g = curves_by_case.get

    def pair(key_a, key_b, lab_a, lab_b, fname, title_line):
        plot_moment_curvature_compare_two(
            g(key_a),
            lab_a,
            g(key_b),
            lab_b,
            compa_dir / fname,
            title_line,
        )

    t_steel = "Moment–curvature (Different steel grade comparison)"
    t_conc = "Moment–curvature (Different concrete strength)"

    pair(
        "Case1_28MPa_SD420",
        "Case2_28MPa_SD550",
        r"$f_c'=28$ MPa, $f_y=420$ MPa",
        r"$f_c'=28$ MPa, $f_y=550$ MPa",
        "Mphi_fc28_fy420_vs_fy550.png",
        t_steel,
    )
    pair(
        "Case3_35MPa_SD420",
        "Case4_35MPa_SD550",
        r"$f_c'=35$ MPa, $f_y=420$ MPa",
        r"$f_c'=35$ MPa, $f_y=550$ MPa",
        "Mphi_fc35_fy420_vs_fy550.png",
        t_steel,
    )
    pair(
        "Case5_42MPa_SD420",
        "Case6_42MPa_SD550",
        r"$f_c'=42$ MPa, $f_y=420$ MPa",
        r"$f_c'=42$ MPa, $f_y=550$ MPa",
        "Mphi_fc42_fy420_vs_fy550.png",
        t_steel,
    )
    pair(
        "Case1_28MPa_SD420",
        "Case3_35MPa_SD420",
        r"$f_c'=28$ MPa, $f_y=420$ MPa",
        r"$f_c'=35$ MPa, $f_y=420$ MPa",
        "Mphi_fy420_fc28_vs_fc35.png",
        t_conc,
    )
    pair(
        "Case1_28MPa_SD420",
        "Case5_42MPa_SD420",
        r"$f_c'=28$ MPa, $f_y=420$ MPa",
        r"$f_c'=42$ MPa, $f_y=420$ MPa",
        "Mphi_fy420_fc28_vs_fc42.png",
        t_conc,
    )
    pair(
        "Case3_35MPa_SD420",
        "Case5_42MPa_SD420",
        r"$f_c'=35$ MPa, $f_y=420$ MPa",
        r"$f_c'=42$ MPa, $f_y=420$ MPa",
        "Mphi_fy420_fc35_vs_fc42.png",
        t_conc,
    )


def export_xtract_mphi_concrete_strength_three_way(compa_dir: Path, df_excel: pd.DataFrame):
    """
    固定 fy、比較 f'c = 28/35/42 MPa 之三曲線圖；檔名與 HW2_test.export_mphi_concrete_strength_three_way 一致。
    """
    compa_dir.mkdir(parents=True, exist_ok=True)
    g = build_xtract_curves_by_case(df_excel).get
    t_conc = "Moment–curvature (Different concrete strength)"
    plot_moment_curvature_compare_three(
        g("Case1_28MPa_SD420"),
        r"$f_c'=28$ MPa, $f_y=420$ MPa",
        g("Case3_35MPa_SD420"),
        r"$f_c'=35$ MPa, $f_y=420$ MPa",
        g("Case5_42MPa_SD420"),
        r"$f_c'=42$ MPa, $f_y=420$ MPa",
        compa_dir / "Mphi_fy420_fc28_fc35_fc42.png",
        t_conc,
    )
    plot_moment_curvature_compare_three(
        g("Case2_28MPa_SD550"),
        r"$f_c'=28$ MPa, $f_y=550$ MPa",
        g("Case4_35MPa_SD550"),
        r"$f_c'=35$ MPa, $f_y=550$ MPa",
        g("Case6_42MPa_SD550"),
        r"$f_c'=42$ MPa, $f_y=550$ MPa",
        compa_dir / "Mphi_fy550_fc28_fc35_fc42.png",
        t_conc,
    )


def _style_mphi_figure():
    """僅 Xtract 單線圖用：不顯示 legend（標題已含案別與材料）。"""
    plt.xlabel(r"Curvature $\phi$ (m$^{-1}$)", fontsize=12)
    plt.ylabel("Moment (kN·m)", fontsize=12)
    plt.grid(True)


def plot_xtract_all_cases(df_excel: pd.DataFrame, out_png: Path):
    """
    六案 Xtract 曲線同圖；格式對齊 HW2_test.plot_moment_curvature_combined
    （figsize、色票、格線、座標刻度、legend、dpi）。
    """
    colors = [
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#d62728",
        "#9467bd",
        "#8c564b",
    ]
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(11, 5.8))
    for i, (((col_m, col_k), _), label) in enumerate(
        zip(zip(COL_PAIRS, CASE_FC_FY_MP), MPHI_LEGEND_LABELS)
    ):
        kappa, m_knm = load_mphi_series_xtract(df_excel, col_m, col_k)
        if len(kappa) == 0:
            continue
        c = colors[i % len(colors)]
        ax.plot(kappa, m_knm, "-", linewidth=1.5, color=c, label=label)
    ax.set_xlabel(r"Curvature $\phi$ (m$^{-1}$)")
    ax.set_ylabel("Moment (kN·m)")
    ax.set_title("Moment–curvature (Xtract) all cases")
    ax.grid(True, linestyle=":", alpha=0.7)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=8, prune="both"))
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.3f"))
    ax.legend(loc="best", fontsize=7.5)
    fig.tight_layout()
    fig.savefig(out_png, dpi=240, bbox_inches="tight")
    plt.close(fig)
    print(f"已儲存（Xtract 合圖）：{out_png}")


def plot_op_vs_xtract_all_cases(df_excel: pd.DataFrame):
    """六案：OpenSees 與 Xtract 兩條曲線同圖。"""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    for idx, (
        ((col_m, col_k), (fc, fy)),
        (sub_dir, csv_name),
    ) in enumerate(zip(zip(COL_PAIRS, CASE_FC_FY_MP), OPENSEES_CASE_CSV), start=1):
        kappa_x, m_x = load_mphi_series_xtract(df_excel, col_m, col_k)
        csv_path = MAINSCRIPT_RESULTS / sub_dir / csv_name
        kappa_o, m_o = load_opensees_mphi_csv(csv_path)

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
            rf"Moment–curvature Case {idx}: $f_c'={fc}$ MPa, $f_y={fy}$ MPa (OpenSees vs Xtract)",
            fontsize=14,
        )
        plt.legend(loc="best", fontsize=11)
        plt.grid(True)
        plt.xlim(0, x_max)
        plt.ylim(0, y_max)
        out = RESULTS_DIR / f"MomentCurvature_Op_vs_Xtract_case{idx}.png"
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"已儲存（疊圖）：{out}")


def main(save_excel_knm: bool = False):
    if not EXCEL_PATH.is_file():
        raise FileNotFoundError(f"找不到 Excel：{EXCEL_PATH}")

    df = pd.read_excel(EXCEL_PATH, sheet_name=SHEET_NAME, header=None)
    df, converted = normalize_excel_moment_to_knm(df)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    for idx, ((col_m, col_k), (fc, fy)) in enumerate(zip(COL_PAIRS, CASE_FC_FY_MP), start=1):
        kappa, m_knm = load_mphi_series_xtract(df, col_m, col_k)
        plt.figure(figsize=(10, 6))
        plt.plot(kappa, m_knm, linewidth=2, color="#2E86AB")
        plt.title(
            rf"Moment–curvature (Xtract) Case {idx}: $f_c'={fc}$ MPa, $f_y={fy}$ MPa",
            fontsize=14,
        )
        _style_mphi_figure()
        out = RESULTS_DIR / f"MomentCurvature_Xtract_case{idx}.png"
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"已儲存：{out} (n={len(kappa)} points)")

    plot_xtract_all_cases(df, RESULTS_DIR / "MomentCurvature_Xtract_all_cases.png")
    export_xtract_mphi_comparison_charts(COMPA_XTRACT_DIR, df)
    export_xtract_mphi_concrete_strength_three_way(COMPA_XTRACT_DIR, df)
    plot_op_vs_xtract_all_cases(df)

    if save_excel_knm:
        if converted:
            try:
                with pd.ExcelWriter(EXCEL_PATH, engine="openpyxl") as w:
                    df.to_excel(w, sheet_name=SHEET_NAME, header=False, index=False)
                print(f"已將 Myy 改為 kN·m 並寫回：{EXCEL_PATH}")
            except PermissionError:
                print(
                    f"無法寫入（檔案可能被 Excel 或其他程式開啟）：{EXCEL_PATH}\n"
                    "請關閉後再執行 --save-excel-knm。"
                )
        else:
            print("Excel 之 Myy 已是 kN·m，未變更檔案。")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="繪製 Xtract / OpenSees M–φ 曲線（彎矩單位 kN·m）")
    ap.add_argument(
        "--save-excel-knm",
        action="store_true",
        help="將 hw2 data_knm.xlsx 內 N·m 之 Myy 轉為 kN·m 並存檔（請先關閉 Excel）",
    )
    args = ap.parse_args()
    main(save_excel_knm=args.save_excel_knm)
