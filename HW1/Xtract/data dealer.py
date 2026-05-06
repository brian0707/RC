# -*- coding: utf-8 -*-
"""
從 case1.xlsx 工作表「ALL」讀取四組 case 之 P–strain（由左至右：case1～case4），
工作表「case0」讀取 baseline。

綜合圖：
  - Xtract_case0_to_case4.png：case0～case4 五條曲線
  - Xtract_case1_to_case4.png：case1～case4 四條曲線
  - Xtract data.png：與 case1～case4 綜合圖相同（沿用舊檔名）

工作表「case0」另輸出單線圖 case0.png（樣式同上）。

另將各 case 與 OpenSees（宥宏/RRRRRR/PDeltaData_*.txt）疊圖對比，輸出至 HW1/Comparison/。

Case1～4 OpenSees 與 Xtract 八條曲線同圖：Case1~4 Op_vs_Xt.png（同 Case 同色；OpenSees 實線、Xtract 虛線）。
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_SCRIPT_DIR = Path(__file__).resolve().parent
_HW1_DIR = _SCRIPT_DIR.parent
EXCEL_PATH = _SCRIPT_DIR / "case1.xlsx"
OUT_PATH = _SCRIPT_DIR / "Xtract data.png"
OUT_PATH_CASE0 = _SCRIPT_DIR / "case0.png"
OUT_PATH_CASE0_TO_CASE4 = _SCRIPT_DIR / "Xtract_case0_to_case4.png"
OUT_PATH_CASE1_TO_CASE4 = _SCRIPT_DIR / "Xtract_case1_to_case4.png"
OUT_PATH_OP_VS_XT = _SCRIPT_DIR / "Case1~4 Op_vs_Xt.png"
# OpenSees 曲線資料（與 HW1_YH 輸出檔名一致）
RRRRRR_DIR = _HW1_DIR / "宥宏" / "RRRRRR"
COMPARISON_DIR = _HW1_DIR / "Comparison"

# 與 HW1_YH.py 圖表一致之標籤
LABELS = [
    "Case1_28MPa_SD420",
    "Case2_28MPa_SD550",
    "Case3_70MPa_SD420",
    "Case4_70MPa_SD550",
]
# ALL：第 6 欄起每兩欄一組 P、strain（0-based），中間空欄略過
COL_PAIRS = [(6, 7), (9, 10), (12, 13), (15, 16)]
# 試算表第 11 列（index 10）起為數值列（首列可為 0,0）
DATA_START_ROW = 10
# 工作表 case0：第 7 列（index 6）起為 P、strain 數值
CASE0_DATA_START_ROW = 6

# 對比圖：Xtract / OpenSees 線色
COLOR_XTRACT = "#2E86AB"
COLOR_OPENSEES = "#E94F37"
# Case1～4 同圖：各 Case 一色（與 matplotlib 預設四色一致）
CASE_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]


def load_xtract_series(df, col_p, col_strain):
    blk = df.iloc[DATA_START_ROW:, [col_p, col_strain]].copy()
    blk = blk.apply(pd.to_numeric, errors="coerce").dropna()
    load = blk.iloc[:, 0].values
    strain = blk.iloc[:, 1].values
    return strain, load


def load_case0_series(df_case0):
    """工作表 case0：欄 0 = P (kN)，欄 1 = Axial Strain。"""
    blk = df_case0.iloc[CASE0_DATA_START_ROW:, [0, 1]].copy()
    blk = blk.apply(pd.to_numeric, errors="coerce").dropna()
    load = blk.iloc[:, 0].values
    strain = blk.iloc[:, 1].values
    return strain, load


def load_opensees_txt(path: Path):
    """PDeltaData_*.txt：第一列為表頭，欄位為 Axial_Strain、Axial_Load(kN)"""
    if not path.is_file():
        raise FileNotFoundError(path)
    data = np.loadtxt(path, skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data[:, 0], data[:, 1]


def _max_load_case(df, case_index):
    """單一 case：Xtract 與 OpenSees 兩條線之軸力最大值（kN）。"""
    label = LABELS[case_index]
    col_p, col_strain = COL_PAIRS[case_index]
    _, load_x = load_xtract_series(df, col_p, col_strain)
    _, load_o = load_opensees_txt(RRRRRR_DIR / f"PDeltaData_{label}.txt")
    return max(float(np.max(load_x)), float(np.max(load_o)))


def comparison_ylim_top_from_case34(df):
    """
    四張對比圖共用 Y 軸上限：取 Case3、Case4 各自（Xtract+OpenSees）軸力上限再取 max，
    並略留邊距。
    """
    m3 = _max_load_case(df, 2)  # Case3
    m4 = _max_load_case(df, 3)  # Case4
    top = max(m3, m4)
    return top * 1.02


def _style_xtract_combined_figure():
    plt.xlabel("Axial Strain", fontsize=12)
    plt.ylabel("Axial Load (kN)", fontsize=12)
    plt.legend()
    plt.grid(True)


def plot_xtract_case1_to_case4(df):
    """Xtract：case1～case4 綜合圖（與原 Xtract data.png 內容相同）。"""
    plt.figure(figsize=(10, 6))
    for label, (col_p, col_strain) in zip(LABELS, COL_PAIRS):
        strain, load = load_xtract_series(df, col_p, col_strain)
        plt.plot(strain, load, label=label, linewidth=2)
    plt.title("Monotonic Axial Response of RC Columns (Xtract)", fontsize=14)
    _style_xtract_combined_figure()
    plt.savefig(OUT_PATH_CASE1_TO_CASE4, dpi=300, bbox_inches="tight")
    plt.savefig(OUT_PATH, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"已儲存：{OUT_PATH_CASE1_TO_CASE4}")
    print(f"已儲存：{OUT_PATH}")


def plot_xtract_case0_to_case4(df_all):
    """Xtract：case0～case4 綜合圖（case0 來自工作表 case0，case1～4 來自 ALL）。"""
    df0 = pd.read_excel(EXCEL_PATH, sheet_name="case0", header=None)
    strain0, load0 = load_case0_series(df0)
    plt.figure(figsize=(10, 6))
    plt.plot(strain0, load0, label="case0", linewidth=2)
    for label, (col_p, col_strain) in zip(LABELS, COL_PAIRS):
        strain, load = load_xtract_series(df_all, col_p, col_strain)
        plt.plot(strain, load, label=label, linewidth=2)
    plt.title("Monotonic Axial Response of RC Columns (Xtract)", fontsize=14)
    _style_xtract_combined_figure()
    plt.savefig(OUT_PATH_CASE0_TO_CASE4, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"已儲存：{OUT_PATH_CASE0_TO_CASE4}")


def plot_case0():
    """工作表 case0：單條 P–strain 曲線，樣式與 Xtract data.png 一致。"""
    df0 = pd.read_excel(EXCEL_PATH, sheet_name="case0", header=None)
    strain, load = load_case0_series(df0)
    plt.figure(figsize=(10, 6))
    plt.plot(strain, load, label="case0", linewidth=2, color=COLOR_XTRACT)
    plt.xlabel("Axial Strain", fontsize=12)
    plt.ylabel("Axial Load (kN)", fontsize=12)
    plt.title("Monotonic Axial Response of RC Columns (case0)", fontsize=14)
    plt.legend()
    plt.grid(True)
    plt.savefig(OUT_PATH_CASE0, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"已儲存：{OUT_PATH_CASE0}")


def plot_comparisons(df):
    COMPARISON_DIR.mkdir(parents=True, exist_ok=True)
    y_top = comparison_ylim_top_from_case34(df)
    for idx, (label, (col_p, col_strain)) in enumerate(zip(LABELS, COL_PAIRS), start=1):
        strain_x, load_x = load_xtract_series(df, col_p, col_strain)
        ospath = RRRRRR_DIR / f"PDeltaData_{label}.txt"
        strain_o, load_o = load_opensees_txt(ospath)

        plt.figure(figsize=(10, 6))
        plt.plot(strain_x, load_x, color=COLOR_XTRACT, linewidth=2, label="Xtract")
        plt.plot(strain_o, load_o, color=COLOR_OPENSEES, linewidth=2, label="OpenSees")
        plt.xlabel("Axial Strain", fontsize=12)
        plt.ylabel("Axial Load (kN)", fontsize=12)
        plt.title(f"{label} — Xtract vs OpenSees", fontsize=14)
        plt.legend()
        plt.grid(True)
        plt.ylim(0, y_top)
        out = COMPARISON_DIR / f"case{idx}_comparison.png"
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"已儲存：{out}")


def plot_op_vs_xt_all_cases(df):
    """Case1～4：OpenSees 與 Xtract 共八條曲線同圖；同 Case 同色，OpenSees 實線、Xtract 虛線。"""
    plt.figure(figsize=(10, 6))
    y_max = 0.0
    for i, (label, (col_p, col_strain)) in enumerate(zip(LABELS, COL_PAIRS)):
        color = CASE_COLORS[i]
        strain_x, load_x = load_xtract_series(df, col_p, col_strain)
        ospath = RRRRRR_DIR / f"PDeltaData_{label}.txt"
        strain_o, load_o = load_opensees_txt(ospath)
        y_max = max(y_max, float(np.max(load_x)), float(np.max(load_o)))
        plt.plot(
            strain_o,
            load_o,
            color=color,
            linewidth=2,
            linestyle="-",
            label=f"OpenSees {label}",
        )
        plt.plot(
            strain_x,
            load_x,
            color=color,
            linewidth=2,
            linestyle="--",
            label=f"Xtract {label}",
        )
    plt.xlabel("Axial Strain", fontsize=12)
    plt.ylabel("Axial Load (kN)", fontsize=12)
    plt.title("Monotonic Axial Response of RC Columns (OpenSees vs Xtract)", fontsize=14)
    plt.legend(loc="best", fontsize=8, ncol=2)
    plt.grid(True)
    plt.ylim(0, y_max * 1.02)
    plt.savefig(OUT_PATH_OP_VS_XT, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"已儲存：{OUT_PATH_OP_VS_XT}")


def main():
    df = pd.read_excel(EXCEL_PATH, sheet_name="ALL", header=None)
    plot_xtract_case1_to_case4(df)
    plot_xtract_case0_to_case4(df)
    plot_case0()
    plot_comparisons(df)
    plot_op_vs_xt_all_cases(df)


if __name__ == "__main__":
    main()
