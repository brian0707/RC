# Introduction ------------------------------------------------------------------------------------------------------------------------------------------------#
# HW3：RC 柱 Pushover；材料組合 1–4（Concrete01／Concrete04 × Steel01／Steel02）。
# 混凝土強度目前為 28／35／42 MPa（已暫時移除 70 MPa；若要加回請指示「幫我加入70MPa混凝土」）。

# i import DLC ------------------------------------------------------------------------------------------------------------------------------------------------#
import math
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import openseespy.opensees as ops

# unit_SI、rebar：以 unit_SI.py 實際路徑為準
def _ensure_vendor_path() -> None:
    here = Path(__file__).resolve().parent
    candidates = [
        here,
        here / "宥宏",
        here.parent / "宥宏",
        here.parent / "HW1" / "宥宏",
        here.parent.parent / "HW1" / "宥宏",
    ]
    for d in candidates:
        if (d / "unit_SI.py").is_file():
            sys.path.insert(0, str(d))
            return
    for root in (here.parent, here.parent.parent):
        hw1 = root / "HW1"
        if hw1.is_dir():
            for found in hw1.rglob("unit_SI.py"):
                if found.name == "unit_SI.py":
                    sys.path.insert(0, str(found.parent))
                    return
    raise ImportError(
        "找不到 unit_SI.py。預期位置：腳本同層、同層「宥宏」、或專案內「HW1/宥宏」。"
    )


_ensure_vendor_path()

# 與腳本同層「results_test」：Pushover 曲線、HW3_case{N}_Data 細部輸出
_DATA_DIR = Path(__file__).resolve().parent
DATA_RESULTS_DIR = _DATA_DIR / "results_test"

# ii User-define function import ------------------------------------------------------------------------------------------------------------------------------#

# Unit, Rebar
# 這章節程式碼定義了使用毫米（unit_SI.unit_SI.mm）、千牛頓（unit_SI.kN）和秒（sec）作為單位的系統。
# 包含了長度、面積、力和應力的轉換關係，並設置了一些常數，如圓周率（PI）和重力加速度（g）、極大值（Ubig）極小值（Usmall）等。
# 可用單位：
    # 長度：unit_SI.unit_SI.mm, inch, ft, cm, m
    # 面積：inch2, unit_SI.unit_SI.mm2, cm2, m2
    # 慣性矩：inch4
    # 力：kip, N, MN, tf
    # 應力：unit_SI.MPa, GPa, ksi, psi, lbf, pcf
import unit_SI
unit_SI.define_units()

# 可使用鋼筋號數：
    # D6, D10, D13, D16, D19, D22, D25, D29, D32, D36, D38, D43, D50, D57, D64
# 可使用鋼筋面積：
    # A6, A10, A13, A16, A19, A22, A25, A29, A32, A36, A38, A43, A50, A57, A64
import rebar
rebar.define_rebar()

# 柱有效長度（與 Pushover 建模節點間距一致）
COL_EFFECTIVE_LENGTH = 3200 * unit_SI.mm



# iii User-define function ------------------------------------------------------------------------------------------------------------------------------------#

# 計算圍束混凝土強度
def func_get_mander_params(fco, fyh, BCol, HCol, clear_cover, n_long_y, n_long_z, d_long, d_tie, s_spacing):
    """
    計算 Mander 受限混凝土參數
    fco: 未圍束混凝土強度 (MPa)
    fyh: 箍筋屈服強度 (MPa)
    B, H: 斷面寬、高 (mm)
    cover: 淨保護層厚度 (mm)
    n_long_y, n_long_z: Y軸與 Z軸向的縱筋支數
    d_long: 縱筋直徑 (mm)
    d_tie: 箍筋直徑 (mm)
    s_spacing: 箍筋垂直間距 (mm)
    """
    # 1. 計算核心尺寸 (至箍筋中心線) [cite: 132, 233]
    bc = BCol - 2 * clear_cover - d_tie
    dc = HCol - 2 * clear_cover - d_tie
    Acc = bc * dc # 核心總面積 [cite: 132, 234]
    
    # 2. 計算縱筋淨距之平方和 (w'i^2) 以求 ke [cite: 226, 231, 236]
    # 假設鋼筋均勻分佈，計算相鄰鋼筋間的淨距離
    wy = (bc-(n_long_y*d_long)) / (n_long_y-1)
    wz = (dc-(n_long_z*d_long)) / (n_long_z-1)
    sum_wi_sq = (n_long_y - 1) * 2 * (wy**2) + (n_long_z - 1) * 2 * (wz**2)
    
    # 3. 計算圍束有效係數 ke [cite: 137, 236, 577]
    s_prime = s_spacing - d_tie # 箍筋間淨距 [cite: 171, 726]
    rho_cc = ( (n_long_y*2 + (n_long_z-2)*2) * (0.25 * math.pi * d_long**2) ) / Acc # 縱筋比 [cite: 142, 143]
    
    ke = (1 - sum_wi_sq / (6 * bc * dc)) * (1 - s_prime / (2 * bc)) * (1 - s_prime / (2 * dc)) / (1 - rho_cc)
    
    # 4. 計算側向圍束應力 fl' [cite: 134, 219, 255]
    # 假設為 4 肢箍筋配置 (你的圖面配置)
    Asx = 4 * (0.25 * math.pi * d_tie**2) # X向總箍筋面積
    Asy = 4 * (0.25 * math.pi * d_tie**2) # Y向總箍筋面積
    flx = (Asx / (s_spacing * dc)) * fyh # 全域側壓力 [cite: 248, 249]
    fly = (Asy / (s_spacing * bc)) * fyh
    
    flx_prime = ke * flx # 有效側壓力 [cite: 255, 256]
    fly_prime = ke * fly
    
    # 5. 計算受限強度提升 fcc' [cite: 91, 266, 267]
    # 使用 Mander 針對不同兩向側壓力的圖解簡化公式或威廉-瓦克準則
    # 此處取兩向平均值進行強度計算 (對稱斷面常用)
    fl_avg = (flx_prime + fly_prime) / 2
    fcc = fco * (-1.254 + 2.254 * math.sqrt(1 + 7.94 * fl_avg / fco) - 2 * fl_avg / fco)
    
    # 6. 計算對應應變 ecc [cite: 51, 94, 579]
    eco = 0.002 # 未圍束峰值應變 [cite: 117]
    ecc = eco * (1 + 5 * (fcc / fco - 1))
    
    # 7. 計算形狀參數 r 與 Ec [cite: 118, 120, 581]
    fco_MPa = fco / unit_SI.MPa # 轉換為 MPa 單位
    Ec = 5000 * math.sqrt(fco_MPa)*unit_SI.MPa # MPa [cite: 120]
    Esec = fcc / ecc
    r = Ec / (Ec - Esec)
    
    return fcc, ecc, Ec, r


def _ec_sqrt_fc_si(fco_abs_mpa: float) -> float:
    """與 func_get_mander_params 一致：Ec = 5000√f'c (MPa)。"""
    return 5000.0 * math.sqrt(max(fco_abs_mpa, 1e-9)) * unit_SI.MPa


def _ec_popovics_match_mander_si(fcc_abs_mpa: float) -> float:
    """Concrete04：OpenSees 手冊建議 Ec = 57000√|fc|（fc 以 psi）時壓縮包絡等同 Mander (1988)。"""
    fc_psi = max(fcc_abs_mpa, 1e-9) * 145.0377377
    ec_psi = 57000.0 * math.sqrt(fc_psi)
    ec_mpa = ec_psi / 145.0377377
    return ec_mpa * unit_SI.MPa


def _confined_crushing_strain_eps_cu(fpc_uc_signed: float) -> float:
    """延續原腳本 Concrete01 圍束混凝土峰後應變慣例（負值）。"""
    f_abs_mpa = abs(fpc_uc_signed) / unit_SI.MPa
    if f_abs_mpa <= 35.0 + 1e-9:
        return -0.02
    return -0.015


def _concrete04_uc_unconfined(fpc_uc_signed: float, mat_tag: int) -> None:
    """保護層：未圍束 Concrete04（Popovics），峰值應變 −0.002、極限壓應變 −0.004。"""
    f_abs_mpa = abs(fpc_uc_signed) / unit_SI.MPa
    fc = fpc_uc_signed
    ec = -0.002
    ecu = -0.004
    Ec = _ec_sqrt_fc_si(f_abs_mpa)
    ft = 0.05 * math.sqrt(max(f_abs_mpa, 1e-9)) * unit_SI.MPa
    et = 0.0001
    beta = 0.1
    ops.uniaxialMaterial("Concrete04", mat_tag, fc, ec, ecu, Ec, ft, et, beta)


def _concrete04_cc_mander_popovics(
    fcc_pos: float,
    ecc_pos: float,
    eps_cu_neg: float,
    mat_tag: int,
) -> None:
    """核心圍束：Concrete04，Ec 取 Mander 匹配；fcc、ecc 來自 func_get_mander_params。"""
    f_abs_mpa = abs(fcc_pos) / unit_SI.MPa
    fc = -abs(fcc_pos)
    ec = -abs(ecc_pos)
    ecu = eps_cu_neg if eps_cu_neg < ec - 1e-15 else ec * 1.5
    Ec = _ec_popovics_match_mander_si(f_abs_mpa)
    ft = 0.05 * math.sqrt(max(f_abs_mpa, 1e-9)) * unit_SI.MPa
    et = 0.0001
    beta = 0.1
    ops.uniaxialMaterial("Concrete04", mat_tag, fc, ec, ecu, Ec, ft, et, beta)


# 材料組合（每次 wipe 後呼叫一次）
# 1 baseline：Concrete01（保護層、核心）+ Steel01
# 2：Concrete01（保護層）+ Concrete04 核心（Mander fcc/ecc + Popovics Ec）+ Steel01
# 3：Concrete04（保護層未圍束）+ Concrete04 核心 + Steel02（Menegotto–Pinto／Filippou）
# 4：Concrete01（保護層）+ Concrete04 核心 + Steel02
def func_MatDef(BCol, HCol, D_Longi, D_StirUp, material_combo: int = 1):
    fyh = 420.0 * unit_SI.MPa
    clear_cover = 40.0 * unit_SI.mm
    n_long_y = 5
    n_long_z = 5
    s_spacing = 100.0 * unit_SI.mm

    ID28MPaUC = 4001
    ID35MPaUC = 4002
    ID42MPaUC = 4003
    ID28MPaCC = 4101
    ID35MPaCC = 4102
    ID42MPaCC = 4103

    fpc28MPaUC = -28.0 * unit_SI.MPa
    fpc35MPaUC = -35.0 * unit_SI.MPa
    fpc42MPaUC = -42.0 * unit_SI.MPa

    uc_specs = (
        (ID28MPaUC, fpc28MPaUC),
        (ID35MPaUC, fpc35MPaUC),
        (ID42MPaUC, fpc42MPaUC),
    )

    # ----- 保護層（4001–4003；若加回 70 MPa 見檔頭說明）-----
    for mid, fpc in uc_specs:
        if material_combo == 3:
            _concrete04_uc_unconfined(fpc, mid)
        else:
            ops.uniaxialMaterial(
                "Concrete01",
                mid,
                fpc,
                -0.002,
                0.0,
                -0.004,
            )

    # ----- 核心圍束（4101–4103），Mander 參數 -----
    fcc_28, ecc_28, Ec_28, r_28 = func_get_mander_params(
        -fpc28MPaUC, fyh, BCol, HCol, clear_cover, n_long_y, n_long_z, D_Longi, D_StirUp, s_spacing
    )
    fcc_35, ecc_35, Ec_35, r_35 = func_get_mander_params(
        -fpc35MPaUC, fyh, BCol, HCol, clear_cover, n_long_y, n_long_z, D_Longi, D_StirUp, s_spacing
    )
    fcc_42, ecc_42, Ec_42, r_42 = func_get_mander_params(
        -fpc42MPaUC, fyh, BCol, HCol, clear_cover, n_long_y, n_long_z, D_Longi, D_StirUp, s_spacing
    )

    cc_blocks = (
        (ID28MPaCC, fcc_28, ecc_28, fpc28MPaUC),
        (ID35MPaCC, fcc_35, ecc_35, fpc35MPaUC),
        (ID42MPaCC, fcc_42, ecc_42, fpc42MPaUC),
    )

    for cc_tag, fcc_i, ecc_i, fpc_uc_i in cc_blocks:
        if material_combo == 1:
            eps_cu_mander01 = -0.02 if abs(fpc_uc_i) <= 35.0 * unit_SI.MPa + 1e-9 else -0.015
            ops.uniaxialMaterial(
                "Concrete01",
                cc_tag,
                -fcc_i,
                -ecc_i,
                -0.2 * fcc_i,
                eps_cu_mander01,
            )
        else:
            ecu = _confined_crushing_strain_eps_cu(fpc_uc_i)
            _concrete04_cc_mander_popovics(fcc_i, ecc_i, ecu, cc_tag)

    # ----- 鋼筋 4201 SD420、4202 SD550 -----
    IDSD420 = 4201
    IDSD550 = 4202
    fySD420 = 420 * unit_SI.MPa
    fySD550 = 550 * unit_SI.MPa
    Es = 200 * unit_SI.GPa
    b_hard = 0.018

    if material_combo in (1, 2):
        ops.uniaxialMaterial("Steel01", IDSD420, fySD420, Es, b_hard)
        ops.uniaxialMaterial("Steel01", IDSD550, fySD550, Es, b_hard)
    else:
        R0 = 18.0
        cR1 = 0.925
        cR2 = 0.15
        ops.uniaxialMaterial("Steel02", IDSD420, fySD420, Es, b_hard, R0, cR1, cR2)
        ops.uniaxialMaterial("Steel02", IDSD550, fySD550, Es, b_hard, R0, cR1, cR2)

    return fcc_28, ecc_28, fcc_35, ecc_35, fcc_42, ecc_42

# 斷面定義參數
def func_sectionDef(IDSec, IDUC, IDCC, IDMat_Longi, BCol, HCol, DLongi, ALongi, DStirUp):
    # Y軸為橫軸，Z軸為縱軸

    # Define Section
    ops.section('Fiber', IDSec)


    # cover
    cover = 50*unit_SI.mm + DStirUp + DLongi/2

    # edge
    yE1 = -BCol/2
    yE2 = BCol/2
    zE1 = -HCol/2
    zE2 = HCol/2
    # core
    yC1 = yE1 + cover - DLongi/2
    yC2 = -yC1
    zC1 = zE1 + cover - DLongi/2
    zC2 = -zC1

    # Longitudinal
    spacing = (BCol-2*cover)/4
    yL1 = yE1 + cover + spacing*0
    yL2 = yE1 + cover + spacing*1
    yL3 = yE1 + cover + spacing*2
    yL4 = yE1 + cover + spacing*3
    yL5 = yE1 + cover + spacing*4

    zL1 = zE1 + cover + spacing*0
    zL2 = zE1 + cover + spacing*1
    zL3 = zE1 + cover + spacing*2
    zL4 = zE1 + cover + spacing*3
    zL5 = zE1 + cover + spacing*4

    # 核心混凝土
    ops.patch('rect', IDCC, 20, 20, yC1, zC1, yC2, zC2)

    # 保護層混凝土 (上, 下, 左, 下)
    ops.patch('rect', IDUC, 24, 2, yE1, zC2, yE2, zE2)
    ops.patch('rect', IDUC, 24, 2, yE1, zE1, yE2, zC1)
    ops.patch('rect', IDUC, 2, 24,  yE1, zC1, yC1, zC2)
    ops.patch('rect', IDUC, 2, 24,  yC2, zC1, yE2, zC2)

    # Longitudinal bars (左至右直欄：5, 2, 2, 2, 5)
    ops.layer('straight', IDMat_Longi, 5, ALongi, yL1, zL1, yL1, zL5)
    ops.layer('straight', IDMat_Longi, 2, ALongi, yL2, zL1, yL2, zL5)
    ops.layer('straight', IDMat_Longi, 2, ALongi, yL3, zL1, yL3, zL5)
    ops.layer('straight', IDMat_Longi, 2, ALongi, yL4, zL1, yL4, zL5)
    ops.layer('straight', IDMat_Longi, 5, ALongi, yL5, zL1, yL5, zL5)


    print("Define Fiber Section")
    print("-------------------------------------------------------------------------------------")


def _static_analyze_step_pushover() -> int:
    """側推一步靜力分析：Newton 不收斂時改試 NewtonLineSearch（纖維軟化後較常需要）。"""
    ok = ops.analyze(1)
    if ok == 0:
        return 0
    ops.algorithm("NewtonLineSearch")
    ok = ops.analyze(1)
    ops.algorithm("Newton")
    return ok


def func_run_pushover_analysis(
    case_num,
    id_sec,
    b_col,
    h_col,
    id_uc,
    id_cc,
    id_mat_longi,
    d_longi,
    a_longi,
    d_stirrup,
    material_combo: int = 1,
    results_root=None,
):
    """重力載重（軸壓 0.2×對應案別之圍束混凝土 f'_cc×斷面積）後，位移控制側推（頂點水平自由度 1）。
    results_root：若指定，輸出至該資料夾底下 HW3_case{n}_Data（供四種材料組分分開存放）。
    回傳 (curve_xy, fcc_28, ...)：curve 兩欄為 [頂點水平位移 mm, 估算基底剪力 kN]。"""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)

    root = Path(results_root) if results_root is not None else DATA_RESULTS_DIR
    root.mkdir(parents=True, exist_ok=True)
    data_dir = str(root / f"HW3_case{case_num}_Data")
    os.makedirs(data_dir, exist_ok=True)

    l_col = COL_EFFECTIVE_LENGTH

    ops.node(3001, 0, 0)
    ops.node(3002, 0, l_col / 2)
    ops.node(3003, 0, l_col)
    ops.fix(3001, 1, 1, 1)

    fcc_28, ecc_28, fcc_35, ecc_35, fcc_42, ecc_42 = func_MatDef(
        b_col, h_col, d_longi, d_stirrup, material_combo
    )

    ag = b_col * h_col
    fcc_by_case = {
        1: fcc_28,
        2: fcc_28,
        3: fcc_35,
        4: fcc_35,
        5: fcc_42,
        6: fcc_42,
    }
    fcc_case = fcc_by_case[case_num]
    p_gravity = -0.2 * abs(fcc_case) * ag

    func_sectionDef(id_sec, id_uc, id_cc, id_mat_longi, b_col, h_col, d_longi, a_longi, d_stirrup)

    col_transf_tag = 6901
    ops.geomTransf("Linear", col_transf_tag)
    num_int_pts = 5
    id_integration_tag = 1
    max_iter = 50
    tol = 1e-7
    ops.beamIntegration("Lobatto", id_integration_tag, id_sec, num_int_pts)
    ops.element(
        "forceBeamColumn",
        6001,
        3001,
        3002,
        col_transf_tag,
        id_integration_tag,
        "iter",
        max_iter,
        tol,
    )
    ops.element(
        "forceBeamColumn",
        6002,
        3002,
        3003,
        col_transf_tag,
        id_integration_tag,
        "iter",
        max_iter,
        tol,
    )

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3003, 0.0, p_gravity, 0.0)

    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormUnbalance", 1e-7, 80)
    ops.algorithm("Newton")
    n_grav = 12
    ops.integrator("LoadControl", 1.0 / n_grav)
    ops.analysis("Static")
    for gi in range(n_grav):
        ok = ops.analyze(1)
        if ok != 0:
            print(f"Case {case_num}: 重力載重不收斂於步 {gi}")
            break

    ops.loadConst("-time", 0.0)

    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(3003, 1.0, 0.0, 0.0)

    du = 2.0 * unit_SI.mm
    ops.integrator("DisplacementControl", 3003, 1, du)
    ops.analysis("Static")

    plot_data = []
    max_steps = 1200
    for i in range(max_steps):
        ok = _static_analyze_step_pushover()
        if ok != 0:
            print(
                f"Case {case_num} Pushover 停止於第 {i} 步（可能軟化或不收斂）"
                f"（頂點位移約 {ops.nodeDisp(3003, 1):.2f} mm）"
            )
            break
        # OpenSees 須先 reactions()，nodeReaction 才有值（否則常為全 0，圖會貼在橫軸）
        ops.reactions()
        ux = ops.nodeDisp(3003, 1)
        vb = -ops.nodeReaction(3001, 1)
        plot_data.append([ux, vb])

    if plot_data:
        curve_path = os.path.join(data_dir, "Pushover_disp_vs_base_shear.txt")
        np.savetxt(
            curve_path,
            np.asarray(plot_data, dtype=float),
            fmt="%12.6e",
            header="top_disp_x_mm base_shear_kN",
            comments="",
        )

    print(f"Case {case_num} Pushover 完成，輸出於 {data_dir}.")
    return (
        np.asarray(plot_data, dtype=float),
        fcc_28,
        ecc_28,
        fcc_35,
        ecc_35,
        fcc_42,
        ecc_42,
    )


# HW3 main program --------------------------------------------------------------------------------------------------------------------------------------------#


DLongi = rebar.D25
ALongi = rebar.A25
DStirUp = rebar.D13

BCol = 600 * unit_SI.mm
HCol = 600 * unit_SI.mm

ID28MPaUC = 4001
ID35MPaUC = 4002
ID42MPaUC = 4003
ID28MPaCC = 4101
ID35MPaCC = 4102
ID42MPaCC = 4103
IDSD420 = 4201
IDSD550 = 4202
IDsec1 = 5001
IDsec2 = 5002
IDsec3 = 5003
IDsec4 = 5004
IDsec5 = 5005
IDsec6 = 5006

# --- Pushover：材料組合 × 六案（28／35／42 MPa × SD420／SD550）---
# 改為 RUN_MATERIAL_COMBOS = (2,) 等可只跑單一組合。
RUN_MATERIAL_COMBOS = (1, 2, 3, 4)
MATERIAL_COMBO_TITLE = {
    1: "Concrete01+Concrete01+Steel01 (baseline)",
    2: "Concrete01(uc)+Concrete04(cc,Mander Popovics Ec)+Steel01",
    3: "Concrete04(uc)+Concrete04(cc)+Steel02",
    4: "Concrete01(uc)+Concrete04(cc)+Steel02",
}
COMBO_DIRNAME = {
    1: "combo01_baseline_C01_S01",
    2: "combo02_C01_C04_S01",
    3: "combo03_C04_C04_S02",
    4: "combo04_C01_C04_S02",
}
# 各材料組合一律輸出「初始段放大圖」之橫軸上限（mm）；若該組最大位移較短則自動縮為資料上限
PUSHOVER_ZOOM_INITIAL_X_MM = 600.0

CASE_RUNS = (
    (1, IDsec1, ID28MPaUC, ID28MPaCC, IDSD420, "Case1_28MPa_SD420"),
    (2, IDsec2, ID28MPaUC, ID28MPaCC, IDSD550, "Case2_28MPa_SD550"),
    (3, IDsec3, ID35MPaUC, ID35MPaCC, IDSD420, "Case3_35MPa_SD420"),
    (4, IDsec4, ID35MPaUC, ID35MPaCC, IDSD550, "Case4_35MPa_SD550"),
    (5, IDsec5, ID42MPaUC, ID42MPaCC, IDSD420, "Case5_42MPa_SD420"),
    (6, IDsec6, ID42MPaUC, ID42MPaCC, IDSD550, "Case6_42MPa_SD550"),
)
# func_run_pushover_analysis 回傳項索引：1→fcc_28、3→fcc_35、5→fcc_42
IDX_RET_FCC_FOR_CASE = {1: 1, 2: 1, 3: 3, 4: 3, 5: 5, 6: 5}

DATA_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
Ag = BCol * HCol

for combo_id in RUN_MATERIAL_COMBOS:
    combo_root = DATA_RESULTS_DIR / COMBO_DIRNAME[combo_id]
    combo_root.mkdir(parents=True, exist_ok=True)
    push_dir = combo_root / "PushoverData"
    push_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 80)
    print(f"材料組合 {combo_id}：{MATERIAL_COMBO_TITLE[combo_id]}")
    print(f"輸出根目錄：{combo_root}")
    print("=" * 80)

    results = {}
    for cn, idsec, iduc, idcc, idsteel, label in CASE_RUNS:
        ret = func_run_pushover_analysis(
            cn,
            idsec,
            BCol,
            HCol,
            iduc,
            idcc,
            idsteel,
            DLongi,
            ALongi,
            DStirUp,
            material_combo=combo_id,
            results_root=combo_root,
        )
        results[label] = ret

    print("\n" + "-" * 80)
    print(
        f"[組合 {combo_id}] Pushover 重力階段軸力（kN；壓力為負）："
        "P = 0.2×|f'_cc,case|×Ag"
    )
    print("-" * 80)
    for cn, _, _, _, _, lab in CASE_RUNS:
        ret = results[lab]
        fi = IDX_RET_FCC_FOR_CASE[cn]
        fcc_case = ret[fi]
        _P_kN = -0.2 * abs(fcc_case) * Ag
        _fcc_MPa = abs(fcc_case) / unit_SI.MPa
        print(
            f"  Case {cn} ({lab}): P = {_P_kN:.6f} kN"
            f"  |  f'cc = {_fcc_MPa:.2f} MPa  |  |P| = {abs(_P_kN):.6f} kN"
        )
    print("-" * 80 + "\n")

    print(f"[組合 {combo_id}] 正在輸出 Pushover 彙整 TXT …")
    print(f"PushoverData：{push_dir}")
    for label, (
        data,
        fcc28,
        ecc28,
        fcc35,
        ecc35,
        fcc42,
        ecc42,
    ) in results.items():
        if len(data) > 0:
            fcc_28MPa = fcc28 / unit_SI.MPa
            fcc_35MPa = fcc35 / unit_SI.MPa
            fcc_42MPa = fcc42 / unit_SI.MPa
            hdr = (
                f"material_combo={combo_id}\n"
                f"f'cc28MPa = {fcc_28MPa:.2f} MPa, ecc28MPa = {ecc28:.4f}\n"
                f"f'cc35MPa = {fcc_35MPa:.2f} MPa, ecc35MPa = {ecc35:.4f}\n"
                f"f'cc42MPa = {fcc_42MPa:.2f} MPa, ecc42MPa = {ecc42:.4f}\n"
                "top_disp_x_mm base_shear_kN"
            )
            fname = str(push_dir / f"PushoverData_{label}.txt")
            np.savetxt(fname, data, fmt="%12.6e", header=hdr, comments="")
            print(f"已儲存: {fname}")

    curves_xy = []
    for label, val in results.items():
        data = val[0]
        if len(data) > 0:
            curves_xy.append((label, np.asarray(data, dtype=float)))

    plt.figure(figsize=(10, 6))
    for label, xy in curves_xy:
        plt.plot(xy[:, 0], xy[:, 1], label=label, linewidth=2)

    plt.xlabel("Top lateral displacement (mm)", fontsize=12)
    plt.ylabel("Base shear (kN)", fontsize=12)
    plt.title(
        f"Pushover — combo {combo_id}: {MATERIAL_COMBO_TITLE[combo_id]}",
        fontsize=10,
    )
    plt.legend(fontsize=7)
    plt.grid(True)

    _push_png = combo_root / "Pushover_all_cases.png"
    plt.savefig(_push_png, dpi=300, bbox_inches="tight")
    print(f"[組合 {combo_id}] 已輸出合圖：{_push_png}")
    plt.close()

    if curves_xy:
        xmax_each = [float(np.max(xy[:, 0])) for _, xy in curves_xy]
        global_max_x = max(xmax_each)
        short_min_x = min(xmax_each)

        x_ini_hi = float(min(PUSHOVER_ZOOM_INITIAL_X_MM, global_max_x))
        if x_ini_hi > 15.0:
            plt.figure(figsize=(9, 5.5))
            for label, xy in curves_xy:
                m = xy[:, 0] <= x_ini_hi
                plt.plot(xy[m, 0], xy[m, 1], label=label, linewidth=2)
            plt.xlabel("Top lateral displacement (mm)", fontsize=12)
            plt.ylabel("Base shear (kN)", fontsize=12)
            plt.title(
                f"Pushover initial zoom (0–{x_ini_hi:.0f} mm) — combo {combo_id}",
                fontsize=10,
            )
            plt.xlim(0.0, x_ini_hi)
            plt.legend(fontsize=7)
            plt.grid(True)
            _ini_png = combo_root / "Pushover_all_cases_zoom_initial.png"
            plt.savefig(_ini_png, dpi=300, bbox_inches="tight")
            print(f"[組合 {combo_id}] 已輸出初始段放大圖：{_ini_png}")
            plt.close()

        # SD420 + 較高 f'c 常提早不收斂，曲線僅數 ten mm，在全尺度圖上像「消失」—另存更窄原點圖
        if global_max_x > 80.0 and short_min_x < 0.25 * global_max_x:
            x_zoom_hi = float(min(150.0, max(70.0, short_min_x * 4.0 + 30.0)))
            plt.figure(figsize=(9, 5.5))
            for label, xy in curves_xy:
                m = xy[:, 0] <= x_zoom_hi
                plt.plot(xy[m, 0], xy[m, 1], label=label, linewidth=2)
            plt.xlabel("Top lateral displacement (mm)", fontsize=12)
            plt.ylabel("Base shear (kN)", fontsize=12)
            plt.title(
                f"Pushover origin zoom (0–{x_zoom_hi:.0f} mm) — combo {combo_id}",
                fontsize=10,
            )
            plt.xlim(0.0, x_zoom_hi)
            plt.legend(fontsize=7)
            plt.grid(True)
            _zoom_png = combo_root / "Pushover_all_cases_zoom_origin.png"
            plt.savefig(_zoom_png, dpi=300, bbox_inches="tight")
            print(f"[組合 {combo_id}] 已輸出原點放大圖（short curves）：{_zoom_png}")
            plt.close()
