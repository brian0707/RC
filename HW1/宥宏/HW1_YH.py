# Introduction ------------------------------------------------------------------------------------------------------------------------------------------------#


# i import DLC ------------------------------------------------------------------------------------------------------------------------------------------------#
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import os
import math

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

# 4. 定義材料
def func_MatDef(BCol, HCol, D_Longi, D_StirUp):
    # 材料定義參數
    # unconfined concrete (28MPa)
    ID28MPaUC = 4001
    fpc28MPaUC = -28.0 * unit_SI.MPa
    epsc028MPaUC = -0.002
    fpcu28MPaUC = 0.0
    epsU28MPaUC = -0.004
    ops.uniaxialMaterial("Concrete01", ID28MPaUC, fpc28MPaUC, epsc028MPaUC, fpcu28MPaUC, epsU28MPaUC)


    # unconfined concrete (70MPa)
    ID70MPaUC = 4002
    fpc70MPaUC = -70 * unit_SI.MPa
    epsc070MPaUC = -0.002
    fpcu70MPaUC = 0.0
    epsU70MPaUC = -0.004
    ops.uniaxialMaterial("Concrete01", ID70MPaUC, fpc70MPaUC, epsc070MPaUC, fpcu70MPaUC, epsU70MPaUC)


    # confined concrete (28MPa)
    ID28MPaCC = 4101
    fyh = 420.0 * unit_SI.MPa
    clear_cover = 40.0 * unit_SI.mm
    n_long_y = 5
    n_long_z = 5
    s_spacing = 100.0 * unit_SI.mm
    fcc_28, ecc_28, Ec_28, r_28 = func_get_mander_params(-fpc28MPaUC, fyh, BCol, HCol, clear_cover, n_long_y, n_long_z, D_Longi, D_StirUp, s_spacing)
    # Concrete04: fpc, epsc, fpcu, epsu, ec, r
    # ops.uniaxialMaterial('Concrete04', ID28MPaCC, -fcc_28, -ecc_28, -0.2*fcc_28, -0.02, Ec_28, r_28)
    ops.uniaxialMaterial('Concrete01', ID28MPaCC, -fcc_28, -ecc_28, -0.2*fcc_28, -0.02)
    
    # confined concrete (70MPa)
    ID70MPaCC = 4102
    fyh = 420.0 * unit_SI.MPa
    clear_cover = 40.0 * unit_SI.mm
    n_long_y = 5
    n_long_z = 5
    s_spacing = 100.0 * unit_SI.mm
    fcc_70, ecc_70, Ec_70, r_70 = func_get_mander_params(-fpc70MPaUC, fyh, BCol, HCol, clear_cover, n_long_y, n_long_z, D_Longi, D_StirUp, s_spacing)
    # ops.uniaxialMaterial('Concrete04', ID70MPaCC, -fcc_70, -ecc_70, -0.2*fcc_70, -0.015, Ec_70, r_70)
    ops.uniaxialMaterial('Concrete01', ID70MPaCC, -fcc_70, -ecc_70, -0.2*fcc_70, -0.015)

    # rebar (SD420)
    IDSD420 = 4201
    fySD420 = 420 * unit_SI.MPa
    fuSD420 = fySD420*1.25
    EsSD420 = 200 * unit_SI.GPa
    EshSD420 = EsSD420/30
    bSD420 = 0.018
    R0SD420 = 0.13
    ops.uniaxialMaterial('ReinforcingSteel', IDSD420, fySD420, fuSD420, EsSD420, EshSD420, bSD420, R0SD420)

    # rebar (SD550)
    IDSD550 = 4202
    fySD550 = 550 * unit_SI.MPa
    fuSD550 = fySD550*1.25
    EsSD550 = 200 * unit_SI.GPa
    EshSD550 = EsSD550/30
    bSD550 = 0.018
    R0SD550 = 0.13
    ops.uniaxialMaterial('ReinforcingSteel', IDSD550, fySD550, fuSD550, EsSD550, EshSD550, bSD550, R0SD550)

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

# 執行單調軸壓分析的獨立函數
def func_monotonic_axial_analysis(LCol, handles, caseNum):
    """
    執行單調軸壓分析的獨立函數
    """
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(3003, 0.0, -1.0, 0.0) # 施加單位軸壓

    ops.system('BandGeneral')
    ops.numberer('Plain')
    ops.constraints('Plain')
    ops.integrator('DisplacementControl', 3003, 2, -0.04) # 每次推進 -0.04
    ops.algorithm('Newton')
    ops.analysis('Static')

    plot_data = []
    
    for i in range(1000):
        ok = ops.analyze(1)
        if ok == 0:
            # 讀取節點位移與元素力量
            force_ele1 = ops.eleResponse(6001, 'globalForce')
            force_ele2 = ops.eleResponse(6002, 'globalForce')
            
            d = [
                ops.nodeDisp(3001, 2), force_ele1[1],
                ops.nodeDisp(3002, 2), force_ele1[4],
                ops.nodeDisp(3003, 2), force_ele2[4]
            ]
            
            # 手動寫入檔案
            for j in range(6):
                handles[j].write(f"{d[j]:12.6e}\n")
            
            # 紀錄繪圖數據 [軸向應變, 軸向力量(kN)]
            plot_data.append([abs(d[4]) / LCol, abs(d[5])])
        else:
            print(f"Case {caseNum} 分析停止於第 {i} 步 (可能是壓碎或不收斂)")
            break
            
    return plot_data

# 主函數：執行柱子分析
def func_run_column_analysis(caseNum, IDSec, BCol, HCol, IDUC, IDCC, IDMat_Longi, DLongi, ALongi,  DStirUp):

    # 1. initialize system setting --------------------------------------------------------------------------------------------------------------------------------#
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    
    dataDir = os.path.join(OUTPUT_DIR, f"HW1_case{caseNum}_Data")
    if not os.path.exists(dataDir):
        os.makedirs(dataDir)


    # 2. Geometry -------------------------------------------------------------------------------------------------------------------------------------------------#
    LCol = 3200 *unit_SI.mm
    BCol =  600 * unit_SI.mm
    HCol =  600 * unit_SI.mm

    # 3. Node -----------------------------------------------------------------------------------------------------------------------------------------------------#
    ops.node(3001, 0, 0)
    ops.node(3002, 0, LCol/2)
    ops.node(3003, 0, LCol)
    ops.fix(3001, 1, 1, 1)

    # 4. Materials ------------------------------------------------------------------------------------------------------------------------------------------------#
    func_MatDef(BCol, HCol, DLongi, DStirUp)
        
    # 5. Section --------------------------------------------------------------------------------------------------------------------------------------------------#
    func_sectionDef(IDSec, IDUC, IDCC, IDMat_Longi, BCol, HCol, DLongi, ALongi, DStirUp)
    # 6. Elements -------------------------------------------------------------------------------------------------------------------------------------------------#
    
    ColTransfTag = 6901
    ColTransType = 'Linear'
    ops.geomTransf(ColTransType, ColTransfTag)

    numIntPts = 5
    IDIntegrationTag = 1
    maxIter = 50
    tol = 1e-7
    ops.beamIntegration('Lobatto', IDIntegrationTag, IDSec, numIntPts)

    # element('forceBeamColumn', eleTag, iNode, jNode, transfTag, integrationTag)
    ops.element('forceBeamColumn', 6001, 3001, 3002, ColTransfTag, IDIntegrationTag, 'iter', maxIter, tol)
    ops.element('forceBeamColumn', 6002, 3002, 3003, ColTransfTag, IDIntegrationTag, 'iter', maxIter, tol)

    # 7. recorder -------------------------------------------------------------------------------------------------------------------------------------------------#
    file_names = ['Disp1.out', 'Force1.out', 'Disp2.out', 'Force2.out', 'Disp3.out', 'Force3.out']
    handles = [open(os.path.join(dataDir, name), "w") for name in file_names]

    # 8. gravity analysis -----------------------------------------------------------------------------------------------------------------------------------------#

    
    # 9. Monotonic Axial Analysis----------------------------------------------------------------------------------------------------------------------------------#

    plot_data = func_monotonic_axial_analysis(LCol, handles, caseNum)
    
    for h in handles: h.close()
    print (f"Case {caseNum} analysis completed. Data saved in {dataDir}.")
    return np.array(plot_data)
    



# HW1 main program --------------------------------------------------------------------------------------------------------------------------------------------#

# 匯出數據目錄（與本腳本同層之 RRRRRR）
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(_SCRIPT_DIR, "RRRRRR")
os.makedirs(OUTPUT_DIR, exist_ok=True)

DLongi = rebar.D25
ALongi = rebar.A25
DStirUp = rebar.D13

BCol = 600 * unit_SI.mm
HCol = 600 * unit_SI.mm
LCol = 3200 * unit_SI.mm

case1 = 1 ; case2 = 2 ; case3 = 3 ; case4 = 4
ID28MPaUC = 4001 ; ID70MPaUC = 4002
ID28MPaCC = 4101 ; ID70MPaCC = 4102
IDSD420   = 4201 ; IDSD550   = 4202
IDsec1 = 5001 ; IDsec2 = 5002 ; IDsec3 = 5003 ; IDsec4 = 5004

# --- Terminal：各 Case 之 Concrete01 f_cu (=|fpcu|)、eps_cu (=|epsU|)（與 func_MatDef 一致）---
_fpc28_uc = -28.0 * unit_SI.MPa
_fpc70_uc = -70.0 * unit_SI.MPa
_fyh_mander = 420.0 * unit_SI.MPa
_clear_cc = 40.0 * unit_SI.mm
_nly = _nlz = 5
_s_sp = 100.0 * unit_SI.mm
_fcc28, _ecc28, _, _ = func_get_mander_params(
    -_fpc28_uc, _fyh_mander, BCol, HCol, _clear_cc, _nly, _nlz, DLongi, DStirUp, _s_sp
)
_fcc70, _ecc70, _, _ = func_get_mander_params(
    -_fpc70_uc, _fyh_mander, BCol, HCol, _clear_cc, _nly, _nlz, DLongi, DStirUp, _s_sp
)
_f_cu_uc_MPa = 0.0
_eps_cu_uc = 0.004
_f_cu_cc28_MPa = 0.2 * _fcc28 / unit_SI.MPa
_eps_cu_cc28 = 0.02
_f_cu_cc70_MPa = 0.2 * _fcc70 / unit_SI.MPa
_eps_cu_cc70 = 0.015
print("=" * 72)
print("Concrete01：f_cu = |fpcu|（殘餘壓力強度，MPa），eps_cu = |epsU|（極限壓應變）")
print("-" * 72)
for _cn, _lbl, _fc in [
    (1, "28 MPa 混凝土 / SD420", (_f_cu_uc_MPa, _eps_cu_uc, _f_cu_cc28_MPa, _eps_cu_cc28)),
    (2, "28 MPa 混凝土 / SD550", (_f_cu_uc_MPa, _eps_cu_uc, _f_cu_cc28_MPa, _eps_cu_cc28)),
    (3, "70 MPa 混凝土 / SD420", (_f_cu_uc_MPa, _eps_cu_uc, _f_cu_cc70_MPa, _eps_cu_cc70)),
    (4, "70 MPa 混凝土 / SD550", (_f_cu_uc_MPa, _eps_cu_uc, _f_cu_cc70_MPa, _eps_cu_cc70)),
]:
    _fuc, _euc, _fcc, _ecc = _fc
    print(f"Case {_cn} ({_lbl})")
    print(f"  未圍束(保護層): f_cu = {_fuc:.4f} MPa,  eps_cu = {_euc:.4f}")
    print(f"  圍束(核心):     f_cu = {_fcc:.4f} MPa,  eps_cu = {_ecc:.4f}")
print("=" * 72)

# --- 執行四種組合 ---
print("Running Column Analyses for All Cases...")
PDeltaCurve1 = func_run_column_analysis(case1, IDsec1, BCol, HCol,  ID28MPaUC, ID28MPaCC, IDSD420, DLongi, ALongi, DStirUp)
PDeltaCurve2 = func_run_column_analysis(case2, IDsec2, BCol, HCol,  ID28MPaUC, ID28MPaCC, IDSD550, DLongi, ALongi, DStirUp)
PDeltaCurve3 = func_run_column_analysis(case3, IDsec3, BCol, HCol,  ID70MPaUC, ID70MPaCC, IDSD420, DLongi, ALongi, DStirUp)
PDeltaCurve4 = func_run_column_analysis(case4, IDsec4, BCol, HCol,  ID70MPaUC, ID70MPaCC, IDSD550, DLongi, ALongi, DStirUp)


results = {
        'Case1_28MPa_SD420': PDeltaCurve1,
        'Case2_28MPa_SD550': PDeltaCurve2,
        'Case3_70MPa_SD420': PDeltaCurve3,
        'Case4_70MPa_SD550': PDeltaCurve4
    }
 
# --- 輸出四個 Case 的曲線數據到 TXT 檔 ---
print("正在輸出數據檔案...")
for label, data in results.items():
    if len(data) > 0:
        # 儲存為 txt 檔，第一欄是 Strain，第二欄是 Axial Load (kN)
        filename = os.path.join(OUTPUT_DIR, f"PDeltaData_{label}.txt")
        np.savetxt(filename, data, fmt='%12.6e', header='Axial_Strain Axial_Load(kN)', comments='')
        print(f"已儲存: {filename}")

# --- 繪圖 ---
plt.figure(figsize=(10,6))
for label, data in results.items():
    if len(data) > 0:
        plt.plot(data[:,0], data[:,1], label=label, linewidth=2)
        
plt.xlabel('Axial Strain', fontsize=12)
plt.ylabel('Axial Load (kN)', fontsize=12)
plt.title('Monotonic Axial Response of RC Columns', fontsize=14)
plt.legend()
plt.grid(True)

_curve_png = os.path.join(OUTPUT_DIR, "PDelta_Curve_Result.png")
plt.savefig(_curve_png, dpi=300, bbox_inches='tight')
print(f"已成功輸出圖片：{_curve_png}")

plt.show()



