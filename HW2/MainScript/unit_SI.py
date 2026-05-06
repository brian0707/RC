#============================================================================================================================================#
# Define Unit

# 這份程式碼定義了使用毫米（mm）、千牛頓（kN）和秒（sec）作為單位的系統。
# 包含了長度、面積、力和應力的轉換關係，並設置了一些常數，如圓周率（PI）和重力加速度（g）、極大值（Ubig）極小值（Usmall）等。
# 主單位： mm, kN, sec
# 可用單位：
    # 長度：mm, inch, ft, cm, m
    # 面積：inch2, mm2, cm2, m2
    # 慣性矩：inch4
    # 力：kip, N, MN, tf
    # 應力：MPa, ksi, psi, lbf, pcf
#============================================================================================================================================#
# !!!! 需要在主程式中匯入此程式包：
# !!!! import unit_SI
# !!!! unit_SI.define_units()
# !!!! 使用單位時需引用此檔名，例：unit_SI.mm, unit_SI.kN
#============================================================================================================================================#
import math

def define_units():
    global mm, kN, sec
    global LunitTXT, FunitTXT, TunitTXT
    global inch, ft, cm, m
    global inch2, mm2, cm2, m2
    global inch4
    global kip, N, MN, tf
    global MPa, GPa, ksi, psi, lbf, pcf
    global PI, g, Ubig, Usmall

    # Base units
    mm = 1.0
    kN = 1.0
    sec = 1.0

    LunitTXT = "mm"
    FunitTXT = "kN"
    TunitTXT = "sec"

    # Length
    inch = 25.4 * mm
    ft = 12.0 * inch
    cm = 10.0 * mm
    m = 1000.0 * mm

    # Area
    inch2 = inch * inch
    mm2 = mm * mm
    cm2 = cm * cm
    m2 = m * m

    # Moment of inertia
    inch4 = inch**4

    # Force
    kip = 4.448 * kN
    N = kN / 1000.0
    MN = 1000.0 * kN
    tf = 9.807 * kN

    # Stress
    MPa = N / mm2
    GPa = 1000*MPa
    ksi = kip / (inch**2)
    psi = ksi / 1000.0
    lbf = psi * inch * inch
    pcf = lbf / (ft**3)

    # Constants
    PI = 2.0 * math.asin(1.0)
    g = 32.2 * ft / (sec**2)  # gravitational acceleration
    Ubig = 1.0e10
    Usmall = 1.0 / Ubig

    print("-------------------------------------------------------------------------------------")
    print("Define Unit_SI done")
    print("-------------------------------------------------------------------------------------")