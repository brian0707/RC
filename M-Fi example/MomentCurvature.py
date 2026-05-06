import csv
from pathlib import Path

import openseespy.opensees as ops


def _section_moment(ele_tag=1, sec_num=1):
    """2D fiber section: eleResponse returns [axial N, bending M]."""
    sf = ops.eleResponse(ele_tag, "section", sec_num, "force")
    return float(sf[1])


def MomentCurvature(secTag, axialLoad, maxK, numIncr=100):

    # Define two nodes at (0,0)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, 0.0)

    # Fix all degrees of freedom except axial and bending
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 0)

    # Define element
    #                             tag ndI ndJ  secTag
    ops.element("zeroLengthSection", 1, 1, 2, secTag)

    # Define constant axial load
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, axialLoad, 0.0, 0.0)

    # Define analysis parameters
    ops.integrator("LoadControl", 0.0)
    ops.system("SparseGeneral", "-piv")
    ops.test("NormUnbalance", 1e-9, 10)
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.algorithm("Newton")
    ops.analysis("Static")

    # Do one analysis for constant axial load
    ops.analyze(1)

    # Initial state on M–φ path (after axial load, before curvature cycles)
    kappa_hist = [ops.nodeDisp(2, 3)]
    moment_hist = [_section_moment()]

    # Define reference moment
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(2, 0.0, 0.0, 1.0)

    # Compute curvature increment
    dK = maxK / numIncr

    # Use displacement control at node 2 for section analysis
    ops.integrator("DisplacementControl", 2, 3, dK, 1, dK, dK)

    # Section analysis: one static step per increment so we can record M–φ
    for _ in range(numIncr):
        ops.analyze(1)
        kappa_hist.append(ops.nodeDisp(2, 3))
        moment_hist.append(_section_moment())

    return kappa_hist, moment_hist


def export_moment_curvature(
    kappa_hist,
    moment_hist,
    csv_path=None,
    plot_path=None,
):
    """Write curvature–moment series to CSV and optionally plot PNG."""
    if csv_path is None:
        csv_path = Path(__file__).resolve().parent / "moment_curvature.csv"
    else:
        csv_path = Path(csv_path)

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["curvature", "moment"])
        for k, m in zip(kappa_hist, moment_hist):
            w.writerow([k, m])

    if plot_path is None:
        plot_path = csv_path.with_suffix(".png")
    else:
        plot_path = Path(plot_path)

    plot_saved = None
    try:
        import matplotlib.pyplot as plt
        import matplotlib.ticker as mticker

        # 橫軸：φ·10^6，數值較大、較易比對差異（CSV 仍為原始曲率）
        x_plot = [k * 1e6 for k in kappa_hist]

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(x_plot, moment_hist, "-", linewidth=1.2)
        ax.set_xlabel(r"Curvature $\phi$ ($\times 10^{-6}$)")
        ax.set_ylabel("Moment M")
        ax.set_title("Moment–curvature")
        ax.grid(True, linestyle=":", alpha=0.7)
        ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=7, prune="both"))
        ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f"))
        fig.tight_layout()
        fig.savefig(plot_path, dpi=150)
        plt.close(fig)
        plot_saved = str(plot_path)
    except ImportError:
        pass

    return str(csv_path), plot_saved


ops.wipe()
print("Start MomentCurvature.py example")

# Define model builder
# --------------------
ops.model("basic", "-ndm", 2, "-ndf", 3)

# Define materials for nonlinear columns
# ------------------------------------------
# CONCRETE                  tag   f'c        ec0   f'cu        ecu
# Core concrete (confined)
ops.uniaxialMaterial("Concrete01", 1, -6.0, -0.004, -5.0, -0.014)

# Cover concrete (unconfined)
ops.uniaxialMaterial("Concrete01", 2, -5.0, -0.002, 0.0, -0.006)

# STEEL
# Reinforcing steel
fy = 60.0  # Yield stress
E = 30000.0  # Young's modulus

#                        tag  fy E0    b
ops.uniaxialMaterial("Steel01", 3, fy, E, 0.01)

# Define cross-section for nonlinear columns
# ------------------------------------------

# set some paramaters
colWidth = 15
colDepth = 24

cover = 1.5
As = 0.60  # area of no. 7 bars

# some variables derived from the parameters
y1 = colDepth / 2.0
z1 = colWidth / 2.0


ops.section("Fiber", 1)

# Create the concrete core fibers
ops.patch("rect", 1, 10, 1, cover - y1, cover - z1, y1 - cover, z1 - cover)

# Create the concrete cover fibers (top, bottom, left, right)
ops.patch("rect", 2, 10, 1, -y1, z1 - cover, y1, z1)
ops.patch("rect", 2, 10, 1, -y1, -z1, y1, cover - z1)
ops.patch("rect", 2, 2, 1, -y1, cover - z1, cover - y1, z1 - cover)
ops.patch("rect", 2, 2, 1, y1 - cover, cover - z1, y1, z1 - cover)

# Create the reinforcing fibers (left, middle, right)
ops.layer("straight", 3, 3, As, y1 - cover, z1 - cover, y1 - cover, cover - z1)
ops.layer("straight", 3, 2, As, 0.0, z1 - cover, 0.0, cover - z1)
ops.layer("straight", 3, 3, As, cover - y1, z1 - cover, cover - y1, cover - z1)

# Estimate yield curvature
# (Assuming no axial load and only top and bottom steel)
# d -- from cover to rebar
d = colDepth - cover
# steel yield strain
epsy = fy / E
Ky = epsy / (0.7 * d)

# Print estimate to standard output
print("Estimated yield curvature: ", Ky)

# Set axial load
P = -180.0

# Target ductility for analysis
mu = 15.0

# Number of analysis increments
numIncr = 100

# Call the section analysis procedure
kappa_hist, moment_hist = MomentCurvature(1, P, Ky * mu, numIncr)

csv_file, plot_file = export_moment_curvature(kappa_hist, moment_hist)
print(f"Moment–curvature data: {csv_file}")
if plot_file:
    print(f"Moment–curvature figure: {plot_file}")

results = open("results.out", "a+")

u = ops.nodeDisp(2, 3)
if abs(u - 0.00190476190476190541) < 1e-12:
    results.write("PASSED : MomentCurvature.py\n")
    print("Passed!")
else:
    results.write("FAILED : MomentCurvature.py\n")
    print("Failed!")

results.close()

print("==========================")
