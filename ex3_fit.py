import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import least_squares
import uproot, math

# Input/output
file = "fitInputs.root"
hdata_name = "hdata"
hbkg_name = "hbkg"
out_pdf = "ex3.pdf"

# Read 2D histograms
def read_hist(f, name):
    h = f[name]
    vals, xedges, yedges = h.to_numpy()
    return vals, xedges, yedges

f = uproot.open(file)
hdata, xedges, yedges = read_hist(f, hdata_name)
hbkg, _, _ = read_hist(f, hbkg_name)

nx, ny = hdata.shape
xcent = 0.5 * (xedges[:-1] + xedges[1:])
ycent = 0.5 * (yedges[:-1] + yedges[1:])
Xc, Yc = np.meshgrid(xcent, ycent, indexing="xy")
dx, dy = np.diff(xedges)[0], np.diff(yedges)[0]
dA = dx * dy

x, y = Xc.ravel(), Yc.ravel()
data, bkg = hdata.ravel(), hbkg.ravel()
sigma = np.sqrt(np.maximum(data, 1))
mask = np.isfinite(sigma)
x, y, data, bkg, sigma = x[mask], y[mask], data[mask], bkg[mask], sigma[mask]

# Model
def signal(p, x, y):
    A, mu1, mu2, s1, s2 = p
    return A * np.exp(-((x - mu1) ** 2) / s1 ** 2) * np.exp(-((y - mu2) ** 2) / s2 ** 2)

def model(p, x, y, bkg):
    A, mu1, mu2, s1, s2, B = p
    return signal([A, mu1, mu2, s1, s2], x, y) * dA + B * bkg

def residuals(p, x, y, bkg, data, sigma):
    if p[3] <= 0 or p[4] <= 0: return 1e6 * np.ones_like(data)
    return (data - model(p, x, y, bkg)) / sigma

# Initial guesses
i = np.argmax(data)
p0 = [np.max(data), x[i], y[i], (xedges[-1]-xedges[0])*0.1, (yedges[-1]-yedges[0])*0.1, 1]
lower = [0, xedges[0], yedges[0], 1e-6, 1e-6, -10]
upper = [np.inf, xedges[-1], yedges[-1], 10, 10, 10]

# Fit
res = least_squares(residuals, p0, args=(x, y, bkg, data, sigma), bounds=(lower, upper))
popt = res.x
J = res.jac
chi2 = np.sum(res.fun ** 2)
ndof = len(x) - len(popt)
cov = np.linalg.inv(J.T @ J) * (chi2 / ndof)
err = np.sqrt(np.diag(cov))

# Signal yield
A, mu1, mu2, s1, s2, B = popt
Nsig = A * math.pi * s1 * s2
g = np.array([math.pi * s1 * s2, math.pi * A * s2, math.pi * A * s1])
idx = [0, 3, 4]
cov_sub = cov[np.ix_(idx, idx)]
err_Nsig = math.sqrt(g @ cov_sub @ g)

# Prepare plots
model_all = model(popt, Xc.ravel(), Yc.ravel(), hbkg.ravel()).reshape(nx, ny)
resid = hdata - model_all
data_minus_bkg = hdata - B * hbkg

def lego(ax, X, Y, Z, title):
    ax.bar3d(X.ravel(), Y.ravel(), np.zeros_like(Z.ravel()),
             dx, dy, Z.ravel(), shade=True)
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("counts")

# 2x2 figure
fig = plt.figure(figsize=(14, 10))
ax1 = fig.add_subplot(221, projection="3d"); lego(ax1, Xc, Yc, hdata, "Data")
ax2 = fig.add_subplot(222, projection="3d"); lego(ax2, Xc, Yc, model_all, "Fit (Signal+Background)")
ax3 = fig.add_subplot(223, projection="3d"); lego(ax3, Xc, Yc, resid, "Residuals (Data-Fit)")
ax4 = fig.add_subplot(224, projection="3d"); lego(ax4, Xc, Yc, data_minus_bkg, "Data - Best-fit Background")
plt.tight_layout()

with PdfPages(out_pdf) as pdf:
    pdf.savefig(fig)
    fig2 = plt.figure(figsize=(8.5, 11))
    text = f"""
Fit results:
A = {popt[0]:.4g} ± {err[0]:.4g}
mu1 = {popt[1]:.4g} ± {err[1]:.4g}
mu2 = {popt[2]:.4g} ± {err[2]:.4g}
sigma1 = {popt[3]:.4g} ± {err[3]:.4g}
sigma2 = {popt[4]:.4g} ± {err[4]:.4g}
B = {popt[5]:.4g} ± {err[5]:.4g}

Chi2/ndof = {chi2:.2f}/{ndof} = {chi2/ndof:.3f}
Signal yield = {Nsig:.4g} ± {err_Nsig:.4g}
"""
    plt.text(0.05, 0.95, text, fontsize=12, va="top")
    pdf.savefig(fig2)
plt.close("all")

print("Fit done. Results saved in", out_pdf)

