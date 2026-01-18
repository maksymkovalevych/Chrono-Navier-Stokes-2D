import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft2, ifft2, fftfreq

# --- 1. CONFIGURATION & PARAMETERS ---
N = 128           # Resolution (128 for speed, 256 for quality)
L = 2 * np.pi     # Domain size
Re = 10000.0      # Reynolds number
nu = 1.0 / Re     # Viscosity
dt = 0.001        # Time step
T_end = 2.0       # Duration
n_steps = int(T_end / dt)

# Grid setup
k = fftfreq(N) * N
kx, ky = np.meshgrid(k, k)
K2 = kx**2 + ky**2
K2[0, 0] = 1e-12  # Avoid division by zero
inv_K2 = 1.0 / K2
K = np.sqrt(K2)

# Dealiasing mask (2/3 rule)
mask = (np.abs(kx) < (N/3)) * (np.abs(ky) < (N/3))

# --- 2. CORE FUNCTIONS ---

def get_rhs(w_hat_local, gamma, c_light=100.0):
    """
    Computes the Right-Hand Side (RHS) of the equations.
    gamma: Chrono-coupling coefficient (0 = Standard, >0 = Chrono)
    """
    # 1. Streamfunction & Velocity
    psi_hat = -w_hat_local * inv_K2
    u = np.real(ifft2(1j * ky * psi_hat))
    v = np.real(ifft2(-1j * kx * psi_hat))
    w = np.real(ifft2(w_hat_local))
    
    # 2. Gradients of vorticity
    wx = np.real(ifft2(1j * kx * w_hat_local))
    wy = np.real(ifft2(1j * ky * w_hat_local))
    
    # 3. Nonlinear Term (Advection)
    adv = u * wx + v * wy
    adv_hat = fft2(adv) * mask
    
    # 4. Viscous Term (Diffusion) - Always acts in Lab Time t
    diff_hat = -nu * K2 * w_hat_local
    
    # 5. CHRONO-MODIFICATION
    if gamma > 0:
        Ko = gamma / (c_light**2)
        dtau_dt = np.exp(-Ko * (w**2))
        # Advection scales with local time, Diffusion does not
        rhs_hat = fft2(dtau_dt * (-adv)) * mask + diff_hat
    else:
        # Standard NS
        rhs_hat = -adv_hat + diff_hat
        
    return rhs_hat

def run_simulation(gamma_val):
    """Runs the simulation for a specific gamma."""
    print(f"Starting run for Gamma={gamma_val}...")
    
    # Initial Condition: Decaying Turbulence
    np.random.seed(42) # Fixed seed for fair comparison
    w0_hat = fft2(np.random.randn(N, N))
    w0_hat *= K * np.exp(-K2/20) * mask  # Band-limited noise
    
    w_hat = w0_hat.copy()
    
    # RK4 Time Stepping
    for step in range(n_steps):
        k1 = get_rhs(w_hat, gamma_val)
        k2 = get_rhs(w_hat + 0.5*dt*k1, gamma_val)
        k3 = get_rhs(w_hat + 0.5*dt*k2, gamma_val)
        k4 = get_rhs(w_hat + dt*k3, gamma_val)
        
        w_hat += (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4)
        
    return np.real(ifft2(w_hat))

# --- 3. RUNNING EXPERIMENTS ---
print("Running Standard NS...")
w_std = run_simulation(gamma_val=0.0)

print("Running Chrono-NS...")
w_chrono = run_simulation(gamma_val=10.0) # Higher gamma to see effect at low Re/Resolution

# --- 4. VISUALIZATION (THE DASHBOARD) ---
def plot_comparison(w_std, w_chrono):
    
    # Prep Statistics
    w_std_flat = w_std.flatten()
    w_chrono_flat = w_chrono.flatten()
    
    min_w = min(w_std_flat.min(), w_chrono_flat.min())
    max_w = max(w_std_flat.max(), w_chrono_flat.max())
    bins = np.linspace(min_w, max_w, 100)
    centers = 0.5 * (bins[1:] + bins[:-1])
    
    hist_std, _ = np.histogram(w_std_flat, bins=bins, density=True)
    hist_chrono, _ = np.histogram(w_chrono_flat, bins=bins, density=True)
    
    # Spectrum Helper
    def get_spectrum(w_field):
        w_hat = fft2(w_field)
        energy_dens = 0.5 * (np.abs(w_hat)**2) / (K2 + 1e-12)
        energy_dens[0,0] = 0
        k_bins = np.arange(0.5, N//3, 1)
        E_k = np.zeros(len(k_bins)-1)
        for i in range(len(k_bins)-1):
            m = (K >= k_bins[i]) & (K < k_bins[i+1])
            E_k[i] = np.sum(energy_dens[m])
        return 0.5*(k_bins[1:]+k_bins[:-1]), E_k

    k_ax, E_std = get_spectrum(w_std)
    _, E_chrono = get_spectrum(w_chrono)

    # Plotting
    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(2, 2)
    fig.suptitle(f"Chrono-NS vs Standard NS: Immersive Armor (Re={int(Re)})", fontsize=16, fontweight='bold')
    
    # Maps
    x = np.linspace(0, L, N)
    X, Y = np.meshgrid(x, x)
    
    ax1 = fig.add_subplot(gs[0, 0])
    im1 = ax1.imshow(w_std, cmap='RdBu_r', origin='lower', extent=[0, L, 0, L])
    ax1.contour(X, Y, w_std, levels=8, colors='k', alpha=0.3, linewidths=0.5)
    ax1.set_title('Standard NS (Chaos)', fontsize=12)
    plt.colorbar(im1, ax=ax1)
    
    ax2 = fig.add_subplot(gs[0, 1])
    im2 = ax2.imshow(w_chrono, cmap='RdBu_r', origin='lower', extent=[0, L, 0, L])
    ax2.contour(X, Y, w_chrono, levels=8, colors='k', alpha=0.3, linewidths=0.5)
    ax2.set_title('Chrono-NS (Armored Coherence)', fontsize=12)
    plt.colorbar(im2, ax=ax2)
    
    # PDF
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.semilogy(centers, hist_std, 'r-', alpha=0.6, lw=2, label='Standard')
    ax3.semilogy(centers, hist_chrono, 'b-', lw=2, label='Chrono-NS')
    ax3.set_title('PDF: Safety Cut-off')
    ax3.set_xlabel('Vorticity |ω|')
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    # Spectrum
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.loglog(k_ax, E_std, 'r-', alpha=0.6, label='Standard')
    ax4.loglog(k_ax, E_chrono, 'b-', label='Chrono-NS')
    
    # Ref line fitting
    if len(E_std) > 10:
        mid = len(E_std)//2
        ref = E_std[mid] * (k_ax/k_ax[mid])**(-3)
        ax4.loglog(k_ax, ref, 'k--', label='$k^{-3}$')
        
    ax4.set_title('Energy Spectrum: Adaptive Filter')
    ax4.set_xlabel('Wavenumber k')
    ax4.grid(True, alpha=0.3)
    ax4.legend()
    
    plt.tight_layout()
    plt.savefig('ChronoNS_Dashboard.png', dpi=150)
    plt.show()

# Run Viz
plot_comparison(w_std, w_chrono)
