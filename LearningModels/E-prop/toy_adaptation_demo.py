"""
Toy Adaptation Demo Script
===========================
A comparison of spike-frequency adaptation (SFA) implementations:

1. E-PROP / ALIF MODEL (Your current implementation)
   ------------------------------------------------
   Membrane voltage equation:
       tau * dV/dt = (E_L - V) - R * g_ad * (V - E_k) + I_input

   Adaptation conductance equation:
       tau_ad * dg_ad/dt = -g_ad

   On spike:
       V -> V_reset
       g_ad -> g_ad + g_inc

2. NEURON / AHP MODEL (After-Hyperpolarization)
   ---------------------------------------------
   NEURON typically implements adaptation via calcium-activated K+ currents.
   
   Membrane voltage equation:
       C * dV/dt = g_L*(E_L - V) - g_AHP*[Ca]*(V - E_K) + I_input

   Calcium dynamics (increments on spike, decays exponentially):
       d[Ca]/dt = -[Ca] / tau_Ca

   On spike:
       V -> V_reset
       [Ca] -> [Ca] + Ca_inc

   The key difference: NEURON explicitly models calcium as the "memory" of
   recent spiking, and the AHP current is proportional to calcium concentration.

   Common NEURON mechanisms for adaptation:
   - IKahp.mod: Calcium-dependent K+ current (SK channels)
   - IKm.mod: M-type slow K+ current (voltage-dependent)
   - CaHVA/CaLVA: Calcium channels that feed into AHP

   Both approaches are mathematically equivalent when:
       g_ad ≈ g_AHP * [Ca]
       tau_ad ≈ tau_Ca
       g_inc ≈ g_AHP * Ca_inc
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for saving figures
import matplotlib.pyplot as plt

# =============================================================================
# NEURON PARAMETERS (matching your Lif_Neuron.py defaults where applicable)
# =============================================================================

# Membrane parameters
C = 0.1          # Membrane capacitance [nF]
g_L = 1/200      # Leak conductance [uS]
R = 1/g_L        # Membrane resistance [MOhm] = 200 MOhm
tau = C * R      # Membrane time constant [ms] = 20 ms
E_L = -65        # Resting potential [mV]
V_thresh = -47   # Spike threshold [mV]
V_reset = -54    # Reset potential [mV]

# Adaptation parameters (E-prop style)
E_k = -80        # Adaptation reversal potential [mV] (hyperpolarizing)
tau_ad = 50      # Adaptation time constant [ms] - controls how fast adaptation decays
g_inc = 0.005    # Adaptation increment on spike [uS] - controls adaptation strength

# AHP parameters (NEURON style) - for comparison
# These are typical values from NEURON's IKahp mechanism
g_AHP = 0.01     # AHP conductance [uS] - max conductance of AHP channel
tau_Ca = 50      # Calcium decay time constant [ms] - same as tau_ad for comparison
Ca_inc = 0.5     # Calcium increment per spike [arbitrary units]
E_K_ahp = -80    # K+ reversal for AHP [mV] - same as E_k

# =============================================================================
# SIMULATION PARAMETERS
# =============================================================================

dt = 0.05        # Time step [ms]
T_total = 500    # Total simulation time [ms]
n_steps = int(T_total / dt)
time = np.arange(0, T_total, dt)

# Input current - constant injection to drive spiking
I_tonic = 0.15   # [nA] - enough to drive repeated spiking

# =============================================================================
# SIMULATION FUNCTIONS
# =============================================================================

def simulate_lif_with_adaptation(tau_ad_val, g_inc_val, label=""):
    """
    Simulate an LIF neuron with spike-frequency adaptation.
    
    Parameters:
    -----------
    tau_ad_val : float
        Adaptation time constant [ms]
    g_inc_val : float
        Adaptation increment per spike [uS]
    label : str
        Label for plotting
        
    Returns:
    --------
    V_trace : array
        Voltage trace over time
    g_ad_trace : array
        Adaptation conductance trace over time
    spike_times : list
        Times of spikes [ms]
    """
    
    # Initialize state variables
    V = E_L           # Start at resting potential
    g_ad = 0.0        # Start with no adaptation
    
    # Storage
    V_trace = np.zeros(n_steps)
    g_ad_trace = np.zeros(n_steps)
    I_adapt_trace = np.zeros(n_steps)
    spike_times = []
    
    for i in range(n_steps):
        # Store current state
        V_trace[i] = V
        g_ad_trace[i] = g_ad
        
        # Compute adaptation current (this opposes depolarization)
        # I_adapt = g_ad * (V - E_k)
        # Since E_k = -80 < V, and g_ad >= 0, this term is positive
        # In the ODE: dV/dt = ... - R * g_ad * (V - E_k) ...
        # So this hyperpolarizes the neuron
        I_adapt = g_ad * (V - E_k)
        I_adapt_trace[i] = I_adapt
        
        # Membrane voltage ODE (Euler method)
        # dV/dt = [(E_L - V) - R * g_ad * (V - E_k) + R * I_tonic] / tau
        dV_dt = ((E_L - V) - R * g_ad * (V - E_k) + R * I_tonic) / tau
        
        # Adaptation conductance ODE
        # dg_ad/dt = -g_ad / tau_ad
        dg_ad_dt = -g_ad / tau_ad_val
        
        # Euler step
        V_new = V + dt * dV_dt
        g_ad_new = g_ad + dt * dg_ad_dt
        
        # Check for spike (threshold crossing)
        if V_new >= V_thresh:
            spike_times.append(time[i])
            # Reset voltage
            V_new = V_reset
            # Increment adaptation conductance
            g_ad_new = g_ad_new + g_inc_val
        
        # Update state
        V = V_new
        g_ad = g_ad_new
    
    return V_trace, g_ad_trace, I_adapt_trace, spike_times


def compute_isi(spike_times):
    """Compute inter-spike intervals."""
    if len(spike_times) < 2:
        return []
    return np.diff(spike_times)


def simulate_neuron_ahp(g_AHP_val, tau_Ca_val, Ca_inc_val, label=""):
    """
    Simulate an LIF neuron with NEURON-style AHP (After-Hyperpolarization).
    
    This mimics NEURON's IKahp.mod mechanism where:
    - Calcium accumulates with each spike
    - AHP current is proportional to calcium: I_AHP = g_AHP * [Ca] * (V - E_K)
    - Calcium decays exponentially: d[Ca]/dt = -[Ca] / tau_Ca
    
    Parameters:
    -----------
    g_AHP_val : float
        Maximum AHP conductance [uS]
    tau_Ca_val : float
        Calcium decay time constant [ms]
    Ca_inc_val : float
        Calcium increment per spike [a.u.]
    label : str
        Label for plotting
        
    Returns:
    --------
    V_trace : array
        Voltage trace over time
    Ca_trace : array
        Calcium concentration trace over time
    I_ahp_trace : array
        AHP current trace over time
    spike_times : list
        Times of spikes [ms]
    """
    
    # Initialize state variables
    V = E_L           # Start at resting potential
    Ca = 0.0          # Start with no calcium
    
    # Storage
    V_trace = np.zeros(n_steps)
    Ca_trace = np.zeros(n_steps)
    I_ahp_trace = np.zeros(n_steps)
    spike_times = []
    
    for i in range(n_steps):
        # Store current state
        V_trace[i] = V
        Ca_trace[i] = Ca
        
        # Compute AHP current (NEURON style)
        # I_AHP = g_AHP * [Ca] * (V - E_K)
        # This is the key difference: current scales with BOTH conductance AND calcium
        I_ahp = g_AHP_val * Ca * (V - E_K_ahp)
        I_ahp_trace[i] = I_ahp
        
        # Membrane voltage ODE (using conductance form like NEURON)
        # C * dV/dt = g_L*(E_L - V) - g_AHP*[Ca]*(V - E_K) + I_input
        # Rearranging: dV/dt = [g_L*(E_L - V) - g_AHP*Ca*(V - E_K) + I_input] / C
        dV_dt = (g_L * (E_L - V) - g_AHP_val * Ca * (V - E_K_ahp) + I_tonic) / C
        
        # Calcium dynamics ODE
        # d[Ca]/dt = -[Ca] / tau_Ca
        dCa_dt = -Ca / tau_Ca_val
        
        # Euler step
        V_new = V + dt * dV_dt
        Ca_new = Ca + dt * dCa_dt
        
        # Check for spike (threshold crossing)
        if V_new >= V_thresh:
            spike_times.append(time[i])
            # Reset voltage
            V_new = V_reset
            # Increment calcium (spike-triggered calcium influx)
            Ca_new = Ca_new + Ca_inc_val
        
        # Update state
        V = V_new
        Ca = max(0, Ca_new)  # Calcium can't go negative
    
    return V_trace, Ca_trace, I_ahp_trace, spike_times


# =============================================================================
# RUN SIMULATIONS
# =============================================================================

print("=" * 70)
print("Toy Adaptation Demo - E-prop vs NEURON AHP Comparison")
print("=" * 70)

# --- E-PROP STYLE SIMULATIONS ---
print("\n--- E-PROP / ALIF Style (Your Implementation) ---")

# Simulation 1: No adaptation (baseline)
V1, g_ad1, I_adapt1, spikes1 = simulate_lif_with_adaptation(
    tau_ad_val=50, g_inc_val=0.0, label="No Adaptation"
)

# Simulation 2: Weak adaptation
V2, g_ad2, I_adapt2, spikes2 = simulate_lif_with_adaptation(
    tau_ad_val=50, g_inc_val=0.003, label="Weak Adaptation"
)

# Simulation 3: Strong adaptation (your default)
V3, g_ad3, I_adapt3, spikes3 = simulate_lif_with_adaptation(
    tau_ad_val=50, g_inc_val=0.005, label="Strong Adaptation"
)

# Simulation 4: Very strong + fast decay
V4, g_ad4, I_adapt4, spikes4 = simulate_lif_with_adaptation(
    tau_ad_val=20, g_inc_val=0.008, label="Strong + Fast Decay"
)

# --- NEURON AHP STYLE SIMULATIONS ---
print("\n--- NEURON / AHP Style (Calcium-based) ---")

# Simulation 5: NEURON AHP - No adaptation
V5, Ca5, I_ahp5, spikes5 = simulate_neuron_ahp(
    g_AHP_val=0.0, tau_Ca_val=50, Ca_inc_val=0.5, label="No AHP"
)

# Simulation 6: NEURON AHP - Weak
V6, Ca6, I_ahp6, spikes6 = simulate_neuron_ahp(
    g_AHP_val=0.005, tau_Ca_val=50, Ca_inc_val=0.5, label="Weak AHP"
)

# Simulation 7: NEURON AHP - Strong (matched to E-prop strong)
V7, Ca7, I_ahp7, spikes7 = simulate_neuron_ahp(
    g_AHP_val=0.01, tau_Ca_val=50, Ca_inc_val=0.5, label="Strong AHP"
)

# Simulation 8: NEURON AHP - Fast calcium decay
V8, Ca8, I_ahp8, spikes8 = simulate_neuron_ahp(
    g_AHP_val=0.015, tau_Ca_val=20, Ca_inc_val=0.5, label="Fast Ca Decay"
)

# =============================================================================
# PRINT RESULTS
# =============================================================================

print("\n" + "=" * 70)
print("SPIKE COUNT COMPARISON")
print("=" * 70)

print("\nE-PROP / ALIF Style:")
print(f"  No Adaptation (g_inc=0):        {len(spikes1)} spikes")
print(f"  Weak Adaptation (g_inc=0.003):  {len(spikes2)} spikes")
print(f"  Strong Adaptation (g_inc=0.005): {len(spikes3)} spikes")
print(f"  Strong+Fast (g_inc=0.008, τ=20): {len(spikes4)} spikes")

print("\nNEURON / AHP Style:")
print(f"  No AHP (g_AHP=0):               {len(spikes5)} spikes")
print(f"  Weak AHP (g_AHP=0.005):         {len(spikes6)} spikes")
print(f"  Strong AHP (g_AHP=0.01):        {len(spikes7)} spikes")
print(f"  Fast Ca Decay (τ_Ca=20):        {len(spikes8)} spikes")

print("\n--- Inter-Spike Interval (ISI) Analysis ---")
print("\nE-PROP Style:")
for label, spikes in [("No Adapt", spikes1), ("Weak", spikes2), 
                       ("Strong", spikes3), ("Fast Decay", spikes4)]:
    isi = compute_isi(spikes)
    if len(isi) > 0:
        print(f"  {label:12s}: First ISI={isi[0]:.2f}ms, Last ISI={isi[-1]:.2f}ms, "
              f"Ratio={isi[-1]/isi[0]:.2f}x")
    else:
        print(f"  {label:12s}: Not enough spikes for ISI analysis")

print("\nNEURON AHP Style:")
for label, spikes in [("No AHP", spikes5), ("Weak AHP", spikes6), 
                       ("Strong AHP", spikes7), ("Fast Ca", spikes8)]:
    isi = compute_isi(spikes)
    if len(isi) > 0:
        print(f"  {label:12s}: First ISI={isi[0]:.2f}ms, Last ISI={isi[-1]:.2f}ms, "
              f"Ratio={isi[-1]/isi[0]:.2f}x")
    else:
        print(f"  {label:12s}: Not enough spikes for ISI analysis")

# =============================================================================
# PLOTTING - COMPARISON FIGURE
# =============================================================================

fig, axes = plt.subplots(3, 2, figsize=(14, 10))

# Left column: E-PROP style
# Right column: NEURON AHP style

# Row 1: Voltage traces
ax1 = axes[0, 0]
ax1.plot(time, V1, 'b-', alpha=0.7, label='No Adaptation')
ax1.plot(time, V3, 'r-', alpha=0.7, label='Strong (g_inc=0.005)')
ax1.axhline(V_thresh, color='k', linestyle='--', alpha=0.3)
ax1.set_ylabel('Voltage [mV]')
ax1.set_ylim(-75, -40)
ax1.legend(loc='upper right', fontsize=8)
ax1.set_title('E-PROP / ALIF Style\n(Your Implementation)')

ax2 = axes[0, 1]
ax2.plot(time, V5, 'b-', alpha=0.7, label='No AHP')
ax2.plot(time, V7, 'r-', alpha=0.7, label='Strong AHP (g_AHP=0.01)')
ax2.axhline(V_thresh, color='k', linestyle='--', alpha=0.3)
ax2.set_ylabel('Voltage [mV]')
ax2.set_ylim(-75, -40)
ax2.legend(loc='upper right', fontsize=8)
ax2.set_title('NEURON / AHP Style\n(Calcium-based)')

# Row 2: Adaptation variable
ax3 = axes[1, 0]
ax3.plot(time, g_ad1 * 1000, 'b-', alpha=0.7, label='No Adapt')
ax3.plot(time, g_ad2 * 1000, 'g-', alpha=0.7, label='Weak')
ax3.plot(time, g_ad3 * 1000, 'r-', alpha=0.7, label='Strong')
ax3.set_ylabel('g_ad [nS]')
ax3.legend(loc='upper right', fontsize=8)
ax3.set_xlabel('Time [ms]')

ax4 = axes[1, 1]
ax4.plot(time, Ca5, 'b-', alpha=0.7, label='No AHP')
ax4.plot(time, Ca6, 'g-', alpha=0.7, label='Weak AHP')
ax4.plot(time, Ca7, 'r-', alpha=0.7, label='Strong AHP')
ax4.set_ylabel('[Ca] (a.u.)')
ax4.legend(loc='upper right', fontsize=8)
ax4.set_xlabel('Time [ms]')

# Row 3: Spike rasters
ax5 = axes[2, 0]
for idx, (spikes, color, label) in enumerate([
    (spikes1, 'b', 'No Adapt'),
    (spikes2, 'g', 'Weak'),
    (spikes3, 'r', 'Strong'),
    (spikes4, 'm', 'Fast Decay')
]):
    ax5.eventplot(spikes, lineoffsets=idx, colors=color, label=label)
ax5.set_yticks([0, 1, 2, 3])
ax5.set_yticklabels(['None', 'Weak', 'Strong', 'Fast'])
ax5.set_xlabel('Time [ms]')
ax5.set_ylabel('E-PROP')

ax6 = axes[2, 1]
for idx, (spikes, color, label) in enumerate([
    (spikes5, 'b', 'No AHP'),
    (spikes6, 'g', 'Weak AHP'),
    (spikes7, 'r', 'Strong AHP'),
    (spikes8, 'm', 'Fast Ca')
]):
    ax6.eventplot(spikes, lineoffsets=idx, colors=color, label=label)
ax6.set_yticks([0, 1, 2, 3])
ax6.set_yticklabels(['None', 'Weak', 'Strong', 'Fast'])
ax6.set_xlabel('Time [ms]')
ax6.set_ylabel('NEURON AHP')

plt.suptitle('Spike-Frequency Adaptation: E-PROP vs NEURON Comparison', fontsize=14)
plt.tight_layout()
plt.savefig('adaptation_eprop_vs_neuron.png', dpi=150)
print("\nFigure saved: adaptation_eprop_vs_neuron.png")
# plt.show()  # Commented for non-interactive mode

# =============================================================================
# ADDITIONAL ANALYSIS: ISI over time for both models
# =============================================================================

fig2, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(12, 5))

# E-PROP style
for spikes, color, label in [
    (spikes1, 'b', 'No Adaptation'),
    (spikes2, 'g', 'Weak (g_inc=0.003)'),
    (spikes3, 'r', 'Strong (g_inc=0.005)'),
    (spikes4, 'm', 'Fast Decay (τ=20ms)')
]:
    isi = compute_isi(spikes)
    if len(isi) > 1:
        ax_left.plot(range(1, len(isi)+1), isi, 'o-', color=color, label=label, alpha=0.7)

ax_left.set_xlabel('Spike Number')
ax_left.set_ylabel('Inter-Spike Interval [ms]')
ax_left.set_title('E-PROP / ALIF Style')
ax_left.legend(fontsize=8)
ax_left.grid(True, alpha=0.3)

# NEURON AHP style
for spikes, color, label in [
    (spikes5, 'b', 'No AHP'),
    (spikes6, 'g', 'Weak AHP'),
    (spikes7, 'r', 'Strong AHP'),
    (spikes8, 'm', 'Fast Ca Decay')
]:
    isi = compute_isi(spikes)
    if len(isi) > 1:
        ax_right.plot(range(1, len(isi)+1), isi, 'o-', color=color, label=label, alpha=0.7)

ax_right.set_xlabel('Spike Number')
ax_right.set_ylabel('Inter-Spike Interval [ms]')
ax_right.set_title('NEURON / AHP Style')
ax_right.legend(fontsize=8)
ax_right.grid(True, alpha=0.3)

plt.suptitle('ISI vs Spike Number - Both Models Show Similar Adaptation Dynamics', fontsize=12)
plt.tight_layout()
plt.savefig('isi_comparison_both_models.png', dpi=150)
print("Figure saved: isi_comparison_both_models.png")
# plt.show()  # Commented for non-interactive mode

# =============================================================================
# DETAILED COMPARISON PRINTOUT
# =============================================================================

print("\n" + "=" * 70)
print("E-PROP vs NEURON AHP: DETAILED COMPARISON")
print("=" * 70)

print("""
┌─────────────────────────────────────────────────────────────────────┐
│                    E-PROP / ALIF MODEL                              │
│                 (Your current implementation)                        │
├─────────────────────────────────────────────────────────────────────┤
│ State Variable:  g_ad (adaptation conductance)                      │
│                                                                     │
│ Equations:                                                          │
│   dV/dt = [(E_L - V) - R·g_ad·(V - E_k) + R·I_input] / τ           │
│   dg_ad/dt = -g_ad / τ_ad                                           │
│                                                                     │
│ On Spike:                                                           │
│   g_ad → g_ad + g_inc                                               │
│                                                                     │
│ Parameters:                                                         │
│   - E_k:    Reversal potential (typically -80 mV, hyperpolarizing)  │
│   - τ_ad:   Adaptation time constant (your default: 5 ms - FAST!)   │
│   - g_inc:  Conductance increment per spike                         │
│                                                                     │
│ Pros: Simple, differentiable, good for gradient-based learning     │
│ Cons: Abstract, doesn't map directly to ion channels                │
└─────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────┐
│                    NEURON / AHP MODEL                               │
│               (Biophysically-inspired approach)                      │
├─────────────────────────────────────────────────────────────────────┤
│ State Variable:  [Ca] (intracellular calcium concentration)         │
│                                                                     │
│ Equations:                                                          │
│   C·dV/dt = g_L·(E_L - V) - g_AHP·[Ca]·(V - E_K) + I_input         │
│   d[Ca]/dt = -[Ca] / τ_Ca                                           │
│                                                                     │
│ On Spike:                                                           │
│   [Ca] → [Ca] + Ca_inc  (spike-triggered calcium influx)            │
│                                                                     │
│ Common NEURON Mechanisms:                                           │
│   - IKahp.mod:  SK channels (calcium-activated K+)                  │
│   - IKCa.mod:   BK channels (big conductance Ca-activated K+)       │
│   - CaDyn.mod:  Calcium dynamics (influx, buffering, extrusion)     │
│   - mAHP:       Medium AHP (50-200 ms, SK channels)                 │
│   - sAHP:       Slow AHP (1-5 s, different mechanism)               │
│                                                                     │
│ Pros: Biophysically realistic, maps to real channels                │
│ Cons: More complex, harder to differentiate for learning            │
└─────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────┐
│                     MATHEMATICAL EQUIVALENCE                        │
├─────────────────────────────────────────────────────────────────────┤
│ The two models are mathematically equivalent when:                  │
│                                                                     │
│   E-PROP g_ad  ≈  NEURON g_AHP · [Ca]                               │
│   E-PROP τ_ad  ≈  NEURON τ_Ca                                       │
│   E-PROP g_inc ≈  NEURON g_AHP · Ca_inc                             │
│                                                                     │
│ The E-PROP model essentially "pre-multiplies" the conductance       │
│ and calcium into a single variable (g_ad), making it simpler        │
│ for gradient computation but losing the explicit calcium dynamics.  │
└─────────────────────────────────────────────────────────────────────┘
""")

print("\n--- Parameter Mapping ---")
print(f"Your E-PROP: τ_ad = {tau_ad} ms, g_inc = {g_inc}")
print(f"Equivalent NEURON: τ_Ca = {tau_Ca} ms, g_AHP × Ca_inc = {g_AHP * Ca_inc:.4f}")

print("""
┌─────────────────────────────────────────────────────────────────────┐
│                  TYPES OF AHP IN NEURON                             │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│ 1. fAHP (Fast AHP): τ ~ 1-5 ms                                      │
│    - Caused by BK channels (large conductance Ca-activated K+)      │
│    - Immediate repolarization after spike                           │
│    - NOT typically used for spike-frequency adaptation              │
│                                                                     │
│ 2. mAHP (Medium AHP): τ ~ 50-200 ms                                 │
│    - Caused by SK channels (small conductance Ca-activated K+)      │
│    - THIS is what produces spike-frequency adaptation               │
│    - Your τ_ad=5ms is much faster than typical mAHP!                │
│                                                                     │
│ 3. sAHP (Slow AHP): τ ~ 1-5 seconds                                 │
│    - Mechanism debated (possibly KCNQ channels)                     │
│    - Produces very slow adaptation over seconds                     │
│                                                                     │
│ RECOMMENDATION: Set τ_ad = 50-200 ms to match realistic mAHP        │
└─────────────────────────────────────────────────────────────────────┘
""")
