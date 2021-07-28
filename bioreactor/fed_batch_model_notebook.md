---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.7.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

```python
import numpy as np
import matplotlib.pyplot as plt

def placeholder_cfd(our_max):
    '''
    This function is a stand in for a CFD + Paraview simulation that
    would be happening on the HPC.  Given a our_max value from the 
    time-stepping code, it should set the appropriate our_max_inhibited,
    run the OpenFOAM simulation, kick off the post-processing script, and
    finally return the reactor-average our.  Here, a percentage of the 
    input our_max is simply returned, somewhere around 93% since that's
    the average that was found in Hari's original results.
    '''

    #coeff = np.exp(-20.0*rho_furfural/(rho_glucose+rho_xylose))
    coeff = 1.0
    our_max_inh = coeff*our_max
    
    # A dummy operation, return somewhere between 90% and 96% of the our_max_inh
    our_cfd_avg = np.random.uniform(0.9, 0.96)*our_max_inh
    
    return our_cfd_avg

# ================================================================

def simulate_over_window(t_initial, t_final, initial_conditions, params):
    '''
    This function advances the forward-Euler time stepping starting
    from t_initial and ending at t_final using a timestep size dt.
    initial conditions is a dictionary with entries for up to 3 values.
    In the event a key:value is omitted (which happens during the 
    very first call) the initial value is calculated rather than copied
    from the preceding time window. params is a dictionary of constant or
    mutable param values which includes f_l, f_our, and lambda_our.
    '''
    
    # Set a constant timestep size
    dt = 2
    
    # Calculate the number of timesteps
    t_steps = int((t_final - t_initial)/dt)
    
    # Create a list of the timesteps (for visualization purposes)
    tt = np.linspace(t_initial, t_final, t_steps+1)
    
    # Initialize all arrays
    rho_l = np.zeros(t_steps + 1)
    rho_dcw = np.zeros(t_steps + 1)
    our_max = np.zeros(t_steps + 1)
    m_s = np.zeros(t_steps + 1)
    
    # Set initial conditions
    rho_dcw_init = initial_conditions.get('rho_dcw')
    if rho_dcw_init is None:
        raise ValueError('No initial value given for rho_dcw.')
    
    rho_l_init = initial_conditions.get('rho_l')
    if rho_l_init is None:
        print('No initial value given for rho_l, using calculated value.')
        rho_l_init = params['f_l']*rho_dcw_init
        
    our_max_init = initial_conditions.get('our_max')
    if our_max_init is None:
        print('No initial value given for our_max, using calculated value.')
        our_max_init = params['f_our']*params['lambda_our']*rho_dcw_init
        
    m_s_init = initial_conditions.get('m_s')
    if m_s_init is None:
        print('No initial value given for m_s, using calculated value.')
        m_s_init = (20 + 70)*0.3
    
    rho_l[0] = rho_l_init
    rho_dcw[0] = rho_dcw_init
    our_max[0] = our_max_init
    m_s[0] = m_s_init
    
    r_l_coeff = (0.2/1.75)*282.0*0.001
    r_s_coeff = (1.0/1.75)*180.0*0.001
    
    # Advance the equations in time
    for k in range(t_steps):
        r_s = -our_max[k]*r_s_coeff
        m_s[k+1] = m_s[k] + dt*r_s*params['vol']
        
        if m_s[k+1] > 0:
            r_l = our_max[k]*r_l_coeff
        else:
            r_l = 0
            m_s[k+1] = 0 
        
        rho_l[k+1] = rho_l[k] + dt*r_l
        rho_dcw[k+1] = rho_l[k+1]/params['f_l']
        our_max[k+1] = params['f_our']*params['lambda_our']*rho_dcw[k+1]
    
    print_ranges = True
    
    if print_ranges == True:
        print('\nCalculated Ranges for this window')
        print('| rho_l:   [%6.2f, ..., %6.2f]' % (rho_l[0], rho_l[-1]))
        print('| rho_dcw: [%6.2f, ..., %6.2f]' % (rho_dcw[0], rho_dcw[-1]))
        print('| our_max: [%6.2f, ..., %6.2f]' % (our_max[0], our_max[-1]))
        print('| m_s: [%6.2f, ..., %6.2f]' % (m_s[0], m_s[-1]))
        print('| tt:      [%6.2f, ..., %6.2f]\n' % (tt[0], tt[-1]))

    return rho_l, rho_dcw, our_max, m_s, tt
    
# ================================================================

# OURcfd = np.array([4.54, 8.21]) # mol/m^3/h, per Hari's runs
# fOUR = 0.93 # OUR_actual/OUR_max; average of [4.54/4.8, 8.21/8.96]
# Find out where the values 4.8 (ic), 8.96 (?) came from

# Simulation constants and parameters
params = {}
params['lambda_our'] = 0.48
params['f_l'] = 0.5
params['vol'] = 0.3 # reactor volume in liters

# Set the number of CFD evaluations to be performed throughout simulation
num_cfd_evals = 4
t_initial = 12.0
t_final = 96.0
t_split = np.linspace(t_initial, t_final, num_cfd_evals+1)

# ================================================================

# Advance the simulation over the specified time window from the starting point
# Begin with the assumption that f_our = 1.0
# params['f_our'] = 0.93 # Constant based on Hari's CFD
params['f_our'] = 1.0

# Set all initial conditions (initially, only rho_dcw is known)
ic = {}
ic['rho_dcw'] = 10.0

# Initialize values outside the loop
rho_l, rho_dcw, our_max, m_s, tt = simulate_over_window(t_split[0], t_split[0], ic, params)

# Initialize lists to hold results (for later visualization)
rho_l_total = []
rho_dcw_total = []
our_max_total = []
m_s_total = []
tt_total = []

for k in range(num_cfd_evals):
    print('========Simulating from %.1f to %.1f========' %(t_split[k], t_split[k+1]))
    
    # Placeholder for CFD simulation
    our_cfd_avg = placeholder_cfd(our_max[-1])

    # Update the value of f_our based on the latest CFD and our_max calculations
    # params['f_our'] = 0.93 # Constant based on Hari's CFD
    params['f_our'] = our_cfd_avg/our_max[-1]
    print('our_cfd_avg/our_max = f_our = %.2f/%.2f = %.3f' % (our_cfd_avg, our_max[-1], params['f_our']))

    # Update all initial conditions used for this new time window
    ic['rho_l'] = rho_l[-1]
    ic['rho_dcw'] = rho_dcw[-1]
    ic['our_max'] = our_max[-1]
    ic['m_s'] = m_s[-1]

    # Advance the simulation
    rho_l, rho_dcw, our_max, m_s, tt = simulate_over_window(t_split[k], t_split[k+1], ic, params)
    
    # Store the results in a global array (for later visualization)
    rho_l_total.append(rho_l)
    rho_dcw_total.append(rho_dcw)
    our_max_total.append(our_max)
    m_s_total.append(m_s)
    tt_total.append(tt)

```

```python
# Format and plot the data
rho_l_total = np.hstack(rho_l_total)
rho_dcw_total = np.hstack(rho_dcw_total)
our_max_total = np.hstack(our_max_total)
m_s_total = np.hstack(m_s_total)
tt_total = np.hstack(tt_total)

plt.figure(figsize=(7, 5), dpi=100)
plt.plot(tt_total, rho_dcw_total, '.-', label='DCW')
plt.plot(tt_total, rho_l_total, '.-', label='Lipid')
plt.plot(tt_total, our_max_total, '.-', label='OUR')
plt.plot(tt_total, m_s_total, '.-', label='Sugar')
for k, x in enumerate(t_split):
    if k == 0:
        plt.axvline(x, color='k', linestyle=':', zorder=0, label='CFD Eval')
    else:
        plt.axvline(x, color='k', linestyle=':', zorder=0)

# approximate experimental data from Fei et al. (2016), uncomment to include
# texp = np.array([0, 24, 48, 72, 96])
# DCWexp = np.array([0, 20, 33, 51, 53])
# lipexp = np.array([0, 7, 14, 28, 30])

# plt.plot(texp, DCWexp, 'o--C0', mfc="none", label='DCW (exp)')
# plt.plot(texp, lipexp, 'o--C1', mfc="none", label='Lipid (exp)')

plt.xlabel('Time (hr)')
plt.ylabel("DCW, Lipid [g/L]; OUR [mol/($m^3$.h)]; Sugar [g]")
plt.legend()
plt.show()

```

```python

```
