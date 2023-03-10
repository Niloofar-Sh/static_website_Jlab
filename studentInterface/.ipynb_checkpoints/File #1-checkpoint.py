# Size of variable arrays:
sizeAlgebraic = 2
sizeStates = 4
sizeConstants = 4
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "t in component main (second)"
    legend_states[0] = "q_0 in component main (concentrationUnit)"
    legend_states[1] = "q_1 in component main (concentrationUnit)"
    legend_states[2] = "q_2 in component main (concentrationUnit)"
    legend_states[3] = "q_3 in component main (concentrationUnit)"
    legend_algebraic[0] = "v_0 in component main (fluxUnit)"
    legend_algebraic[1] = "v_1 in component main (fluxUnit)"
    legend_constants[0] = "kf_0 in component main (fluxUnit)"
    legend_constants[1] = "kr_0 in component main (fluxUnit)"
    legend_constants[2] = "kf_1 in component main (fluxUnit)"
    legend_constants[3] = "kr_1 in component main (fluxUnit)"
    legend_rates[0] = "d/dt q_0 in component main (concentrationUnit)"
    legend_rates[1] = "d/dt q_1 in component main (concentrationUnit)"
    legend_rates[2] = "d/dt q_2 in component main (concentrationUnit)"
    legend_rates[3] = "d/dt q_3 in component main (concentrationUnit)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 0.105
    states[1] = 0.047
    states[2] = 0.05
    states[3] = 0.003
    constants[0] = 0.647352377784473
    constants[1] = 10.789206255809654
    constants[2] = 1.0537671758687692
    constants[3] = 0.9398323867797539
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[1] = 0.178699*constants[2]*(power(states[0], 2.00000))-constants[3]*(power(states[1], 2.00000))
    rates[0] = -2.00000*algebraic[1]
    rates[1] = (2.00000)*algebraic[1]
    algebraic[0] = constants[0]*(power(states[2], 1.00000))-constants[1]*(power(states[3], 1.00000))
    rates[2] = -1.00000*algebraic[0]
    rates[3] = (1.00000)*algebraic[0]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[1] = 0.178699*constants[2]*(power(states[0], 2.00000))-constants[3]*(power(states[1], 2.00000))
    algebraic[0] = constants[0]*(power(states[2], 1.00000))-constants[1]*(power(states[3], 1.00000))
    return algebraic

def solve_model():
    """Solve model with ODE solver"""
    from scipy.integrate import ode
    # Initialise constants and state variables
    (init_states, constants) = initConsts()

    # Set timespan to solve over
    voi = linspace(0, 10, 500)

    # Construct ODE object to solve
    r = ode(computeRates)
    r.set_integrator('vode', method='bdf', atol=1e-06, rtol=1e-06, max_step=1)
    r.set_initial_value(init_states, voi[0])
    r.set_f_params(constants)

    # Solve model
    states = array([[0.0] * len(voi)] * sizeStates)
    states[:,0] = init_states
    for (i,t) in enumerate(voi[1:]):
        if r.successful():
            r.integrate(t)
            states[:,i+1] = r.y
        else:
            break

    # Compute algebraic variables
    algebraic = computeAlgebraic(constants, states, voi)
    return (voi, states, algebraic)

def plot_model(voi, states, algebraic):
    """Plot variables against variable of integration"""
    import pylab
    (legend_states, legend_algebraic, legend_voi, legend_constants) = createLegends()
    pylab.figure(1)
    pylab.plot(voi,vstack((states,algebraic)).T)
    pylab.xlabel(legend_voi)
    pylab.legend(legend_states + legend_algebraic, loc='best')
    pylab.show()

if __name__ == "__main__":
    (voi, states, algebraic) = solve_model()
    plot_model(voi, states, algebraic)
