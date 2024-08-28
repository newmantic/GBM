# GBM

Geometric Brownian Motion (GBM) is a widely used model in finance for modeling the price dynamics of assets, such as stocks, over time. The GBM model assumes that the logarithm of the asset price follows a Brownian motion (or Wiener process) with drift, making it a continuous-time stochastic process. This model is the foundation of the Black-Scholes option pricing model.


Asset Price (S_t): The asset price at time t is denoted by S_t. The GBM model assumes that S_t evolves over time according to a stochastic differential equation (SDE).

Stochastic Differential Equation (SDE): The dynamics of the asset price S_t under GBM are governed by the following SDE:
dS_t = mu * S_t * dt + sigma * S_t * dW_t
where:
mu is the drift coefficient, representing the expected return of the asset.
sigma is the volatility coefficient, representing the standard deviation of the asset's returns.
W_t is a Wiener process (or Brownian motion), which introduces randomness into the model.
dW_t is the increment of the Wiener process, which is normally distributed with mean 0 and variance dt.

Interpretation:
The term mu * S_t * dt represents the deterministic part of the asset's price movement, contributing to the expected growth of the asset.
The term sigma * S_t * dW_t represents the stochastic part of the asset's price movement, contributing to the randomness in the price due to market fluctuations.

Solution to the SDE: The solution to the SDE can be written as:
S_t = S_0 * exp((mu - 0.5 * sigma^2) * t + sigma * W_t)
where:
S_0 is the initial asset price at time t = 0.
exp represents the exponential function.
W_t is the Wiener process at time t.

Log-Normal Distribution: Under the GBM model, the logarithm of the asset price log(S_t) is normally distributed. This implies that the asset price S_t itself follows a log-normal distribution:
log(S_t) ~ N(log(S_0) + (mu - 0.5 * sigma^2) * t, sigma^2 * t)
where N denotes a normal distribution with the specified mean and variance.

Discrete-Time Approximation
In practice, to simulate asset prices using GBM, we often discretize the continuous-time model. The discrete approximation of the GBM model is:
S_(t+Delta_t) = S_t * exp((mu - 0.5 * sigma^2) * Delta_t + sigma * sqrt(Delta_t) * Z_t)
where:
Delta_t is the time step.
Z_t is a standard normal random variable, representing the random shock at each time step.



Drift (mu):
Represents the expected return of the asset.
If mu is positive, the asset price has an upward drift over time.

Volatility (sigma):
Represents the uncertainty or risk associated with the asset's returns.
Higher sigma implies greater uncertainty in the asset's future price.

Initial Price (S_0):
The starting value of the asset price at time t = 0.


Option Pricing:
GBM is used in the Black-Scholes model to price European call and put options

Stock Price Modeling:
GBM is a common model for simulating stock prices over time, especially in Monte Carlo simulations.

Risk Management:
GBM is used to assess the risk and return profiles of various financial assets.

GBM captures the random fluctuations in asset prices through the Brownian motion term sigma * S_t * dW_t. The asset price S_t grows exponentially over time, driven by the drift mu. The asset price follows a log-normal distribution, making GBM suitable for modeling financial assets whose prices cannot go negative.
