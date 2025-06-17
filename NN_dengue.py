# Required Libraries
!pip install torchdiffeq

import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib.pyplot as plt
from torchdiffeq import odeint
from google.colab import files
np.random.seed(42)
# Parameters
params = {
    'Pi_H': 100, 'Pi_V': 3000,
    'eta1': 1.4, 'eta2': 28,
    'mu_H': 1/(65*365), 'mu_V': 100/365,
    'phi1': 73/365, 'phi2': 73/365,
    'gamma1': 52/365, 'gamma2': 52/365,
    'delta1': 0.18/365, 'delta2': 0.18/365,
    'deltaV1': 0, 'deltaV2': 0,
    'w1': 12/365, 'w2': 12/365,
    'alphaV1': 1000/365, 'alphaV2': 800/365,
    'alphaH1': 10/365, 'alphaH2': 10/365,
    'waning': 0.005,  # waning immunity rate
    'u1': 0.6, 'u2': 0.3, 'epsilon':0.8,
    'W': 100
}

# Initial conditions
y0 = [1e3, 1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e2, 1e4, 1e2, 1e2, 1e2]  # Add V_H as last

# Time setup
T = 100
t_points = 1000
t = torch.linspace(0, T, t_points)

# Control network for phi(t)
class ControlNN(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(1, 128),
            nn.Tanh(),
            nn.Linear(128, 128),
            nn.Tanh(),
            nn.Linear(128, 1),
            nn.Sigmoid()
        )
    def forward(self, t):
        return 0.7 * self.net(t.view(-1, 1)).squeeze()

# Full model with vaccinated host V_H and waning immunity
class FullModel(nn.Module):
    def __init__(self, control_nn, params):
        super().__init__()
        self.control_nn = control_nn
        self.p = params

    def forward(self, t, y):
        (S, E1, I1, R1, C1, E12, I12, E2, I2, R2, C2, E21, I21, R,
         Sv, VI1, VI2, VH) = y

        # Force of infection
# Total host population
        NH = S + E1 + I1 + R1 + C1 + E12 + I12 + E2 + I2 + R2 + C2 + E21 + I21 + R  + Sv + VI1 + VI2 + VH

        # Updated force of infection
        lambdaV1 = self.p['alphaV1'] * VI1 / NH
        lambdaV2 = self.p['alphaV2'] * VI2 / NH

        lambdaH1 = self.p['alphaH1'] * (I1 + self.p['eta1'] * I21) / NH
        lambdaH2 = self.p['alphaH2'] * (I2 + self.p['eta2'] * I12) / NH

        phi = self.control_nn(t)

        dS = self.p['Pi_H'] - (lambdaV1 + lambdaV2) * S - phi * S + self.p['waning'] * VH - self.p['mu_H'] * S
        dE1 = lambdaV1 * S - (self.p['phi1'] + self.p['mu_H']) * E1 + (1-self.p['epsilon'])* lambdaV1 * VH
        dE2 = lambdaV2 * S - (self.p['phi2'] + self.p['mu_H']) * E2 + (1-self.p['epsilon'])* lambdaV2 * VH
        dI1 = self.p['phi1'] * E1 - (self.p['gamma1'] + self.p['delta1'] + self.p['mu_H']) * I1
        dI2 = self.p['phi2'] * E2 - (self.p['gamma2'] + self.p['delta2'] + self.p['mu_H']) * I2
        dR1 = self.p['gamma1'] * I1 - (self.p['w1'] + self.p['mu_H']) * R1
        dR2 = self.p['gamma2'] * I2 - (self.p['w2'] + self.p['mu_H']) * R2
        dC1 = self.p['w1'] * R1 - (self.p['u1'] * lambdaV2 + (1 - self.p['u1']) + self.p['mu_H']) * C1
        dC2 = self.p['w2'] * R2 - (self.p['u2'] * lambdaV1 + (1 - self.p['u2']) + self.p['mu_H']) * C2
        dE21 = self.p['u2'] * lambdaV1 * C2 - (self.p['phi1'] + self.p['mu_H']) * E21
        dE12 = self.p['u1'] * lambdaV2 * C1 - (self.p['phi2'] + self.p['mu_H']) * E12
        dI21 = self.p['phi1'] * E21 - (self.p['gamma1'] + self.p['delta1'] + self.p['mu_H']) * I21
        dI12 = self.p['phi2'] * E12 - (self.p['gamma2'] + self.p['delta2'] + self.p['mu_H']) * I12
        dR = self.p['gamma1'] * I21 + self.p['gamma2'] * I12 + (1 - self.p['u1']) * C1 + (1 - self.p['u2']) * C2 - self.p['mu_H'] * R
        dSv = self.p['Pi_V'] - (lambdaH1 + lambdaH2) * Sv - self.p['mu_V'] * Sv
        dVI1 = lambdaH1 * Sv - (self.p['deltaV1'] + self.p['mu_V']) * VI1
        dVI2 = lambdaH2 * Sv - (self.p['deltaV2'] + self.p['mu_V']) * VI2
        dVH = phi * S -(1-self.p['epsilon'])*(lambdaV1 + lambdaV2) * VH- self.p['waning'] * VH - self.p['mu_H'] * VH

        return torch.stack([dS, dE1, dI1, dR1, dC1, dE12, dI12, dE2, dI2, dR2, dC2,
                            dE21, dI21, dR, dSv, dVI1, dVI2, dVH])

# Model & optimizer
control_nn = ControlNN()
model = FullModel(control_nn, params)
optimizer = optim.Adam(control_nn.parameters(), lr=0.001)

# Training
y0_tensor = torch.tensor(y0, dtype=torch.float32)
losses = []

for epoch in range(2000):
    optimizer.zero_grad()
    y = odeint(model, y0_tensor, t)
    infected = y[:, 2] + y[:, 6] + y[:, 8] + y[:, 12] # All infected compartments
    phi_t = control_nn(t)
    loss = torch.mean(infected + 0.5 * params['W'] * phi_t**2)
    loss.backward()
    optimizer.step()
    losses.append(loss.item())
    if epoch % 100 == 0:
        print(f"Epoch {epoch}: Loss = {loss.item():.4f}")
# Final simulation
y_pred = odeint(model, y0_tensor, t).detach().numpy()
phi_control = control_nn(t).detach().numpy()
t_np = t.detach().numpy()
# Functional J approximation using trapezoidal rule after final simulation
infected_total = y_pred[:, 2] + y_pred[:, 6] + y_pred[:, 8] + y_pred[:, 12]  # I1 + I12 + I2 + I21
phi_sq = 0.5 * params['W'] * phi_control**2
integrand = infected_total + phi_sq

# Use trapezoidal rule to approximate integral over time
dt = t_np[1] - t_np[0]
functional_J = np.trapezoid(integrand, dx=dt)

print(f"\nFinal Functional Value J ≈ {functional_J:.4f}")
# Plotting
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
# Extracting compartments
I1 = y_pred[:, 2]
I21 = y_pred[:, 12]
I2 = y_pred[:, 8]
I12 = y_pred[:, 6]
plt.plot(t,  I1 + I21, label="(I1 + I21)", color='red')
plt.plot(t,  I2 + I12, label="(I2 + I12)", color='blue')
#plt.plot(t,  np.log10(I1 + I21), label="log(I1 + I21)", color='red')
#plt.plot(t,  np.log10(I2 + I12), label="log(I2 + I12)", color='blue')
plt.xlabel("Time/Days")
plt.ylabel("Infected Population")
plt.title("Infected Compartments by Strain")
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(t, phi_control, label="Optimal Control ϕ(t)", color='green')
plt.xlabel("Time/Days")
plt.ylabel("ϕ(t)")
plt.title("Vaccination Control")
plt.legend()
plt.tight_layout()
plt.show()
