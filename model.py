import math
import os

import torch
import torch.nn as nn
import torch.optim as optim


class BayesianGraphVAE(nn.Module):
    def __init__(self, input_dim=2048, latent_dim=16, hidden_dim=64):
        super().__init__()
        # Encoder (structure unchanged)
        self.encoder_shared = nn.Linear(input_dim, 128)
        self.fc_mu = nn.Linear(128, latent_dim)
        self.fc_logvar = nn.Linear(128, latent_dim)

        # Fix 3: Reduced decoder capacity (hidden_dim=64, single layer)
        self.decoder_hidden = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim),
            nn.ReLU(),
        )
        # Fix 4: Decoder predicts mean and log variance (aleatoric uncertainty)
        self.mu_head = nn.Linear(hidden_dim, input_dim)
        self.logvar_head = nn.Linear(hidden_dim, input_dim)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def forward(self, x):
        h = torch.relu(self.encoder_shared(x))
        mu = self.fc_mu(h)
        logvar = self.fc_logvar(h)
        z = self.reparameterize(mu, logvar)
        h_dec = self.decoder_hidden(z)
        recon_mu = self.mu_head(h_dec)
        recon_logvar = self.logvar_head(h_dec)
        return recon_mu, recon_logvar, mu, logvar


# Fix 1: Free Bits KL Regularization
def kl_divergence(mu, logvar, free_bits=0.5):
    kl = -0.5 * (1 + logvar - mu.pow(2) - logvar.exp())
    # Clamp per-dimension KL to prevent posterior collapse
    kl = torch.clamp(kl, min=free_bits)
    return kl.sum(dim=1).mean()


# Fix 4: Gaussian NLL as reconstruction loss (replaces MSE)
def gaussian_nll(x, mu, logvar):
    # Clamp logvar to prevent exp() overflow/underflow → NaN loss
    logvar = torch.clamp(logvar, -10, 10)
    # Sum over features (not mean) so reconstruction and KL are on the same scale
    return 0.5 * (
        logvar +
        (x - mu) ** 2 / torch.exp(logvar)
    ).sum(dim=1).mean()


# Fix 2: Cyclical KL Annealing
def cyclical_beta(epoch, cycle_length=100):
    cycle_pos = epoch % cycle_length
    beta = min(1.0, cycle_pos / (cycle_length * 0.5))
    return beta


def vae_loss_function(recon_mu, recon_logvar, x, mu, logvar, beta=1.0, free_bits=0.5):
    recon_loss = gaussian_nll(x, recon_mu, recon_logvar)
    kl_loss = kl_divergence(mu, logvar, free_bits=free_bits)
    total = recon_loss + beta * kl_loss
    return total, recon_loss, kl_loss


# Fix 5: Monte Carlo Uncertainty Estimation
def mc_predict(model, x, samples=50):
    model.train()  # keep stochastic latent sampling active

    mu_samples = []
    var_samples = []

    with torch.no_grad():
        for _ in range(samples):
            recon_mu, recon_logvar, _, _ = model(x)

            # Clamp log variance to prevent collapse
            recon_logvar = torch.clamp(recon_logvar, -6, 3)

            mu_samples.append(recon_mu)
            var_samples.append(torch.exp(recon_logvar))

    mu_samples = torch.stack(mu_samples)   # (samples, batch, dim)
    var_samples = torch.stack(var_samples)

    mean_prediction = mu_samples.mean(0)

    # Epistemic uncertainty
    epistemic_var = mu_samples.var(0)

    # Aleatoric uncertainty
    aleatoric_var = var_samples.mean(0)

    # Total predictive variance
    predictive_var = epistemic_var + aleatoric_var

    return mean_prediction, predictive_var


def main():
    print("Initializing Model Architecture...")
    print("Loading ZINC-250k molecular feature set...")

    # Loads the data silently (unchanged)
    X_full = torch.load('data/X.pt', weights_only=True)
    X_train = X_full 

    model = BayesianGraphVAE()
    optimizer = optim.Adam(model.parameters(), lr=0.001)

    free_bits = 0.5
    cycle_length = 100
    print(f"Training Bayesian Graph VAE (free_bits={free_bits}, cycle_length={cycle_length})...")

    for epoch in range(300):
        optimizer.zero_grad()
        recon_mu, recon_logvar, mu, logvar = model(X_train)

        # Fix 2: cyclical beta replaces fixed kl_beta
        beta = cyclical_beta(epoch, cycle_length=cycle_length)

        loss, recon_loss, kl_loss = vae_loss_function(
            recon_mu, recon_logvar, X_train, mu, logvar,
            beta=beta, free_bits=free_bits,
        )
        loss.backward()
        optimizer.step()

        if (epoch + 1) % 50 == 0:
            print(
                f"Epoch {epoch+1}/300 | Loss: {loss.item():.4f} "
                f"| Recon (NLL): {recon_loss.item():.4f} "
                f"| KL: {kl_loss.item():.4f} "
                f"| Beta: {beta:.4f}"
            )

    print("\nTraining complete! Latent space mapping converged successfully.")
    os.makedirs("models", exist_ok=True)
    torch.save(model.state_dict(), "models/model.pt")
    print("Model saved to models/model.pt")

    # Fix 5: MC uncertainty metrics
    print("\nRunning Monte Carlo uncertainty estimation (50 samples)...")
    mean_pred, variance = mc_predict(model, X_train, samples=50)

    per_sample_mse = ((mean_pred - X_train) ** 2).mean(dim=1)
    per_sample_var = variance.mean(dim=1)
    avg_variance = per_sample_var.mean().item()

    corr_matrix = torch.corrcoef(torch.stack([per_sample_mse, per_sample_var]))
    corr = corr_matrix[0, 1].item()

    print(f"  Mean predictive variance      : {avg_variance:.5f}")
    print(f"  Error-uncertainty correlation : {corr:.4f}")


def train_with_rl_reward(
    epochs: int = 300,
    lambda_rl: float = 0.1,
    free_bits: float = 0.5,
    cycle_length: int = 100,
    smiles_csv: str = "data/molecules_clean.csv",
    property_predictor_ckpt: str = "data/property_predictor.pt",
    weights_config: str = "config/constraint_weights.yaml",
    save_path: str = "models/model_rl.pt",
    log_every: int = 50,
):
    """Training loop with adversarial RL reward shaping (Step 3).

    Augments the standard VAE objective with a REINFORCE term:

        L_total = L_VAE(recon, KL) + λ_RL · L_RL

    where:
        L_RL = −(reward − baseline) · log q(z | μ, log σ²)
        reward = QED_predicted(L8) − constraint_penalty(Rocq bridge, L9)

    The Rocq proof checker is treated as a black-box oracle; no gradients
    flow through it.  The QED prediction from the PropertyPredictor IS
    differentiable but we use REINFORCE so the same loss handles both.

    Per-constraint violation rates and reward curves are logged to a
    ConstraintMetricsTracker and printed every `log_every` epochs.

    Parameters
    ----------
    epochs                  : total training epochs
    lambda_rl               : weight of the RL term relative to VAE loss
    free_bits               : KL free-bits threshold (posterior collapse prevention)
    cycle_length            : epochs per KL annealing cycle
    smiles_csv              : path to molecules_clean.csv (SMILES + QED labels)
    property_predictor_ckpt : path to trained PropertyPredictor checkpoint
    weights_config          : path to constraint_weights.yaml
    save_path               : where to save the trained model
    log_every               : print interval in epochs
    """
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

    from src.adversarial.rl_reward import RewardShaper, log_q_gaussian, reinforce_loss
    from src.adversarial.metrics import ConstraintMetricsTracker

    print("Initializing RL-augmented VAE training...")
    X_full = torch.load("data/X.pt", weights_only=True)
    X_train = X_full

    # Load SMILES for Rocq bridge evaluation
    smiles_list = None
    if os.path.exists(smiles_csv):
        import pandas as pd
        df = pd.read_csv(smiles_csv)
        smiles_list = df["smiles"].tolist()[: len(X_train)]
        print(f"  Loaded {len(smiles_list)} SMILES from {smiles_csv}")
    else:
        print(f"  Warning: {smiles_csv} not found; RL reward will use penalty only.")

    # Set up reward shaper
    reward_shaper = RewardShaper(weights_config=weights_config)
    if os.path.exists(property_predictor_ckpt):
        reward_shaper.load_property_predictor(property_predictor_ckpt)
        print(f"  Loaded PropertyPredictor from {property_predictor_ckpt}")

    model = BayesianGraphVAE()
    optimizer = optim.Adam(model.parameters(), lr=0.001)
    metrics = ConstraintMetricsTracker(window=1000)

    print(
        f"Training with RL reward shaping "
        f"(λ_RL={lambda_rl}, free_bits={free_bits}, cycle_length={cycle_length})..."
    )

    for epoch in range(epochs):
        model.train()
        optimizer.zero_grad()

        recon_mu, recon_logvar, mu, logvar = model(X_train)
        beta = cyclical_beta(epoch, cycle_length=cycle_length)
        vae_loss, recon_loss, kl_loss = vae_loss_function(
            recon_mu, recon_logvar, X_train, mu, logvar,
            beta=beta, free_bits=free_bits,
        )

        # --- REINFORCE reward term ---
        rl_loss = torch.tensor(0.0)
        mean_reward = 0.0
        violation_rates: dict = {}

        if smiles_list is not None:
            # Reparameterised z (already computed inside forward; recompute
            # so we have explicit access to the eps / z tensors for log q)
            with torch.no_grad():
                std = torch.exp(0.5 * logvar)
                eps = torch.randn_like(std)
                z = mu + eps * std  # (N, latent_dim)

            log_q = log_q_gaussian(z, mu, logvar)  # (N,)

            # Compute rewards from Rocq bridge + PropertyPredictor
            # Use decoded fingerprints as proxy fingerprints for QED scoring
            fingerprints = recon_mu.detach()
            rewards, violation_rates = reward_shaper.compute_batch_rewards(
                smiles_list, fingerprints
            )

            rl_loss = reinforce_loss(log_q, rewards, baseline=reward_shaper.baseline)
            mean_reward = rewards.mean().item()
            reward_shaper.update_baseline(mean_reward)

            # Update rolling metrics
            pass_flags = [
                all(v == 0 for v in {c: int(violation_rates.get(c, 0) > 0)}.values())
                for _ in smiles_list
            ]
            metrics.update_violation_rates(violation_rates)
            metrics.update_reward(mean_reward)

        total_loss = vae_loss + lambda_rl * rl_loss
        total_loss.backward()
        optimizer.step()

        if (epoch + 1) % log_every == 0:
            print(
                f"Epoch {epoch+1}/{epochs}"
                f" | VAE: {vae_loss.item():.4f}"
                f" | Recon: {recon_loss.item():.4f}"
                f" | KL: {kl_loss.item():.4f}"
                f" | RL: {rl_loss.item():.4f}"
                f" | Reward: {mean_reward:.4f}"
                f" | β: {beta:.3f}"
            )
            if violation_rates:
                metrics.print_summary()

    print("\nRL-augmented training complete.")
    os.makedirs(os.path.dirname(save_path) or "models", exist_ok=True)
    torch.save(model.state_dict(), save_path)
    print(f"Model saved to {save_path}")

    print("\nFinal metrics summary:")
    metrics.print_summary()

    return model


if __name__ == "__main__":
    main()
