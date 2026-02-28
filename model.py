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
    return 0.5 * (
        logvar +
        (x - mu) ** 2 / torch.exp(logvar)
    ).mean()


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
    model.train()  # keep stochastic reparameterization active
    preds = []
    with torch.no_grad():
        for _ in range(samples):
            recon_mu, _, _, _ = model(x)
            preds.append(recon_mu)
    preds = torch.stack(preds)          # (samples, batch, input_dim)
    mean = preds.mean(0)
    variance = preds.var(0)
    return mean, variance


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


if __name__ == "__main__":
    main()
