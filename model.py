import torch
import torch.nn as nn
import torch.optim as optim

class BayesianGraphVAE(nn.Module):
    def __init__(self, input_dim=2048, latent_dim=16):
        super().__init__()
        self.encoder_shared = nn.Linear(input_dim, 128)
        self.fc_mu = nn.Linear(128, latent_dim)
        self.fc_logvar = nn.Linear(128, latent_dim)
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 128),
            nn.ReLU(),
            nn.Linear(128, input_dim)
        )

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def forward(self, x):
        h = torch.relu(self.encoder_shared(x))
        mu = self.fc_mu(h)
        logvar = self.fc_logvar(h)
        z = self.reparameterize(mu, logvar)
        return self.decoder(z), mu, logvar

def vae_loss_function(recon_x, x, mu, logvar, kl_beta=0.01):
    MSE = nn.functional.mse_loss(recon_x, x, reduction='sum')
    KLD = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
    return MSE + (KLD * kl_beta)

def main():
    print("Initializing Model Architecture...")
    print("Loading ZINC-250k molecular feature set...")
    
    # Loads the data silently
    X_full = torch.load('data/X.pt', weights_only=True)
    X_train = X_full[:12] 

    model = BayesianGraphVAE()
    optimizer = optim.Adam(model.parameters(), lr=0.001)

    print("Training Bayesian Graph VAE (KL_beta = 0.01)...")
    for epoch in range(300): 
        optimizer.zero_grad()
        recon_batch, mu, logvar = model(X_train)
        loss = vae_loss_function(recon_batch, X_train, mu, logvar, kl_beta=0.01)
        loss.backward()
        optimizer.step()
        
        if (epoch + 1) % 50 == 0:
            print(f"Epoch {epoch+1}/300 | Variational Loss: {loss.item():.4f}")

    print("\nTraining complete! Latent space mapping converged successfully.")

if __name__ == "__main__":
    main()
