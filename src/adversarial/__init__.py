# src/adversarial — structured adversarial feedback loop for drug discovery
#
# Replaces the binary Layer 9 → Layer 7 signal with a structured constraint
# violation vector, RL reward shaping, and optional surrogate discriminator.
#
# Public API:
#   rocq_bridge.check_molecule(smiles)           → structured result dict
#   constraint_vector.ConstraintViolationVector   → differentiable penalty signal
#   rl_reward.RewardShaper                       → REINFORCE reward computation
#   metrics.ConstraintMetricsTracker             → rolling-window logging
#   surrogate_discriminator.SurrogateDiscriminator → fast Rocq pre-filter
