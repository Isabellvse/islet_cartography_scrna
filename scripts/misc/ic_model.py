# ------ Description -------------------------------------------------------------
## This file contrains scripts for our vae model 

# ------ Load library ------------------------------------------------------------
# Model-related libraries
from sklearn.preprocessing import LabelEncoder  # For encoding categorical variables

# AnnData and single-cell analysis libraries
import scanpy as sc       # For preprocessing and visualization of single-cell data
import anndata as ad      # For handling AnnData objects

# Numerical and tensor operations
import numpy as np        # For numerical computations and array manipulations
import pandas as pd
import torch              # Core PyTorch library for tensor operations
import torch.nn as nn     # Neural network components
from torch.optim import Adam  # Optimizer for training models
import torch.nn.functional as F
from torch.distributions import Distribution, Gamma
from torch.distributions import Poisson as PoissonTorch
import warnings
from typing import Optional, Union
from torch.distributions.utils import (
    broadcast_all,
    lazy_property,
    logits_to_probs,
    probs_to_logits,
)
import copy               # For deep copying model states

# ------ Misc Functions ------------------------------------------------------------

def log_nb_positive(
    x: Union[torch.Tensor],
    mu: Union[torch.Tensor],
    theta: Union[torch.Tensor],
    eps: float = 1e-8,
    log_fn: callable = torch.log,
    lgamma_fn: callable = torch.lgamma,
):
    """
    Computes the log-likelihood of the Negative Binomial distribution.
    
    Parameters:
        x (Tensor): Observed counts.
        mu (Tensor): Mean of the distribution.
        theta (Tensor): Dispersion parameter.
        eps (float): Small constant for numerical stability.
        log_fn (callable): Log function.
        lgamma_fn (callable): Log-gamma function.
    
    Returns:
        Tensor: Log-likelihood values.
    """

    log = log_fn
    lgamma = lgamma_fn
    log_theta_mu_eps = log(theta + mu + eps)
    res = (
        theta * (log(theta + eps) - log_theta_mu_eps)
        + x * (log(mu + eps) - log_theta_mu_eps)
        + lgamma(x + theta)
        - lgamma(theta)
        - lgamma(x + 1)
    )

    return res

def _convert_mean_disp_to_counts_logits(mu, theta, eps=1e-6):
    """
    Converts mean and dispersion parameters to total_count and logits for Negative Binomial.

    
    Parameters:
        mu (Tensor): Mean of the distribution.
        theta (Tensor): Dispersion parameter.
        eps (float): Small constant for numerical stability.

    Returns:
        Tuple[Tensor, Tensor]: total_count and logits.

    
    """
    if not (mu is None) == (theta is None):
        raise ValueError(
            "If using the mu/theta NB parameterization, both parameters must be specified"
        )
    logits = (mu + eps).log() - (theta + eps).log()
    total_count = theta
    return total_count, logits


def _convert_counts_logits_to_mean_disp(total_count, logits):
    theta = total_count
    mu = logits.exp() * theta
    return mu, theta


def _gamma(theta, mu):
    concentration = theta
    rate = theta / mu
    # Important remark: Gamma is parametrized by the rate = 1/scale!
    gamma_d = Gamma(concentration=concentration, rate=rate)
    return gamma_d

class NegativeBinomial(Distribution):
    def __init__(
        self,
        total_count: Optional[torch.Tensor] = None,
        probs: Optional[torch.Tensor] = None,
        logits: Optional[torch.Tensor] = None,
        mu: Optional[torch.Tensor] = None,
        theta: Optional[torch.Tensor] = None,
        scale: Optional[torch.Tensor] = None,
        validate_args: bool = False,
    ):
        self._eps = 1e-8
        if (mu is None) == (total_count is None):
            raise ValueError(
                "Please use one of the two possible parameterizations. Refer to the documentation for more information."
            )

        using_param_1 = total_count is not None and (
            logits is not None or probs is not None
        )
        if using_param_1:
            logits = logits if logits is not None else probs_to_logits(probs)
            total_count = total_count.type_as(logits)
            total_count, logits = broadcast_all(total_count, logits)
            mu, theta = _convert_counts_logits_to_mean_disp(total_count, logits)
        else:
            mu, theta = broadcast_all(mu, theta)
        self.mu = mu
        self.theta = theta
        self.scale = scale
        super().__init__(validate_args=validate_args)

    @property
    def mean(self):
        return self.mu

    @property
    def variance(self):
        return self.mean + (self.mean**2) / self.theta

    @torch.inference_mode()
    def sample(
        self,
        sample_shape: Optional[Union[torch.Size, tuple]] = None,
    ) -> torch.Tensor:
        """Sample from the distribution."""
        sample_shape = sample_shape or torch.Size()
        gamma_d = self._gamma()
        p_means = gamma_d.sample(sample_shape)

        # Clamping as distributions objects can have buggy behaviors when
        # their parameters are too high
        l_train = torch.clamp(p_means, max=1e8)
        counts = PoissonTorch(
            l_train
        ).sample()  # Shape : (n_samples, n_cells_batch, n_vars)
        return counts

    def log_prob(self, value: torch.Tensor) -> torch.Tensor:
        return log_nb_positive(value, mu=self.mu, theta=self.theta, eps=self._eps)

    def _gamma(self):
        return _gamma(self.theta, self.mu)

# ------ Encoder ---------------------------------------------------------------
class Encoder(nn.Module):
    def __init__(
        self,
        input_dim: int,
        layer_sizes: list,
        embedding_classes: list,
        embedding_dims: list,
        use_batch_norm: bool = True,
        use_layer_norm: bool = False,
        encoder_bias: bool = True,
        activation_function: nn.Module = nn.LeakyReLU(),
        latent_size: int = 25,
        dropout_rate: float = 0.1,    
     ):
        super().__init__()

        ## Encoder
        encoder_layers = [input_dim + np.sum(embedding_dims)] + layer_sizes
        self.encoder = nn.ModuleList()
        for i in range(len(encoder_layers)-1):
            self.encoder.append(nn.Linear(encoder_layers[i], encoder_layers[i+1], bias = encoder_bias)),
            self.encoder.append(nn.BatchNorm1d(encoder_layers[i+1], affine=True) if use_batch_norm else nn.Identity())
            self.encoder.append(nn.LayerNorm(encoder_layers[i+1], elementwise_affine=False) if use_layer_norm else nn.Identity())
            self.encoder.append(activation_function)
            self.encoder.append(nn.Dropout(p=dropout_rate) if dropout_rate > 0 else nn.Identity())

        # Mean and variance encoder
        self.mean_encoder = nn.Linear(encoder_layers[-1], latent_size)
        self.var_encoder = nn.Linear(encoder_layers[-1], latent_size)
    
    def forward(self, x: torch.Tensor, e: list):
        x = torch.cat([x, e], dim = -1)
        for i in range(len(self.encoder)):
            x = self.encoder[i](x)
        mu = self.mean_encoder(x)
        log_var = self.var_encoder(x)
        return mu, log_var

# ------ Decoder ---------------------------------------------------------------
class Decoder(nn.Module):
    def __init__(
        self,
        input_dim: int,
        layer_sizes: list,
        embedding_classes: list,
        embedding_dims: list,
        use_batch_norm: bool = True,
        use_layer_norm: bool = False,
        decoder_bias: bool = True,
        activation_function: nn.Module = nn.LeakyReLU(),
        latent_size: int = 25,
        dropout_rate: float = 0.1,     
     ):
        super().__init__()

        ## Decoder
        decoder_layers = [latent_size + np.sum(embedding_dims)] + layer_sizes + [input_dim]
        self.decoder = nn.ModuleList()
        for i in range(len(decoder_layers)-1):
            self.decoder.append(nn.Linear(decoder_layers[i], decoder_layers[i+1], bias = decoder_bias)),
            self.decoder.append(nn.BatchNorm1d(decoder_layers[i+1], affine=True) if use_batch_norm & (i <= (len(decoder_layers)-3)) else nn.Identity())
            self.decoder.append(nn.LayerNorm(decoder_layers[i+1], elementwise_affine=False) if use_layer_norm & (i <= (len(decoder_layers)-3)) else nn.Identity())
            self.decoder.append(activation_function if  (i <= (len(decoder_layers)-3)) else nn.Identity())
            self.decoder.append(nn.Dropout(p=dropout_rate) if (dropout_rate > 0) & (i <= (len(decoder_layers)-3)) else nn.Identity())
            #self.decoder.append(nn.Softmax(dim=-1) if i > (len(decoder_layers)-3) else nn.Identity())
    
    def forward(self, x: torch.Tensor, e: list):
        x = torch.cat([x, e], dim = -1)
        for i in range(len(self.decoder)):
            x = self.decoder[i](x)
        return x

# ------ Module ---------------------------------------------------------------
class VAE(nn.Module):
    def __init__(
        self,
        input_dim: int,
        layer_sizes: list,
        embedding_classes: list,
        embedding_dims: list,
        latent_size: int = 25,
        use_batch_norm: bool = True,
        use_layer_norm: bool = False,
        encoder_bias: bool = True,
        decoder_bias: bool = True,
        dropout_rate: float = 0.1,
        activation_function: nn.Module = nn.LeakyReLU(),
    ):
        super().__init__()
        self.dropout_rate = dropout_rate
        
        ## Encoders
        self.encoder = Encoder(input_dim = input_dim, layer_sizes = layer_sizes, embedding_classes = embedding_classes, embedding_dims = embedding_dims, use_batch_norm = use_batch_norm, use_layer_norm = use_layer_norm, encoder_bias = encoder_bias, activation_function = activation_function, latent_size = latent_size, dropout_rate = dropout_rate)

        ## Decoder 
        self.decoder = Decoder(input_dim = input_dim, layer_sizes = layer_sizes, embedding_classes = embedding_classes, embedding_dims = embedding_dims, use_batch_norm = use_batch_norm, use_layer_norm = use_layer_norm, decoder_bias = decoder_bias, activation_function = activation_function, latent_size = latent_size, dropout_rate = dropout_rate)
            
        ## Embedings and conditional
        self.embeddings = nn.ModuleList()
        for i in range(len(embedding_classes)):
            self.embeddings.append(nn.Embedding(embedding_classes[i], embedding_dims[i], max_norm=1))

        ## Distribution parameters
        self.theta = torch.nn.Parameter(torch.randn(input_dim))

    def reparameterize(self, mu, log_var):
        std = torch.exp(0.5 * log_var)
        eps = torch.randn_like(std)
        z = mu + eps * std
        return z

    def forward(self, x: torch.Tensor, study: torch.Tensor, library: torch.Tensor, source: torch.tensor, donor: torch.tensor, use_attention = False, pretrain = False):
        ## Setup
        x_ = x
        x = torch.log1p(x)

        ## Embeddings
        e_source = self.embeddings[0](source)
        e_library = self.embeddings[1](library)
        e_study = self.embeddings[2](study)
        e_donor = self.embeddings[3](donor)
        embeddings = torch.cat([e_source, e_library, e_study, e_donor], dim = -1)
        
        ## Encoding
        mu, log_var = self.encoder(x, embeddings)

        ## Reparameterization trick
        if self.training:
            z = self.reparameterize(mu, log_var)
        else:
            z = mu

        ## Decode
        xhat = F.softmax(self.decoder(z, embeddings), dim=1)
       
        # Distribution parameters
        rate = torch.exp(torch.log(x_.sum(1)).unsqueeze(1)) * xhat
        px = NegativeBinomial(mu=rate, theta=torch.exp(self.theta), scale=xhat)

        return px, rate, xhat, z, log_var, embeddings
        