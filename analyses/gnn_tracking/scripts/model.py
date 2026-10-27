"""Interaction network edge classifier, following DeZoort et al., arXiv:2103.16701, Fig. 6.

  relational model  a'_ij = phi_R1(x_i, x_j, a_ij)            MLP(hidden, hidden) -> msg_dim (4)
  aggregation       m_i   = sum over edges j->i of a'_ij
  object model      x'_i  = phi_O(x_i, m_i)                   MLP(hidden, hidden) -> node_dim_out (3)
  edge weighting    w_ij  = phi_R2(x'_i, x'_j, a'_ij)         MLP(hidden, hidden) -> 1

With iterations > 1 the relational + object step is repeated on (x', a') before
the edge weighting (not in the paper; iterations=1 is the paper's model).
forward() returns logits; apply a sigmoid for the edge weight w in [0, 1].
"""
import torch
from torch import nn


def mlp(n_in, hidden, n_out):
    return nn.Sequential(nn.Linear(n_in, hidden), nn.ReLU(), nn.Linear(hidden, hidden), nn.ReLU(),
                         nn.Linear(hidden, n_out))


class InteractionNetwork(nn.Module):
    def __init__(self, node_dim, edge_dim, hidden=40, msg_dim=4, node_out=3, iterations=1):
        super().__init__()
        self.relational = nn.ModuleList([mlp(2 * node_dim + edge_dim, hidden, msg_dim)])
        self.object = nn.ModuleList([mlp(node_dim + msg_dim, hidden, node_out)])
        for _ in range(iterations - 1):
            self.relational.append(mlp(2 * node_out + msg_dim, hidden, msg_dim))
            self.object.append(mlp(node_out + msg_dim, hidden, node_out))
        self.edge_weight = mlp(2 * node_out + msg_dim, hidden, 1)

    def forward(self, x, edge_index, edge_attr):
        src, dst = edge_index  # messages flow src -> dst (dst is the receiver i)
        a = edge_attr
        for relational, obj in zip(self.relational, self.object):
            a = relational(torch.cat([x[dst], x[src], a], dim=1))
            m = torch.zeros(x.shape[0], a.shape[1], dtype=x.dtype, device=x.device).index_add_(0, dst, a)
            x = obj(torch.cat([x, m], dim=1))
        return self.edge_weight(torch.cat([x[dst], x[src], a], dim=1)).squeeze(1)


def make_bidirectional(edge_index, edge_attr, signed_cols):
    """Append each edge reversed, flipping the sign of the directional features.

    The paper's graphs are bidirectional so that every stub receives messages from
    both its inner and outer neighbours; the score of an edge is then the minimum
    of its two directions (as in the paper's track-building code).
    """
    rev = edge_attr.clone()
    rev[:, signed_cols] = -rev[:, signed_cols]
    return torch.cat([edge_index, edge_index.flip(0)], dim=1), torch.cat([edge_attr, rev], dim=0)


def edge_scores(model, x, edge_index, edge_attr, signed_cols):
    """Per-edge weight in [0, 1] for the original (one-directional) edges."""
    ei, ea = make_bidirectional(edge_index, edge_attr, signed_cols)
    w = torch.sigmoid(model(x, ei, ea))
    n = edge_index.shape[1]
    return torch.minimum(w[:n], w[n:])
