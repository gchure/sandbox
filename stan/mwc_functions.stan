functions {
    real prob_act(real c, real ep_a, real ep_i, real ep_ai, real n_sites) {
        real numer;
        real denom;
        numer = (1 + c * exp(-ep_a))^n_sites;
        denom = numer + exp(-ep_ai) * (1 + c * exp(-ep_i))^n_sites;
        return numer / denom;}

    real foldchange(real pact, real R, real Nns, real ep_r) {
        return (1 + pact * (R / Nns) * exp(-ep_r))^-1;}
  }
