import ssl
ssl._create_default_https_context = ssl._create_unverified_context

import torch
import esm

# 1. Load a pretrained ESM-2 model (this will download weights if they aren’t already cached)
esm_model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
esm_model.eval()

# 2. Prepare a tiny test sequence (e.g. “ACDEFGHIKLMNPQRSTVWY”)
batch_converter = alphabet.get_batch_converter()
dummy_sequence = [("test", "ACDEFGHIKLMNPQRSTVWY")]
_, _, toks = batch_converter(dummy_sequence)

# 3. Run the model (no gradient needed)
with torch.no_grad():
    out = esm_model(toks, repr_layers=[33])

# 4. Extract the CLS token embedding (represents the whole sequence)
cls_emb = out["representations"][33][0, 0].cpu().numpy()

# 5. Print results
print("ESM-2 model loaded:", esm_model.__class__.__name__)
print("Dummy sequence:", dummy_sequence[0][1])
print("CLS embedding shape:", cls_emb.shape)              # should be (1280,)
print("First 5 dimensions of embedding:", cls_emb[:5])