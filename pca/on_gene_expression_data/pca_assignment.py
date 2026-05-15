import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# -----------------------------
# LOAD FILES
# -----------------------------
classes = pd.read_csv(
    "class.tsv",
    sep="\t",
    header=None
)
expr = pd.read_csv(
    "filtered.tsv.gz",
    sep="\t",
    compression="gzip"
)

expr.columns = expr.columns.str.strip()

columns = pd.read_csv(
    "columns.tsv.gz",
    sep="\t",
    compression="gzip",
    comment="#"
)

print("\nFILES LOADED SUCCESSFULLY\n")

# -----------------------------
# FIND GENE IDS
# -----------------------------

xbp1_id = columns.loc[
    columns["GeneSymbol"] == "XBP1",
    "ID"
].values[0]

gata3_id = columns.loc[
    columns["GeneSymbol"] == "GATA3",
    "ID"
].values[0]

print("XBP1 ID =", xbp1_id)
print("GATA3 ID =", gata3_id)

# convert IDs to string
xbp1_id = str(xbp1_id).strip()
gata3_id = str(gata3_id).strip()
# -----------------------------
# EXTRACT EXPRESSIONS
# -----------------------------

xbp1 = expr[xbp1_id]
gata3 = expr[gata3_id]

labels = classes.iloc[:, 0]

print("\nExpression extraction complete\n")

# -----------------------------
# FIGURE 1A
# -----------------------------

colors = labels.map({
    0: "black",
    1: "red"
})

plt.figure(figsize=(6,6))

plt.scatter(
    gata3,
    xbp1,
    c=colors
)

plt.xlabel("GATA3")
plt.ylabel("XBP1")

plt.title("Figure 1A")

plt.savefig("figure1a.png")

plt.show()

print("Figure 1A saved")

# -----------------------------
# PCA
# -----------------------------

X = np.column_stack((gata3, xbp1))

pca = PCA(n_components=2)

X_pca = pca.fit_transform(X)

pc1 = X_pca[:, 0]

print("\nPCA completed\n")

print("Explained variance ratio:")
print(pca.explained_variance_ratio_)

# -----------------------------
# FIGURE 1C
# -----------------------------

plt.figure(figsize=(8,4))

# all samples
plt.scatter(
    pc1,
    np.repeat(2, len(pc1)),
    c=colors
)

# ER-
mask0 = labels == 0

plt.scatter(
    pc1[mask0],
    np.repeat(1, np.sum(mask0)),
    c="black"
)

# ER+
mask1 = labels == 1

plt.scatter(
    pc1[mask1],
    np.repeat(0, np.sum(mask1)),
    c="red"
)

plt.yticks(
    [2,1,0],
    ["All", "ER-", "ER+"]
)

plt.xlabel("Projection onto PC1")

plt.title("Figure 1C")

plt.savefig("figure1c.png")

plt.show()

print("Figure 1C saved")
