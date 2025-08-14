import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("/Users/effieklimi/Documents/projects/pgspilot/pca_model/ref_scores.csv")   # sample, super_pop, PC1…PC4

fig, ax = plt.subplots()
for pop, g in df.groupby("super_pop"):
    ax.scatter(g["PC10"], g["PC11"], s=12, alpha=0.7, label=pop)

ax.set_xlabel("PC10")
ax.set_ylabel("PC11")
ax.legend(title="1000 G super-pop")
fig.tight_layout()
plt.show()

# note: pc6 and pc7 separation need a better check