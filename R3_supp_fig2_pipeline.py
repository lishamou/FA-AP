

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import os

OUT = "../../results/SuppFig2_feature_engineering_pipeline.pdf"

fig, ax = plt.subplots(figsize=(14, 8))
ax.set_xlim(0, 14)
ax.set_ylim(0, 8)
ax.axis('off')


def draw_box(ax, x, y, w, h, text, color='#A8D5BA'):
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.1",
                          facecolor=color, edgecolor='#2E5E4E', linewidth=2)
    ax.add_patch(box)
    ax.text(x + w/2, y + h/2, text, ha='center', va='center',
            fontsize=10, fontweight='bold')

def draw_arrow(ax, x1, y1, x2, y2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color='#2E5E4E', lw=2))


ax.text(7, 7.6, 'Feature Engineering Pipeline', ha='center',
        fontsize=14, fontweight='bold')


draw_box(ax, 0.5, 6, 2.5, 1, '63 FM-AP modules\n(63 GSVA scores\n+ 3 clinical)', color='#FFE4B5')


draw_box(ax, 4, 6, 2.5, 1, 'Redundancy Removal\n(Pearson |r| > 0.8)', color='#FFD580')
draw_arrow(ax, 3, 6.5, 4, 6.5)


draw_box(ax, 7.5, 6, 2.5, 1, 'Random Forest\nImportance (top N)', color='#FFD580')
draw_arrow(ax, 6.5, 6.5, 7.5, 6.5)


draw_box(ax, 11, 6, 2.5, 1, 'Bootstrap Stability\n(100 iter, >80%)', color='#FFD580')
draw_arrow(ax, 10, 6.5, 11, 6.5)


draw_box(ax, 6.5, 4, 3, 0.8, f'Top N = 10-30 features\n(Selected by grid search)', color='#FFA07A')
draw_arrow(ax, 12.25, 6, 9.25, 4.8)
draw_arrow(ax, 8.75, 6, 8.5, 4.8)


draw_box(ax, 2, 2.5, 3, 0.8, 'Interaction Generation\n(Pairwise products\nof top 10 features)', color='#87CEEB')
draw_arrow(ax, 7.5, 4, 4, 3.3)


draw_box(ax, 6, 2.5, 3, 0.8, 'Top 10 interactions\n(RF importance)', color='#87CEEB')
draw_arrow(ax, 5, 2.9, 6, 2.9)


draw_box(ax, 10, 2.5, 3, 0.8, 'Final: 20 features\n(10 original + 10 interactions)', color='#98FB98')
draw_arrow(ax, 9, 2.9, 10, 2.9)
draw_arrow(ax, 7.5, 4, 11, 3.3)


draw_box(ax, 6, 0.5, 4, 0.8, 'SVM_rbf / Stacking / Voting\n(TCGA training, n=235)', color='#FFB6C1')
draw_arrow(ax, 11.5, 2.5, 9, 1.3)


ax.text(2.6, 5.7, 'Step 1', fontsize=9, color='gray')
ax.text(5.6, 5.7, 'Step 2', fontsize=9, color='gray')
ax.text(9.6, 5.7, 'Step 3-4', fontsize=9, color='gray')
ax.text(8.5, 3.7, 'Step 5', fontsize=9, color='gray')
ax.text(4, 2.2, 'Step 6', fontsize=9, color='gray')
ax.text(7.5, 2.2, 'Step 7', fontsize=9, color='gray')
ax.text(12, 2.2, 'Step 8', fontsize=9, color='gray')
ax.text(8.5, 0.2, 'Step 9', fontsize=9, color='gray')

plt.tight_layout()
plt.savefig(OUT, dpi=200, bbox_inches='tight')
plt.close()
print(f"  → {OUT}")
