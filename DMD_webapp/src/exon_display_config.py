"""
Display configuration — re-exports authoritative Dp427m data and map styles.
"""

from __future__ import annotations

from src.dp427m_exon_data import (
    HINGE_GUIDE_AT_EXONS,
    ISOFORM_ANNOTATIONS,
    PHENOTYPE_OPTIONS,
    bar_style_for_phenotype,
    get_domain_map,
    get_exon_table,
)

COLOR_MUTATION_BAR = "#9a9a9a"
COLOR_MUTATION_BAR_BLACK = "#000000"
COLOR_MUTATION_BAR_OUTLINE = "#FFFFFF"

# Cysteine-rich (CR) domain — mutation map palette
COLOR_CR_FILL = "#DAE9F8"
COLOR_CR_BORDER = "#4D93D9"


def pct_dys_color(value: float, vmin: float = 0.0, vmax: float = 340.0) -> str:
    """RdYlGn-style heat color for %Dys(WB) (reference figure)."""
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap("RdYlGn")
    frac = min(max(value, vmin), vmax) / vmax
    r, g, b, _ = cmap(frac)
    return f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}"
