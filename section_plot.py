# Cross-section figure for the trapezoidal-channel normal depth (shared by the
# ipywidgets notebook and the Streamlit app)

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

from normal_depth import area, wetperi, topwidth, normal_depth, critical_depth, froude


def section_figure(b0, s, s0, mn, Q):
    """Return (fig, results) for the channel at normal depth."""
    yn, iters = normal_depth(b0, s, s0, mn, Q)
    A = area(yn, b0, s)
    P = wetperi(yn, b0, s)
    T = topwidth(yn, b0, s)
    V = Q / A
    yc, _ = critical_depth(b0, s, Q)
    Fr = froude(yn, b0, s, Q)
    regime = 'subcritical' if Fr < 1 else 'supercritical' if Fr > 1 else 'critical'
    hbank = 1.5 * max(yn, yc)                      # bank height above the bed
    pad = 0.12 * (b0 + 2 * s * hbank)      # flat ground beyond each bank
    xLt, xRt = -s * hbank, b0 + s * hbank  # bank tops
    xLw, xRw = -s * yn, b0 + s * yn        # water-surface edges
    ybot = -0.35 * hbank

    fig, ax = plt.subplots(figsize=(11, 4.5))

    ax.add_patch(Polygon([(xLt - pad, hbank), (xLt, hbank), (0, 0), (b0, 0),
                          (xRt, hbank), (xRt + pad, hbank), (xRt + pad, ybot),
                          (xLt - pad, ybot)], facecolor='#e8dcc4', edgecolor='none'))
    ax.add_patch(Polygon([(xLw, yn), (0, 0), (b0, 0), (xRw, yn)],
                         facecolor='#8fbbe3', edgecolor='none'))
    ax.plot([xLt - pad, xLt, 0, b0, xRt, xRt + pad],
            [hbank, hbank, 0, 0, hbank, hbank], color='#4a3b26', lw=2)
    ax.plot([xLw, xRw], [yn, yn], color='#1f5f99', lw=1.5)

    # Water-surface symbol
    xc = b0 / 2
    ax.plot(xc, yn + 0.07 * hbank, marker='v', color='#1f5f99', ms=10)
    ax.plot([xc - 0.03 * b0, xc + 0.03 * b0], [yn - 0.05 * hbank] * 2, color='#1f5f99', lw=1)
    ax.plot([xc - 0.018 * b0, xc + 0.018 * b0], [yn - 0.10 * hbank] * 2, color='#1f5f99', lw=1)

    arrow = dict(arrowstyle='<->', color='black', lw=1)

    # Bottom width
    yb = -0.12 * hbank
    ax.annotate('', xy=(0, yb), xytext=(b0, yb), arrowprops=arrow)
    ax.text(b0 / 2, yb - 0.03 * hbank, f'$b_0$ = {b0:.2f} m', ha='center', va='top', fontsize=11)

    # Top width
    yt = hbank + 0.15 * hbank
    ax.plot([xLw, xLw], [yn, yt], color='gray', lw=0.8, ls=':')
    ax.plot([xRw, xRw], [yn, yt], color='gray', lw=0.8, ls=':')
    ax.annotate('', xy=(xLw, yt), xytext=(xRw, yt), arrowprops=arrow)
    ax.text(b0 / 2, yt + 0.04 * hbank, f'$T$ = {T:.2f} m', ha='center', va='bottom', fontsize=11)

    # Normal depth
    xd = 0.8 * b0
    ax.annotate('', xy=(xd, 0), xytext=(xd, yn), arrowprops=arrow)
    ax.text(xd + 0.015 * b0, yn / 2, f'$y_n$ = {yn:.3f} m', va='center', fontsize=12,
            fontweight='bold')

    # Critical depth
    xLc, xRc = -s * yc, b0 + s * yc
    ax.plot([xLc, xRc], [yc, yc], color='#c0392b', lw=1.5, ls='--')
    ax.text(xLc + 0.02 * b0, yc + 0.02 * hbank, f'$y_c$ = {yc:.3f} m', ha='left',
            va='bottom', fontsize=12, color='#c0392b', fontweight='bold')

    # Side-slope triangle (s : 1) on the right bank, above the water line
    if s > 0:
        dv = 0.3 * yn
        y0 = yn + 0.05 * hbank
        x0 = b0 + s * y0
        ax.plot([x0, x0 + s * dv, x0 + s * dv], [y0, y0, y0 + dv], color='black', lw=1)
        ax.text(x0 + s * dv / 2, y0 - 0.02 * hbank, f'{s:g}', ha='center', va='top', fontsize=10)
        ax.text(x0 + s * dv + 0.01 * b0, y0 + dv / 2, '1', va='center', fontsize=10)

    # Results box
    info = (f'Q = {Q:g} m$^3$/s   n = {mn:.3f}   $S_0$ = {s0:g}\n'
            f'A = {A:.2f} m$^2$   P = {P:.2f} m   R = {A / P:.3f} m\n'
            f'V = Q/A = {V:.3f} m/s   Newton iterations = {iters}\n'
            f'$y_c$ = {yc:.3f} m   Fr = {Fr:.3f} ({regime})')
    ax.text(0.01, 0.02, info, transform=ax.transAxes, fontsize=10, va='bottom',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

    xspan = (xRt + pad) - (xLt - pad)
    yspan = (hbank + 0.35 * hbank) - (ybot - 0.45 * hbank)
    ax.set_xlim(xLt - pad, xRt + pad)
    ax.set_ylim(ybot - 0.45 * hbank, hbank + 0.35 * hbank)
    ax.set_xlabel('Horizontal distance (m)')
    ax.set_ylabel('Elevation above bed (m)')

    # Report the vertical exaggeration implied by the axes shape
    bbox = ax.get_position()
    w_in, h_in = fig.get_size_inches()
    ve = (bbox.height * h_in / yspan) / (bbox.width * w_in / xspan)
    ax.set_title(f'Trapezoidal channel at normal depth  (vertical exaggeration ≈ {ve:.1f}×)')

    results = dict(yn=yn, yc=yc, Fr=Fr, regime=regime, iters=iters, A=A, P=P, R=A / P, T=T, V=V)
    return fig, results
