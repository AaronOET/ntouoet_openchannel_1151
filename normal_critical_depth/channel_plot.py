# Cross-section plot and ipywidgets controls for the sections in channel_sections

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import ipywidgets as widgets

from channel_sections import Rectangle, Trapezoid, Triangle, Circle

# Slider settings (description, min, max, step) of the geometric parameters
GEOM_SLIDERS = {
    'b0': ('Bottom width b0 (m)', 1.0, 50.0, 0.1),
    's': ('Side slope s (H:1V)', 0.0, 5.0, 0.1),
    'd0': ('Diameter d0 (m)', 0.2, 5.0, 0.05),
}
# Shape-specific overrides of GEOM_SLIDERS and the discharge range (min, max, step)
SHAPE_SLIDERS = {Triangle: {'s': ('Side slope s (H:1V)', 0.1, 5.0, 0.1)}}
Q_RANGE = {Circle: (0.05, 30.0, 0.05)}
Q_RANGE_DEFAULT = (1.0, 500.0, 1.0)


def annotate_dimensions(ax, sec, yn, hb):
    """Top width T, bottom width b0 and side slope s : 1 of the open sections at depth yn
    (as in Ch2NormalDepthWidget_v4)."""
    if isinstance(sec, Circle):
        return
    arrow = dict(arrowstyle='<->', color='black', lw=1)
    bx, _ = sec.boundary(hb)
    dx = 0.01 * (max(bx) - min(bx))          # horizontal text offset

    # Top width, above the bank tops
    xL, xR = sec.water_edges(yn)
    T = sec.topwidth(yn)
    yt = 1.15 * hb
    ax.plot([xL, xL], [yn, yt], color='gray', lw=0.8, ls=':')
    ax.plot([xR, xR], [yn, yt], color='gray', lw=0.8, ls=':')
    ax.annotate('', xy=(xL, yt), xytext=(xR, yt), arrowprops=arrow)
    label = '$T = B$' if isinstance(sec, Rectangle) else '$T$'
    ax.text((xL + xR) / 2, yt + 0.04 * hb, f'{label} = {T:.2f} m', ha='center', va='bottom',
            fontsize=11)

    # Bottom width, below the bed
    if isinstance(sec, Trapezoid):
        yb = -0.12 * hb
        ax.annotate('', xy=(0, yb), xytext=(sec.b0, yb), arrowprops=arrow)
        ax.text(sec.b0 / 2, yb - 0.03 * hb, f'$b_0$ = {sec.b0:.2f} m', ha='center', va='top',
                fontsize=11)
        ax.set_ylim(-0.35 * hb, 1.35 * hb)
    else:
        ax.set_ylim(-0.08 * hb, 1.35 * hb)

    # Side-slope triangle (s : 1) on the right bank, above the water line
    s = getattr(sec, 's', 0)
    if s > 0:
        dv = 0.3 * yn
        y0 = yn + 0.05 * hb
        x0 = sec.water_edges(y0)[1]
        ax.plot([x0, x0 + s * dv, x0 + s * dv], [y0, y0, y0 + dv], color='black', lw=1)
        ax.text(x0 + s * dv / 2, y0 - 0.02 * hb, f'{s:g}', ha='center', va='top', fontsize=10)
        ax.text(x0 + s * dv + dx, y0 + dv / 2, '1', va='center', fontsize=10)


def plot_section(sec, s0, mn, Q):
    """Cross-section of sec with the normal depth y_n, critical depth y_c and a results box."""
    yn, it_n = sec.normal_depth(s0, mn, Q)
    yc, it_c = sec.critical_depth(Q)
    full = np.isnan(yn)
    hb = sec.bank_height(yc if full else yn, yc)

    fig, (ax, axi) = plt.subplots(1, 2, figsize=(11, 4.5),
                                  gridspec_kw={'width_ratios': [2, 1]})
    up = full or yn >= yc          # y_n label above its line and y_c label below, or vice versa
    if not full:
        ax.add_patch(Polygon(sec.water_polygon(yn), facecolor='#8fbbe3', edgecolor='none'))
        xL, xR = sec.water_edges(yn)
        ax.plot([xL, xR], [yn, yn], color='#1f5f99', lw=1.5)
        ax.text((xL + xR) / 2, yn, f'$y_n$ = {yn:.3f} m', ha='center',
                va='bottom' if up else 'top', color='#1f5f99', fontsize=12, fontweight='bold')
    bx, by = sec.boundary(hb)
    ax.plot(bx, by, color='#4a3b26', lw=2)
    xLc, xRc = sec.water_edges(yc)
    ax.plot([xLc, xRc], [yc, yc], color='#c0392b', lw=1.5, ls='--')
    ax.text((xLc + xRc) / 2, yc, f'$y_c$ = {yc:.3f} m', ha='center',
            va='top' if up else 'bottom', color='#c0392b', fontsize=12, fontweight='bold')

    if not full:
        annotate_dimensions(ax, sec, yn, hb)

    if full:
        info = (f'Q = {Q:g} m$^3$/s exceeds the part-full capacity\n'
                f'Qmax = {sec.max_discharge(s0, mn):.3f} m$^3$/s: no normal depth')
    else:
        A, P, T = sec.area(yn), sec.wetperi(yn), sec.topwidth(yn)
        Fr = sec.froude(yn, Q)
        regime = 'subcritical' if Fr < 1 else 'supercritical' if Fr > 1 else 'critical'
        info = (f'Q = {Q:g} m$^3$/s   n = {mn:.3f}   $S_0$ = {s0:g}\n'
                f'A = {A:.3f} m$^2$   P = {P:.3f} m\n'
                f'R = {A / P:.3f} m   T = {T:.3f} m\n'
                f'V = Q/A = {Q / A:.3f} m/s\n'
                f'Fr = {Fr:.3f} ({regime})\n'
                f'Newton iterations: $y_n$ {it_n}, $y_c$ {it_c}')
    axi.axis('off')
    axi.text(0.0, 0.5, info, transform=axi.transAxes, fontsize=10, va='center',
             linespacing=1.6, bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
    sec.style_axes(ax)
    ax.set_xlabel('Horizontal distance (m)')
    ax.set_ylabel('Elevation above invert (m)')
    ax.set_title(f'{sec.name}: normal depth $y_n$ and critical depth $y_c$')
    fig.tight_layout()
    plt.show()


def channel_widget(sec, s0=0.001, mn=0.020, Q=40.0):
    """Cross-section plot of sec with parameter controls below it."""
    shape = type(sec)
    style = {'description_width': '170px'}
    layout = widgets.Layout(width='480px')

    geom = {}
    for p, value in sec.geom.items():
        desc, vmin, vmax, step = SHAPE_SLIDERS.get(shape, {}).get(p, GEOM_SLIDERS[p])
        geom[p] = widgets.FloatSlider(value=value, min=vmin, max=vmax, step=step,
                                      description=desc, style=style, layout=layout)
    # Typed-in bed slope; must stay positive for Manning's equation
    w_s0 = widgets.BoundedFloatText(value=s0, min=1e-6, max=1.0, step=0.0001,
                                    description='Bed slope S0', style=style,
                                    layout=widgets.Layout(width='300px'))
    w_mn = widgets.FloatSlider(value=mn, min=0.010, max=0.080, step=0.001,
                               readout_format='.3f',
                               description="Manning's n", style=style, layout=layout)
    qmin, qmax, qstep = Q_RANGE.get(shape, Q_RANGE_DEFAULT)
    w_Q = widgets.FloatSlider(value=Q, min=qmin, max=qmax, step=qstep,
                              description='Discharge Q (m³/s)', style=style, layout=layout)

    def update(s0, mn, Q, **values):
        plot_section(shape(**values), s0, mn, Q)

    out = widgets.interactive_output(update, {**geom, 's0': w_s0, 'mn': w_mn, 'Q': w_Q})
    controls = widgets.VBox([*geom.values(), w_s0, w_mn, w_Q])
    return widgets.VBox([out, controls])
