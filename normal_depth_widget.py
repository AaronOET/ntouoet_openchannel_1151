# Cross-section plot and ipywidgets controls for the trapezoidal-channel normal depth

import matplotlib.pyplot as plt
import ipywidgets as widgets

from section_plot import section_figure


def plot_section(b0, s, s0, mn, Q):
    section_figure(b0, s, s0, mn, Q)
    plt.show()


def channel_widget(b0=20.0, s=2.0, s0=0.001, mn=0.020, Q=40.0):
    """Return the cross-section plot with parameter sliders below it."""
    style = {'description_width': '170px'}
    layout = widgets.Layout(width='480px')
    w_b0 = widgets.FloatSlider(value=b0, min=1.0, max=50.0, step=0.1,
                               description='Bottom width b0 (m)', style=style, layout=layout)
    w_s = widgets.FloatSlider(value=s, min=0.0, max=5.0, step=0.1,
                              description='Side slope s (H:1V)', style=style, layout=layout)
    # Typed-in bed slope; must stay positive for Manning's equation
    w_s0 = widgets.BoundedFloatText(value=s0, min=1e-6, max=1.0, step=0.0001,
                                    description='Bed slope S0', style=style,
                                    layout=widgets.Layout(width='300px'))
    w_mn = widgets.FloatSlider(value=mn, min=0.010, max=0.080, step=0.001,
                               readout_format='.3f',
                               description="Manning's n", style=style, layout=layout)
    w_Q = widgets.FloatSlider(value=Q, min=1.0, max=500.0, step=1.0,
                              description='Discharge Q (m³/s)', style=style, layout=layout)

    out = widgets.interactive_output(plot_section, {'b0': w_b0, 's': w_s, 's0': w_s0,
                                                    'mn': w_mn, 'Q': w_Q})
    controls = widgets.VBox([w_b0, w_s, w_s0, w_mn, w_Q])
    return widgets.VBox([out, controls])
