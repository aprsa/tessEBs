import numpy as np

from .models import Ephemeris

from . import pipeline as pl

from bokeh.themes import Theme
from bokeh import events
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Div, Spinner, Select, Button, HelpButton, InlineStyleSheet, Tooltip, Dialog, TextInput
from bokeh.models.callbacks import CustomJS
from bokeh.layouts import row, column
from bokeh.embed import components


def create_ephemeris_ui(tess_id):
    # Data loading:
    time, flux, ferr = pl.read_from_ascii(tess_id, content='lc')

    # TODO: move this to css if possible:
    # theme = Theme(json={'attrs': {'Title': {'text_color': 'red'}}})
    theme = Theme('catalog/ephemeris_theme.yaml')

    stylesheet = InlineStyleSheet(
        css="""
        .bk-btn {
            background-color: #94b7d7;
            color: black;
            font-size: 11pt;
            font-family: times;
            # font-weight: bold;
            border: 1px solid black;
            border-radius: 4px;
            padding: 4px 6px 6px 6px;
            margin: 2px;
        }
        .bk-btn:hover {
            background-color: #4972ad;
            color: black;
            font-size: 11pt;
            font-family: times;
            border: 1px solid black;
            border-radius: 4px;
            padding: 4px 6px 6px 6px;
            margin: 0px;
        }
        """
    )

    # original, unmodified data:
    lc_data = ColumnDataSource({
        'time': time,
        'phase': pl.bjd2phase(time, 0, 1),
        'flux': flux,
        'dflux': np.diff(flux, append=[flux[-1],]),
    })

    # data that will be displayed:
    displayed_data = ColumnDataSource({
        'time': time,
        'phase': pl.bjd2phase(time, 0, 1),
        'flux': flux,
        'dflux': np.diff(flux, append=[flux[-1],]),
    })

    # periodogram data:
    freq, power, lsfap = pl.read_from_ascii(tess_id, content='spd')
    spd_data = ColumnDataSource({
        'period': 2*np.pi/freq,
        'power': power,
    })

    # figure definitions:
    lc = figure(
        tools='pan,box_zoom,wheel_zoom,save,reset',  # also available: help, hover, box_select, lasso_select, poly_select, tap, crosshair
        active_drag='box_zoom',
        title=f'TIC {tess_id} timeseries',
        x_axis_label='Time [days]',
        y_axis_label='Normalized flux',
        sizing_mode='stretch_width',
        height=300,
    )
    lc.toolbar.logo = None

    spd = figure(
        tools='crosshair,pan,box_zoom,reset',
        active_drag='box_zoom',
        title=f'TIC {tess_id} periodogram',
        x_axis_label='Period [days]',
        y_axis_label='Power',
        sizing_mode='stretch_width',
        height=300,
    )
    spd.toolbar.logo = None

    ph = figure(
        tools='pan,box_zoom,wheel_zoom,save,reset',
        active_drag='box_zoom',
        title=f'TIC {tess_id} phase plot',
        x_axis_label='Phase',
        y_axis_label='Normalized flux',
        sizing_mode='stretch_width',
        height_policy='fixed',
        width=45,
        height=300,
    )
    ph.toolbar.logo = None

    # Callback definitions:
    def switch_yaxis(args):
        callback = CustomJS(
            args=args,
            code="""
                const yaxis = cb_obj.value;
                const data = source.data;

                plots.forEach(plot => {
                    if (yaxis === 'flux') {
                        plot.glyph.y = {'field': 'flux'};
                    } else {
                        plot.glyph.y = {'field': 'dflux'};
                    }
                });

                source.change.emit();
            """
        )
        return callback

    def rephase(args):
        callback = CustomJS(
            args=args,
            code="""
                const t0 = t0_widget.value;
                const P = period_widget.value;

                function mod(n, m) {
                    return ((n % m) + m) % m;
                }

                sources.forEach(source => {
                    const phase = source.data['time'].map(t => -0.5 + mod(t - t0 - 0.5 * P, P) / P);
                    source.data['phase'] = phase;
                    source.change.emit();
                });
            """
        )
        return callback

    def apply_ephem(args):
        callback = CustomJS(
            args=args,
            code="""
                t0_widget.value = t0;
                period_widget.value = period;
            """
        )
        return callback

    def thin_data(args):
        callback = CustomJS(
            args=args,
            code="""
                const thinning = cb_obj.value;
                const full = full_data.data;
                const display = displayed_data.data;

                Object.keys(display).forEach(key => {
                    display[key] = [];
                });

                for (let i = 0; i < full.time.length; i += thinning) {
                    Object.keys(display).forEach(key => {
                        display[key].push(full[key][i]);
                    });
                }

                displayed_data.change.emit();
            """
        )
        return callback

    # widget definitions:
    t0_widget = Spinner(
        title='t0 [days]:',
        low=0,
        high=1e10,
        step=time[1]-time[0],
        value=0.0,
        format='0.00000000',
    )

    period_widget = Spinner(
        title='Period [days]:',
        low=time[1]-time[0],
        high=(time[-1]-time[0])/2,
        step=time[1]-time[0],
        value=1.0,
        format='1.00000000',
    )

    half_period_button = Button(
        label='Half period',
        button_type='default',
        width=120,
        height=35,
    )
    half_period_button.js_on_click(CustomJS(
        args=dict(period=period_widget),
        code="""
            period.value /= 2;
        """
    ))

    double_period_button = Button(
        label='Double period',
        button_type='default',
        width=120,
        height=35,
    )
    double_period_button.js_on_click(CustomJS(
        args=dict(period=period_widget),
        code="""
            period.value *= 2;
        """
    ))

    harmonic_row = row(
        children=[
            half_period_button,
            double_period_button,
        ],
        sizing_mode='stretch_width',
    )

    yaxis_widget = Select(
        title='Displayed quantity:',
        value='flux',
        options=['flux', 'flux differences'],
    )

    thinning_widget = Spinner(
        title='Data thinning:',
        low=1,
        high=100,
        step=1,
        value=1,
        format='1',
    )
    thinning_widget.js_on_change('value', thin_data(args={'full_data': lc_data, 'displayed_data': displayed_data}))

    ephem_rows = []
    for ephem in Ephemeris.objects.filter(tic__tess_id=tess_id):
        button = Button(
            label=f'{ephem.bjd0 or 0:.6f} + E x {ephem.period or 1:.6f}',
            button_type='light',
            width=194,
            height=35,
            stylesheets=[stylesheet]
        )
        button.js_on_click(apply_ephem(
            args={
                't0': ephem.bjd0,
                'period': ephem.period,
                't0_widget': t0_widget,
                'period_widget': period_widget
            }))
        tooltip = Tooltip(content=f'source: {ephem.source.model} {(ephem.source.version) or ""}', position='right')
        help_button = HelpButton(tooltip=tooltip)
        ephem_rows.append(row([button, help_button]))

    add_ephem_widget = Button(
        label='Add new ephemeris',
        button_type='default',
        width=250,
        height=35,
    )

    pmin_widget = Spinner(
        title='Minimum period [days]:',
        low=0.1,
        high=1000,
        step=0.1,
        value=0.1,
        format='0.0',
    )

    pmax_widget = Spinner(
        title='Maximum period [days]:',
        low=0.1,
        high=1000,
        step=0.1,
        value=1,
        format='0.0',
    )

    pstep_widget = Spinner(
        title='Step size [days]:',
        low=0.0007,
        high=1.0,
        step=0.0007,
        value=0.0014,
        format='0.0000',
    )

    run_ls_widget = Button(
        label='Run Lomb-Scargle periodogram',
        button_type='default',
        width=250,
        height=35,
    )

    run_ls_widget.js_on_click(
        CustomJS(
            args=dict(source=displayed_data, spd=spd_data, fig=lc, pmin=pmin_widget, pmax=pmax_widget, pstep=pstep_widget, quantity=yaxis_widget),
            code="""
            // disable the button while the calculation is running:
            cb_obj.disabled = true;

            // allow switching between fluxes and flux differences:
            const qty = quantity.value === 'flux' ? 'flux' : 'dflux';

            // figure out which data to use:
            const xmin = fig.x_range.start;
            const xmax = fig.x_range.end;
            const x = source.data.time;
            const indices = [];
            for (let i = 0; i < x.length; i++) {
                if (x[i] >= xmin && x[i] <= xmax) {
                    indices.push(i);
                }
            }
            const filtered_time = Array.from(indices.map(i => source.data['time'][i]));
            const filtered_flux = Array.from(indices.map(i => source.data[qty][i]));
            console.log('length of the filtered data:', filtered_time.length);

            fetch('/api/lombscargle', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'Accept': 'application/json'
                },
                body: JSON.stringify({
                    'time': filtered_time,
                    'flux': filtered_flux,
                    'pmin': pmin.value,
                    'pmax': pmax.value,
                    'pstep': pstep.value,
                    'npeaks': 0,
                    'test': 'hello from js!'
                }),
            })
            .then(response => {
                if (!response.ok) {
                    console.log(response)
                    throw new Error('Failed to parse API response.');
                }
                return response.json();
            })
            .then(data => {
                if (data.period && data.power) {
                    spd.data = data;
                    spd.change.emit();
                } else {
                    console.log('Failed to assign data to the periodogram.');
                }
            })
            .catch(error => {
                console.error('Error:', error);
            })
            .finally(() => {
                cb_obj.disabled = false;
            });
            """
            )
    )

    ephem_model_widget = TextInput(
        title='Ephemeris model:',
        value='manual',
    )

    ephem_model_version_widget = TextInput(
        title='Ephemeris model version:',
        value='n/a',
    )

    ephem_model_reference_widget = TextInput(
        title='Ephemeris model reference:',
        value='none',
    )

    ephem_attribution_widget = TextInput(
        title='Ephemeris attribution',
        value='self',
    )

    commit_button = Button(label='Commit', button_type='success', width=250, height=35)
    cancel_button = Button(label='Cancel', button_type='danger', width=250, height=35)

    dialog_content = column(
        children=[
            Div(
                text='PLEASE VERIFY CAREFULLY BEFORE COMMITTING!',
                width=510,
                height=35,
            ),
            row([t0_widget, period_widget]),
            row([ephem_model_widget, ephem_model_version_widget]),
            row([ephem_model_reference_widget, ephem_attribution_widget]),
            row([commit_button, cancel_button]),
        ]
    )

    dialog = Dialog(
        title='Add ephemeris to the database',
        collapsible=False,
        maximizable=False,
        minimizable=False,
        pinnable=False,
        closable=True,
        close_action='hide',
        visible=False,
        syncable=False,
        content=dialog_content,
    )

    commit_button.js_on_click(CustomJS(
        args=dict(
            dialog=dialog,
            t0_widget=t0_widget,
            period_widget=period_widget,
            ephem_model_widget=ephem_model_widget,
            ephem_model_version_widget=ephem_model_version_widget,
            ephem_model_reference_widget=ephem_model_reference_widget,
            ephem_attribution_widget=ephem_attribution_widget
        ),
        code="""
            console.log('hiding the dialog');

            fetch('/api/ephem/add', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'Accept': 'application/json'
                },
                body: JSON.stringify({
                    't0': t0_widget.value,
                    'period': period_widget.value,
                    'model': ephem_model_widget.value,
                    'version': ephem_model_version_widget.value,
                    'reference': ephem_model_reference_widget.value,
                    'attribution': ephem_attribution_widget.value,
                }),
            })
            .then(response => {
                if (!response.ok) {
                    console.log(response)
                    throw new Error('Failed to parse API response.');
                }
                return response.json();
            })
            .then(data => {
                if (data.status === 'success') {
                    console.log('Successfully added ephemeris to the database.');
                } else {
                    console.log('Failed to add ephemeris to the database.');
                }
            })
            .catch(error => {
                console.error('Error:', error);
            })
            .finally(() => {
                dialog.visible = false;
            });
        """
    ))

    cancel_button.js_on_click(CustomJS(
        args=dict(dialog=dialog),
        code="""
            console.log('hiding the dialog');
            dialog.visible = false;
            dialog.change.emit();
        """
    ))

    add_ephem_widget.js_on_click(CustomJS(
        args=dict(dialog=dialog),
        code="""
            console.log('showing the dialog');
            dialog.visible = true;
            dialog.change.emit();
        """
    ))

    spd.js_on_event(events.Tap, CustomJS(
        args=dict(period=period_widget),
        code="""
            period.value = cb_obj.x;
        """
    ))

    controls = column(
        children=[
            yaxis_widget,
            thinning_widget,
            t0_widget,
            period_widget,
            harmonic_row,
            pmin_widget,
            pmax_widget,
            pstep_widget,
            run_ls_widget,
            *ephem_rows,
            add_ephem_widget,
        ],
        sizing_mode='stretch_width',
    )

    # Plot data:
    lc_plot = lc.circle('time', 'flux', size=2, color='blue', nonselection_alpha=0.1, source=displayed_data)
    sp_plot = spd.line('period', 'power', color='blue', source=spd_data)
    ph_plot = ph.circle('phase', 'flux', size=2, color='blue', nonselection_alpha=0.1, source=displayed_data)

    # Attach callbacks:
    yaxis_widget.js_on_change('value', switch_yaxis(args={'plots': [lc_plot, ph_plot], 'source': displayed_data}))
    t0_widget.js_on_change(
        'value',
        rephase(args={'plot': ph_plot, 'sources': [lc_data, displayed_data], 't0_widget': t0_widget, 'period_widget': period_widget})
    )
    period_widget.js_on_change(
        'value',
        rephase(args={'plot': ph_plot, 'sources': [lc_data, displayed_data], 't0_widget': t0_widget, 'period_widget': period_widget})
    )

    # UI layout:
    layout = row(
        children=[
            column(
                children=[controls],
                sizing_mode='fixed',
                width=300
            ),
            column(
                children=[lc, spd, ph],
                sizing_mode='stretch_width',
            ),
            # column(
            #     children=[logger],
            #     sizing_mode='fixed',
            #     width=300
            # )
            dialog,
        ],
        sizing_mode='stretch_width',
    )

    return components(layout, theme=theme)
