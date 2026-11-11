# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Alexander Hill

"""MetabFit: fit simulated NEF references to an experimental spectrum in CCPN.

Distributed under the MIT licence; see LICENSE in the repository root.

Written for CCPN 3.4.0 by Alexander Hill and Morgan Hayward for:
A. D. Hill et al., NMR in Biomedicine 39(2), e70215 (2026).
https://doi.org/10.1002/nbm.70215
"""

import csv
import datetime
import random
import sys
from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PyQt5 import QtCore, QtWidgets
from pyqtgraph.exporters import ImageExporter
from scipy.optimize import curve_fit

from ccpn.api import PluginBase, PluginGUIModule

from ccpn.framework.Application import getApplication, getMainWindow, getCurrent, getProject
from ccpn.framework.plugins._GuiBaseFrame import FrameBase
from ccpn._old_plugins.development.pluginAddons import _addRow, _addColumn

from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.LineEdit import LineEdit
from ccpn.ui.gui.widgets.DoubleSpinbox import DoubleSpinbox
from ccpn.ui.gui.widgets.Button import Button
from ccpn.ui.gui.widgets.MessageDialog import showOkCancel, showWarning
from ccpn.ui.gui.widgets.CheckBox import CheckBox
from ccpn.ui.gui.widgets.PulldownList import PulldownList
from ccpn.ui.gui.widgets.Spacer import Spacer
from ccpn.ui.gui.widgets.ScrollArea import ScrollArea
from ccpn.ui.gui.widgets.Frame import Frame

from ccpn.AnalysisMetabolomics.lib.NefDatabaseCaller import NefDatabaseCaller
from ccpn.AnalysisMetabolomics.lib.SpinSystemMatrixCreator import peakListFromSpinSystemMatrix
from ccpn.AnalysisMetabolomics.lib.LineshapeCreator import createLineshape
from ccpn.ui.gui.guiSettings import getColours, CCPNGLWIDGET_HEXBACKGROUND

# Tooltips describe inputs; fitting and export defaults are set in the GUI below.
tooltips = {
    'SpectrumSelection': 'Select spectrum to fit to.',
    'ReferenceOffset': 'Chemical shift offset between the experimental spectra and the database standards.',
    'MultipletOffsetCap': 'Maximum random offset for multiplets when generating simulated database standards.',
    'SpectrometerFrequency': 'The frequency to simulate the reference standards.',
    'DatabasePath': 'The path to the database standards directory.',
    'PlotExportPath': 'Existing output directory or filename. PNG and CSV exports share the same base name.',
    'RunButton': 'Run the fit.',
}

# Default reference database location; editable in the GUI.
# Download database standards from https://github.com/ccpnmr/CASMDB.
nef_database_path = '/Users/Data/NEF_Database'

# Single source for metabolite checkboxes, NEF references, and fitting windows.
# Keep the complete spectrum ID for record lookup; filename lookup strips the
# first three underscore-separated fields and the final simulation identifier.
# Use one spectrum ID per metabolite: additional IDs overwrite its peak list.

metabolite_config = {
    'Glucose': {
        'spectrum_ids': ['nef_nmr_spectrum_alpha-d-glucose_6`1`'],
        'ppm_bounds': (3.19, 3.98),
        'widget': 'glucoseWidget',
    },
    'Citrate': {
        'spectrum_ids': ['nef_nmr_spectrum_citric_acid_8`1`'],
        'ppm_bounds': (2.37, 2.82),
        'widget': 'citrateWidget',
    },
    'Lactate': {
        'spectrum_ids': ['nef_nmr_spectrum_d-lactic_acid_3`1`'],
        'ppm_bounds': (1.17, 1.50),
        'widget': 'lactateWidget',
    },
}


# Fixed tab10 colours keep GUI and PNG traces consistent without Matplotlib.
PLOT_COLOURS = (
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
)


class CurveFitGuiModule(PluginGUIModule):
    """Configure reference simulations, fit selected metabolites, and plot results."""

    def __init__(self, plugin, **kwargs):
        self.plugin = plugin
        self.mainWindow = getMainWindow()
        self.project = getProject()
        self.application = getApplication()
        self.current = getCurrent()
        super(CurveFitGuiModule, self).__init__(plugin, **kwargs)
        # Brand the panel label without changing the plugin's registration name
        # or the Python identifiers used by CCPN to load it.
        self.label.setText('MetabFit')
        self.frameWidgets = FrameBase(parent=self.mainWidget, guiObject=self)
        self.widgetDefinitions = self.getWidgetDefinitions()

        # Scrolling keeps the controls accessible in smaller CCPN module panels.
        self.scrollArea = ScrollArea(self.mainWidget, grid=(0, 0))
        self.scrollArea.setWidgetResizable(True)
        self.scrollAreaLayout = Frame(None, setLayout=True)
        self.scrollArea.setWidget(self.scrollAreaLayout)
        self.scrollAreaLayout.setContentsMargins(20, 20, 20, 20)

        self.signals = {}
        self._plotData = None
        self._exportStem = None

        # Place parameters on the left and the spectrum plot on the right.
        parameterFrame = Frame(self.scrollAreaLayout, setLayout=True, grid=(0, 0))
        plotFrame = Frame(self.scrollAreaLayout, setLayout=True, grid=(0, 1))

        # Each input has a label; CCPN grid positions are (row, column).
        grid = (0, 0)

        Label(parameterFrame, text='Spectra', grid=grid, valign='right')
        grid = _addColumn(grid)
        texts = ['> select-to-add <'] + [spectrum.pid for spectrum in self.project.spectra]
        self.spectrumSelectionWidget = PulldownList(parameterFrame, grid=grid, texts=texts, index=0,
                                                    tipText=tooltips['SpectrumSelection'])

        grid = _addRow(grid)
        Label(parameterFrame, text='Reference Offset', grid=grid, valign='right')
        grid = _addColumn(grid)
        self.referenceOffsetWidget = DoubleSpinbox(parameterFrame, grid=grid, value=0.014, step=0.001, min=-10,
                                                   max=10, suffix='ppm', tipText=tooltips['ReferenceOffset'])

        grid = _addRow(grid)
        Label(parameterFrame, text='Multiplet Offset Cap', grid=grid, valign='right')
        grid = _addColumn(grid)
        self.multipletOffsetCapWidget = DoubleSpinbox(parameterFrame, grid=grid, value=0, step=0.001, min=-0.1,
                                                      max=0.1, suffix='ppm', tipText=tooltips['MultipletOffsetCap'])

        grid = _addRow(grid)
        Label(parameterFrame, text='Spectrometer Frequency', grid=grid, valign='right')
        grid = _addColumn(grid)
        self.spectrometerFrequencyWidget = DoubleSpinbox(parameterFrame, grid=grid, value=80.0, step=10, min=10,
                                                         max=1200, suffix='MHz', tipText=tooltips['SpectrometerFrequency'])

        grid = _addRow(grid)
        Label(parameterFrame, text='Database Path', grid=grid, valign='right')
        grid = _addColumn(grid)
        self.databasePathWidget = LineEdit(parameterFrame, grid=grid,
                                           text=nef_database_path,
                                           tipText=tooltips['DatabasePath'])

        # Generate selections from the same configuration used by the fitter.
        for substance, settings in metabolite_config.items():
            grid = _addRow(grid)
            Label(parameterFrame, text=substance, grid=grid, valign='right')
            grid = _addColumn(grid)
            setattr(self, settings['widget'], CheckBox(
                parameterFrame,
                grid=grid,
                checked=False,
                tipText=f'Simulate {substance}',
            ))

        grid = _addRow(grid)
        Label(parameterFrame, text='Output Path', grid=grid, valign='right')
        grid = _addColumn(grid)
        self.plotExportPathWidget = LineEdit(parameterFrame, grid=grid, tipText=tooltips['PlotExportPath'])

        grid = _addRow(grid)
        self.runButton = Button(parameterFrame, grid=grid, text='Run', callback=self.runButtonCallback,
                                tipText=tooltips['RunButton'])

        grid = _addRow(grid)
        self.saveFigButton = Button(parameterFrame, grid=grid, text='Save Fig',
                                    callback=self.saveFigure)
        grid = _addColumn(grid)
        self.saveCsvButton = Button(parameterFrame, grid=grid, text='Save CSV',
                                    callback=self.saveCsv)
        self.saveFigButton.setEnabled(False)
        self.saveCsvButton.setEnabled(False)

        # Keep controls together when the module grows vertically.
        Spacer(parameterFrame, 5, 5, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding,
               grid=(grid[0] + 1, 10),
               gridSpan=(1, 1))

        # Match CCPN's background and the descending NMR chemical-shift convention.
        backgroundColour = getColours()[CCPNGLWIDGET_HEXBACKGROUND]
        self.plotWidget = pg.PlotWidget(background=backgroundColour)
        self.legend = self.plotWidget.addLegend(offset=(30, 30))
        self.plotWidget.invertX()
        self._createViewControls(plotFrame)
        plotFrame.getLayout().addWidget(self.plotWidget, 2, 0, 1, 1)

    def _createViewControls(self, parent):
        """Provide exact axis limits alongside mouse zoom and pan controls."""
        controls = Frame(parent, setLayout=True, grid=(0, 0))
        self._axisFields = []
        for column, (label, initial) in enumerate((
                ('PPM min', '0'), ('PPM max', '1'),
                ('Intensity min', '0'), ('Intensity max', '1'))):
            Label(controls, text=label, grid=(0, column))
            field = LineEdit(controls, text=initial, grid=(1, column))
            field.setToolTip('Enter an axis limit; scientific notation is accepted. Click Apply View.')
            field.returnPressed.connect(self.applyView)
            self._axisFields.append(field)

        actions = Frame(parent, setLayout=True, grid=(1, 0))
        Button(actions, text='Reset', grid=(0, 0), callback=self.autoScaleView)
        Button(actions, text='Apply View', grid=(0, 1), callback=self.applyView)
        Button(actions, text='Zoom Region', grid=(0, 2), callback=self.zoomRegion)
        Button(actions, text='Pan', grid=(0, 3), callback=self.panView)

        # Keep numeric limits aligned with mouse gestures and automatic scaling.
        view = self.plotWidget.getViewBox()
        view.sigRangeChanged.connect(self._syncViewControls)
        self.panView()
        self._syncViewControls()

    def _syncViewControls(self, *args):
        """Show numeric bounds in ascending order, even with the NMR axis inverted."""
        x_range, y_range = self.plotWidget.getViewBox().viewRange()
        for field, value in zip(self._axisFields, (*x_range, *y_range)):
            field.setText(f'{value:.12g}')

    def applyView(self):
        """Apply finite, increasing bounds without PyQtGraph's default padding."""
        try:
            limits = [float(field.text()) for field in self._axisFields]
        except ValueError:
            showWarning('Invalid axis limits', 'Enter a number in each axis field.')
            return
        x_min, x_max, y_min, y_max = limits
        if not all(np.isfinite(limits)) or x_min >= x_max or y_min >= y_max:
            showWarning('Invalid axis limits',
                        'Each minimum must be smaller than its maximum, and all limits must be finite.')
            return
        self.plotWidget.getViewBox().setRange(
            xRange=(x_min, x_max), yRange=(y_min, y_max),
            padding=0, disableAutoRange=True,
        )

    def autoScaleView(self):
        """Handle Reset by fitting all traces once, preserving later manual zooms."""
        self.plotWidget.getViewBox().autoRange()

    def zoomRegion(self):
        """Enable left-drag rectangle zoom over a region of interest."""
        self.plotWidget.getViewBox().setMouseMode(pg.ViewBox.RectMode)

    def panView(self):
        """Enable left-drag panning; the mouse wheel can still zoom."""
        self.plotWidget.getViewBox().setMouseMode(pg.ViewBox.PanMode)

    def runButtonCallback(self):
        """Prepare reference signals and plot fits for the selected spectrum."""
        self.fit()
        spectraPid = self.spectrumSelectionWidget.getText()
        spectrum = self.project.getByPid(spectraPid)
        x = spectrum.positions
        y = spectrum.intensities
        self.plot(x, y)

    def fit(self):
        """Load NEF spin systems and prepare reference peak lists for fitting.

        Each multiplet receives one shared random shift, preserving its internal
        splitting. The reference offset applies to every diagonal chemical shift.
        All configured references are loaded, including unchecked metabolites;
        parameter optimisation happens later in get_fit().
        """

        # Database paths are supplied by the user in the CCPN environment.
        self.caller = NefDatabaseCaller(project=self.project, directoryPaths=[self._getValue(self.databasePathWidget)])

        for substance, settings in metabolite_config.items():
            spectrumIds = settings['spectrum_ids']
            # Derive the substance filename key, but retain the full NEF record ID.
            specData = {
                key: self.caller.getFileData(
                    substanceName='_'.join(key.split('_')[3:-1])
                )
                for key in spectrumIds
            }
            for specId in spectrumIds:
                fileData = specData[specId]
                multiplets, _, ssm = self.caller.getSimulationData(fileData, specId)
                # Shift diagonal chemical shifts together within each multiplet;
                # off-diagonal coupling values remain unchanged.
                for _, multiplet in multiplets.iterrows():
                    multipletOffset = random.uniform(-self._getValue(self.multipletOffsetCapWidget),
                                                     self._getValue(self.multipletOffsetCapWidget))
                    for index in multiplet['ssm_indices']:
                        newShift = ssm[index][index] + self._getValue(self.referenceOffsetWidget) + multipletOffset
                        ssm[index][index] = newShift
                # Convert the shifted spin system into reference peak signals.
                self.signals[substance] = peakListFromSpinSystemMatrix(spinSystemMatrix=ssm,
                                                                       frequency=self._getValue(self.spectrometerFrequencyWidget),
                                                                       width=1)

    def plot(self, x, y):
        """Fit selected metabolites and retain the displayed data for export."""
        X = x
        Y = y

        arrays = []
        names = ['Experimental', 'Total Fit']
        for name, settings in metabolite_config.items():
            widget = getattr(self, settings['widget'])
            if widget.isChecked():
                fitted_array, _, _ = self.get_fit(X, Y, substance=name)
                arrays.append(fitted_array)
                names.append(name)

        # Keep the plot useful even if no reference compounds were selected.
        y_total = np.sum(arrays, axis=0) if arrays else np.zeros_like(Y)

        selected_substances = names[2:]
        plot_title = 'Metabolite fit'
        if selected_substances:
            plot_title += ' — ' + ', '.join(selected_substances)

        # Plot the experimental data, sum, and components in matching legend order.
        all_arrays = [y, y_total] + arrays
        self.plotWidget.clear()
        self.plotWidget.setTitle(plot_title)
        self.plotWidget.setLabel('bottom', 'Chemical shift', units='ppm')
        self.plotWidget.setLabel('left', 'Intensity', units='a.u.')
        for index, array in enumerate(all_arrays):
            colour = PLOT_COLOURS[index % len(PLOT_COLOURS)]
            pen = pg.functions.mkPen(color=colour, width=2, style=QtCore.Qt.SolidLine)
            self.plotWidget.plot(x=x, y=array, name=names[index], pen=pen)

        # Snapshot the fitted arrays: exporting must not rerun the fit or depend
        # on later changes to checkboxes, spectrum selection, or fitting inputs.
        self._plotData = (
            np.array(x, copy=True),
            list(names),
            [np.array(array, copy=True) for array in all_arrays],
        )
        self._exportStem = 'fitting_plot_' + datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S_%f')
        self.saveFigButton.setEnabled(True)
        self.saveCsvButton.setEnabled(True)
        self.autoScaleView()

    def _exportPath(self, suffix):
        """Resolve a directory or filename to a PNG/CSV path; confirm replacement.

        Directories use the current fit's timestamp. A filename's suffix is
        replaced by the requested format; parent directories must already exist.
        """
        if self._plotData is None:
            showWarning('No fit available', 'Run the fit before saving.')
            return None

        value = self._getValue(self.plotExportPathWidget).strip()
        if not value:
            showWarning('Output Path required', 'Enter an output directory or filename.')
            return None

        path = Path(value).expanduser()
        if path.is_dir():
            path = path / (self._exportStem + suffix)
        else:
            path = path.with_suffix(suffix)
        if not path.parent.is_dir():
            showWarning('Invalid Output Path', 'The output directory must already exist.')
            return None
        if path.exists() and not showOkCancel(
                'File Exists', f'{path}\nFile already exists. Do you want to replace it?'):
            return None
        return path

    def _figureImage(self):
        """Render the live PlotItem for saving, without refitting.

        Do not auto-scale here: the ViewBox's current limits, including mouse
        zooms and explicit axis limits, define the region exported by PyQtGraph.
        """
        exporter = ImageExporter(self.plotWidget.getPlotItem())
        exporter.parameters()['width'] = 2200
        # LegendItem normally ignores transforms to keep a constant screen size.
        # During high-resolution export that leaves the legend at screen-pixel
        # size while the axes and plot grow. Temporarily let the entire legend
        # (text, samples, and spacing) scale with the scene instead.
        ignore_transforms = QtWidgets.QGraphicsItem.ItemIgnoresTransformations
        originally_ignored = bool(self.legend.flags() & ignore_transforms)
        try:
            self.legend.setFlag(ignore_transforms, False)
            return exporter.export(toBytes=True)
        finally:
            self.legend.setFlag(ignore_transforms, originally_ignored)

    def saveFigure(self):
        """Save the visible plot region with the same styling as the CCPN module."""
        path = self._exportPath('.png')
        if path is None:
            return
        try:
            if not self._figureImage().save(str(path), 'PNG'):
                raise OSError(f'Could not save {path}')
        except Exception as error:
            showWarning('Figure export failed', str(error))

    def saveCsv(self):
        """Save full-resolution experimental, component, and total intensities."""
        path = self._exportPath('.csv')
        if path is None:
            return
        x, names, arrays = self._plotData
        # Put the total last, after the separate metabolite fits. Keep the axis
        # in its original order and export all points regardless of GUI zoom.
        columns = [x, arrays[0], *arrays[2:], arrays[1]]
        header = ['Chemical shift (ppm)', names[0], *names[2:], names[1]]
        try:
            with path.open('w', newline='', encoding='utf-8') as output:
                writer = csv.writer(output)
                writer.writerow(header)
                writer.writerows(zip(*columns))
        except OSError as error:
            showWarning('CSV export failed', str(error))

    def get_fit(self, X, Y, substance='Glucose', plot=False):
        """Fit scale, shift (ppm), and linewidth (Hz) within a compound's window.

        Return (full-axis simulation, parameters, covariance), with parameters
        ordered as scale, shift, linewidth. Each compound uses the experimental
        data directly, without subtracting other fits. Only scale is bounded.
        The unused plot argument is retained for compatibility.
        """
        try:
            settings = metabolite_config[substance]
        except KeyError as error:
            raise ValueError(f'No fitting configuration found for {substance!r}') from error

        ppm_lower, ppm_upper = settings['ppm_bounds']
        mask = (X >= min(ppm_lower, ppm_upper)) & (X <= max(ppm_lower, ppm_upper))

        def fit_function(x_array, scale=1.,
                         x_shift=0.0, width=2.0):
            return self.GenData(x_array, scale, x_shift=x_shift, width=width,
                                substance=substance)

        X_RED = X[mask]
        Y_RED = Y[mask]

        # Start from the window's maximum intensity, zero shift, and 2 Hz width.
        # Shift and linewidth remain unbounded; inspect their physical plausibility.
        initial_guesses = [
            np.nanmax(Y_RED),  # Amplitude scale (not a calibrated concentration).
            0,                # Additional chemical shift in ppm.
            2,                # Linewidth in Hz.
        ]

        popt, pcov = curve_fit(fit_function, X_RED, Y_RED, p0=initial_guesses,
                               bounds=([0, -np.inf, -np.inf], np.inf))

        # Evaluate on the full axis so the fitted components can be summed.
        opt_scale, opt_x_shift, opt_width = popt

        Y_fit = fit_function(x_array=X, scale=opt_scale, x_shift=opt_x_shift,
                             width=opt_width)

        return Y_fit, popt, pcov

    def GenData(self, x_array, scale, x_shift=0, substance='Lactate', width=2.):
        """Build a scaled lineshape from the reference peaks prepared by fit().

        createLineshape reconstructs a uniform descending axis from the extrema
        and point count; it does not evaluate arbitrary supplied x coordinates.
        The input must therefore use that grid. Hz / MHz gives linewidth in ppm;
        scale multiplies the reference peak heights and is not a concentration.
        """
        points = len(x_array)
        limits = (np.nanmax(x_array), np.nanmin(x_array))
        frequency = self._getValue(self.spectrometerFrequencyWidget)
        signals = self.signals[substance]
        signalList = [
            (row['chemical_shift'] + x_shift, row['height'], width / frequency)
            for _, row in signals.iterrows()
        ]
        _, y = createLineshape(signalList, points=points, limits=limits)
        y = y * scale
        return y

    def _getValue(self, widget):
        """Read a value from a CCPN input widget or a literal string."""
        if isinstance(widget, str):
            return widget
        for accessor in ('get', 'isChecked', 'value'):
            method = getattr(widget, accessor, None)
            if callable(method):
                return method()
        raise TypeError(f'Unsupported input widget: {type(widget).__name__}')


class MetabFit(PluginBase):
    """Register MetabFit with CCPN using the existing MetabFit identity."""

    def __init__(self, descriptor, application):
        """Attach the GUI using descriptor and application supplied by CCPN."""
        super().__init__(descriptor, application)
        self.ui = CurveFitGuiModule

    def run(self, *args, **kwargs):
        """Log framework calls; the GUI Run button uses runButtonCallback()."""
        sys.stdout.write(f"\n Run method called with arguments: {kwargs}")
