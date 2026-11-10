import datetime
import os
import random
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyqtgraph as pg
from PyQt5 import QtCore, QtWidgets
from scipy.optimize import curve_fit

from ccpn.api import PluginBase, PluginGUIModule

from ccpn.framework.Application import getApplication, getMainWindow, getCurrent, getProject
from ccpn.framework.plugins._GuiBaseFrame import FrameBase
from ccpn._old_plugins.development.pluginAddons import _addRow, _addColumn

from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.LineEdit import LineEdit
from ccpn.ui.gui.widgets.DoubleSpinbox import DoubleSpinbox
from ccpn.ui.gui.widgets.Button import Button
from ccpn.ui.gui.widgets.ButtonList import ButtonList
from ccpn.ui.gui.widgets.MessageDialog import showOkCancel
from ccpn.ui.gui.widgets.CheckBox import CheckBox
from ccpn.ui.gui.widgets.PulldownList import PulldownList
from ccpn.ui.gui.widgets.Spacer import Spacer
from ccpn.ui.gui.widgets.ScrollArea import ScrollArea
from ccpn.ui.gui.widgets.Frame import Frame

from ccpn.AnalysisMetabolomics.lib.NefDatabaseCaller import NefDatabaseCaller
from ccpn.AnalysisMetabolomics.lib.SpinSystemMatrixCreator import peakListFromSpinSystemMatrix
from ccpn.AnalysisMetabolomics.lib.LineshapeCreator import createLineshape
from ccpn.ui.gui.guiSettings import getColours, CCPNGLWIDGET_HEXBACKGROUND

# Global Variables
SettingsWidgetFixedWidths = (200, 350, 350)

# Widget tooltip texts
tooltips = {'SpectrumSelection'    : 'Select spectrum to fit to.',
        'ReferenceOffset'      : 'Chemical shift offset between the experimental spectra and the database standards.',
        'MultipletOffsetCap'   : 'Maximum random offset for multiplets when generating simulated database standards.',
        'SpectrometerFrequency': 'The frequency to simulate the reference standards.',
        'DatabasePath'         : 'The path to the database standards directory.',
        'PlotExportPath'       : 'Path to the file or directory for the matplotlib plot. Leave empty for no export.',
        'RunButton'            : 'Run the fit.'
        }

### Written for CCPN 3.4.0 by Alexander Hill and Morgan Hayward

### Developed for paper  A. D. Hill, G. Travish, M. Phelan, M. Hayward, and C. P. Welsch, “ Detection Limits of Blood Metabolites at Physiological Concentrations Using Benchtop 1
#H NMR,” NMR in Biomedicine 39, no. 2 (2026): e70215, https://doi.org/10.1002/nbm.70215. 

# Define the defualt database path here.
# In practice grab from https://github.com/ccpnmr/CASMDB
nef_database_path = '/Users/alexhill/Software/ccpnmr3.2.0/Extras/NEF_Database'

# Configuration for each metabolite used by the fitting workflow. Add a new
# metabolite here, rather than adding another compound-specific branch below.
# Explanation of spectrum_id name: the whole term (the spectrum ID) is used by the DatabaseCaller to find the specific data inside the file. The trimming of the spectrum ID is an easy way to find the file. While in theory, one can keep just the metabolite name and the sim number, it is safer to keep the whole ID to prevent unknown errors when trying to re-assemble these spectrum IDs. The prefix initially comes from the way Analysis automatically exports the spectrum data to NEF files. So it was quick and safe to keep the prefix, because it was more consistent than trying to remove it.


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


# List of hexcode colours for the plotting widget
tab10_hex_list = [c['color'] for c in plt.rcParams['axes.prop_cycle'] if 'color' in c]



class CurveFitGuiModule(PluginGUIModule):
    """A class to create the GUI element of the plugin as a GuiModule """
    def __init__(self, plugin, **kwargs):
        self.plugin = plugin
        self.mainWindow = getMainWindow()
        self.project = getProject()
        self.application = getApplication()
        self.current = getCurrent()
        super(CurveFitGuiModule, self).__init__(plugin, **kwargs)
        self.frameWidgets = FrameBase(parent=self.mainWidget, guiObject=self)
        self.widgetDefinitions = self.getWidgetDefinitions()

        # Set up the main gui area for the plugin. All widgets go in this GUI frame.
        self.scrollArea = ScrollArea(self.mainWidget, grid=(0, 0))
        self.scrollArea.setWidgetResizable(True)
        self.scrollAreaLayout = Frame(None, setLayout=True)
        self.scrollArea.setWidget(self.scrollAreaLayout)
        self.scrollAreaLayout.setContentsMargins(20, 20, 20, 20)

        self.signals = {}

        # Create the two frame for the parameter widgets (left) and plotting widget (right).
        parameterFrame = Frame(self.scrollAreaLayout, setLayout=True, grid=(0, 0))
        plotFrame = Frame(self.scrollAreaLayout, setLayout=True, grid=(0, 1))

        # Create parameter widgets. All parameter widgets here consist of a label and variable widget (spinbox, text edit etc.).
        # Grid for widget positions (column, row).
        grid = (0, 0)

        # Parameter widget for selecting the experimental spectrum.
        labelWidget = Label(parameterFrame, text='Spectra', grid=grid, valign='right')
        grid = _addColumn(grid)
        texts = ['> select-to-add <'] + [spectrum.pid for spectrum in self.project.spectra]
        self.spectrumSelectionWidget = PulldownList(parameterFrame, grid=grid, texts=texts, index=0,
                                                    tipText=tooltips['SpectrumSelection'])

        # Parameter widget for reference offset.
        grid = _addRow(grid)
        labelWidget = Label(parameterFrame, text='Reference Offset', grid=grid, valign='right')
        grid = _addColumn(grid)
        self.referenceOffsetWidget = DoubleSpinbox(parameterFrame, grid=grid, value=0.014, step=0.001, min=-10,
                                                   max=10, suffix='ppm', tipText=tooltips['ReferenceOffset'])

        # Parameter widget for multiplet offset cap.
        grid = _addRow(grid)
        labelWidget = Label(parameterFrame, text='Multiplet Offset Cap', grid=grid, valign='right')
        grid = _addColumn(grid)
        self.multipletOffsetCapWidget = DoubleSpinbox(parameterFrame, grid=grid, value=0, step=0.001, min=-0.1,
                                                      max=0.1, suffix='ppm', tipText=tooltips['MultipletOffsetCap'])

        # Parameter widget for frequency.
        grid = _addRow(grid)
        labelWidget = Label(parameterFrame, text='Spectrometer Frequency', grid=grid, valign='right')
        grid = _addColumn(grid)
        self.spectrometerFrequencyWidget = DoubleSpinbox(parameterFrame, grid=grid, value=80.0, step=10, min=10,
                                                         max=1200, suffix='MHz', tipText=tooltips['SpectrometerFrequency'])

        # Parameter widget for database path
        grid = _addRow(grid)
        labelWidget = Label(parameterFrame, text='Database Path', grid=grid, valign='right')
        grid = _addColumn(grid)
        self.databasePathWidget = LineEdit(parameterFrame, grid=grid,
                                           text=nef_database_path,
                                           tipText=tooltips['DatabasePath'])

        # Parameter widgets for metabolite selection. These are generated from
        # metabolite_config so the GUI stays in sync with the fit configuration.
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

        # Parameter widget for plot export path
        grid = _addRow(grid)
        labelWidget = Label(parameterFrame, text='Output Path', grid=grid, valign='right')
        grid = _addColumn(grid)
        self.plotExportPathWidget = LineEdit(parameterFrame, grid=grid, tipText=tooltips['PlotExportPath'])

        # The execute button.
        grid = _addRow(grid)
        self.runButton = Button(parameterFrame, grid=grid, text='Run', callback=self.runButtonCallback,
                                tipText=tooltips['RunButton'])

        # Add a spacer at the bottom to keep all the widgets together
        Spacer(parameterFrame, 5, 5, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding,
               grid=(grid[0] + 1, 10),
               gridSpan=(1, 1))

        # Add the plot widget to the right hand frame
        backgroundColour = getColours()[CCPNGLWIDGET_HEXBACKGROUND]
        self.plotWidget = pg.PlotWidget(background=backgroundColour)
        self.legend = self.plotWidget.addLegend(offset=(30, 30))
        self.plotWidget.invertX()
        plotFrame.getLayout().addWidget(self.plotWidget, 2, 0, 1, 1)

    def runButtonCallback(self):
        # The method called when the run button is clicked
        self.fit()
        spectraPid = self.spectrumSelectionWidget.getText()
        spectrum = self.project.getByPid(spectraPid)
        x = spectrum.positions
        y = spectrum.intensities
        self.plot(x, y)

    def fit(self):
        # Method for fitting the sim spectra to the sample spectrum.
        # Variables are collected from the widgets using _getValue() rather than being hard coded.
        # SSMs are put into a dict.

        # Get the caller to retrieve the substance sim data from the NEF database.
        self.caller = NefDatabaseCaller(project=self.project, directoryPaths=[self._getValue(self.databasePathWidget)])

        for substance, settings in metabolite_config.items():
            spectrumIds = settings['spectrum_ids']
            # Get the spectrum sim data from the NEF file using the configured IDs.
            specData = {key: self.caller.getFileData(substanceName=f'{"_".join(key.split("_")[3:-1])}') for key in spectrumIds}
            for specId in spectrumIds:
                fileData = specData[specId]
                multiplets, signals, ssm = self.caller.getSimulationData(fileData, specId)
                # Add any random offset if specified in the multipletOffsetCapWidget.
                for mIndex, multiplet in multiplets.iterrows():
                    multipletOffset = random.uniform(-self._getValue(self.multipletOffsetCapWidget),
                                                     self._getValue(self.multipletOffsetCapWidget))
                    for index in multiplet['ssm_indices']:
                        newShift = ssm[index][index] + self._getValue(self.referenceOffsetWidget) + multipletOffset
                        ssm[index][index] = newShift
                # Generate the spin system matrix.
                self.signals[substance] = peakListFromSpinSystemMatrix(spinSystemMatrix=ssm,
                                                                       frequency=self._getValue(self.spectrometerFrequencyWidget),
                                                                       width=1)

    def plot(self, x, y):
        # The method for plotting the results of the fit.
        # Exports the plot with matplotlib if the Output Path is given.
        # Also plots within the module plotWidget for inspection within the module.
        X = x
        Y = y

        arrays = []
        names = ['Sample', 'Full Sim']
        for name, settings in metabolite_config.items():
            widget = getattr(self, settings['widget'])
            if widget.isChecked():
                fitted_array, _, _ = self.get_fit(X, Y, substance=name)
                arrays.append(fitted_array)
                names.append(name)

        # Keep the plot useful even if no reference compounds were selected.
        y_total = np.sum(arrays, axis=0) if arrays else np.zeros_like(Y)

        selected_substances = [name for name in names[2:]]
        plot_title = 'Metabolite fit'
        if selected_substances:
            plot_title += ' — ' + ', '.join(selected_substances)

        exportPathValue = self._getValue(self.plotExportPathWidget)
        exportPath = Path(exportPathValue) if exportPathValue else None
        isFile = exportPath is not None and os.path.isfile(exportPath)
        isDir = exportPath is not None and os.path.isdir(exportPath)
        if isFile or isDir:
            fig, axs = plt.subplots(1, 1, figsize=(11, 5.5), constrained_layout=True)
            try:
                # Use the same ppm window as the analysis, but calculate the intensity
                # range from the data so prominent resonances are not clipped.
                finite_intensities = np.concatenate((np.asarray(Y), np.asarray(y_total)))
                finite_intensities = finite_intensities[np.isfinite(finite_intensities)]
                if finite_intensities.size:
                    intensity_min = finite_intensities.min()
                    intensity_max = finite_intensities.max()
                    intensity_range = intensity_max - intensity_min
                    padding = max(intensity_range * 0.05, abs(intensity_max) * 0.02, 1.0)
                    axs.set_ylim(intensity_min - padding, intensity_max + padding)

                axs.plot(X, Y, color='#1f77b4', lw=1.1, alpha=0.9,
                         label='Experimental')
                axs.plot(X, y_total, color='#d627c2', lw=1.4, ls='--',
                         label='Fitted simulation')
                axs.set_xlim(-1, 5.4)
                axs.invert_xaxis()
                axs.set_title(plot_title, fontsize=14, fontweight='bold', pad=12)
                axs.set_xlabel('Chemical shift (ppm)', labelpad=8)
                axs.set_ylabel('Intensity (a.u.)', labelpad=8)
                axs.ticklabel_format(axis='y', style='sci', scilimits=(0, 0),
                                    useMathText=True)
                axs.grid(axis='y', color='#b0b0b0', alpha=0.25, linewidth=0.8)
                axs.spines['top'].set_visible(False)
                axs.spines['right'].set_visible(False)
                axs.legend(loc='upper right', frameon=True, framealpha=0.95,
                           edgecolor='#d0d0d0')

                if isFile:
                    dialog = showOkCancel('File Exists', 'File already exists. Do you want to replace?')
                    if dialog:
                        fig.savefig(exportPath, dpi=300, bbox_inches='tight', facecolor='white')
                if isDir:
                    fileName = f'fitting_plot_{datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}.png'
                    fig.savefig(exportPath.joinpath(fileName), dpi=300,
                                bbox_inches='tight', facecolor='white')
            finally:
                plt.close(fig)

        # Plot in the plot widget
        # Get all the arrays we want to plot.
        all_arrays = [y] + [y_total] + arrays
        # Clear the plot widget for a new batch of arrays.
        self.plotWidget.clear()
        self.plotWidget.setTitle(plot_title)
        self.plotWidget.setLabel('bottom', 'Chemical shift', units='ppm')
        self.plotWidget.setLabel('left', 'Intensity', units='a.u.')
        # For each array, get a unique colour, based on matplotlib colours, and plot.
        for index, array in enumerate(all_arrays):
            colour = tab10_hex_list[index]
            pen = pg.functions.mkPen(color=colour, width=2, style=QtCore.Qt.SolidLine)
            line = self.plotWidget.plot(x=x,
                                        y=array,
                                        name=names[index],
                                        pen=pen)

    def get_fit(self, X, Y, substance='Glucose', plot=False):
        try:
            settings = metabolite_config[substance]
        except KeyError as error:
            raise ValueError(f'No fitting configuration found for {substance!r}') from error

        ppm_lower, ppm_upper = settings['ppm_bounds']
        mask = (X >= min(ppm_lower, ppm_upper)) & (X <= max(ppm_lower, ppm_upper))

        def fit_function(x_array=np.linspace(12, -1, 65536), scale=1.,
                         x_shift=0.0, width=2.0):
            return self.GenData(x_array, scale, x_shift=x_shift, width=width,
                                substance=substance)

        X_RED = X[mask]
        Y_RED = Y[mask]

        # Initial guesses
        initial_guesses = [
            np.nanmax(Y_RED),  # scale
            0,  # x_shift
            2  # linewidth
        ]

        # Perform curve fitting
        popt, pcov = curve_fit(fit_function, X_RED, Y_RED, p0=initial_guesses,
                               bounds=([0, -np.inf, -np.inf], np.inf))

        # Get optimised parameters
        opt_scale, opt_x_shift, opt_width = popt

        Y_fit = fit_function(x_array=X, scale=opt_scale, x_shift=opt_x_shift,
                             width=opt_width)

        return Y_fit, popt, pcov

    def GenData(self, x_array, scale, x_shift=0, substance='Lactate', width=2.):
        # Method for generating the simulated compound spectrum.
        # Uses the signals generated from the fit method.
        points = len(x_array)
        limits = (np.nanmax(x_array), np.nanmin(x_array))
        frequency = self._getValue(self.spectrometerFrequencyWidget)
        signals = self.signals[substance]
        signalList = [(row['chemical_shift'] + x_shift, row['height'], width / frequency) for index, row in
                      signals.iterrows()]
        x, y = createLineshape(signalList, points=points, limits=limits)
        y = y * scale
        return y

    def _getValue(self, widget):
        # Get the current value of the widget:
        if isinstance(widget, str):
            value = widget
        if not isinstance(widget, Button) and not isinstance(widget, ButtonList):
            if hasattr(widget, 'get'):
                value = widget.get()
            elif hasattr(widget, 'isChecked'):
                value = widget.isChecked()
            elif hasattr(widget, 'value'):
                value = widget.value()
        return value


class CurveFitModule(PluginBase):
    def __init__(self,  descriptor, application):
        """
             Plugin initialisation. Called automatically by the PluginManager.

             Framework-provided arguments:
                   :param descriptor: Metadata object describing the plugin (name, version, entry point, etc.).
                   :param application: Reference to the main application instance, used to access
                                       shared services such as project data, managers, and configuration.
             Notes:
               - These parameters are injected by the plugin manager and are not intended to be supplied by end users.
        """
        super().__init__(descriptor, application)
        self.ui = CurveFitGuiModule # we need to attach the Gui element here. See class above

    def run(self, *args, **kwargs):
        """
        Default entry point for executing a plugin action.
        This method is called automatically when the user clicks *Run* in the plugin’s GUI if there is a GUI element.
        By default this method does nothing. Subclasses should override it
        to implement the plugin’s behaviour. Typical uses include:

          • Running an external command or script via the provided process runner
          • Attaching a file watcher to monitor results and feed them back into the project
          • Performing a direct action on the current project (e.g. adding, updating, or analysing data)

        The keyword arguments (`kwargs`) are automatically populated from
        the plugin’s GUI layer, based on the widgets defined in
        `getWidgetDefinitions()`. Each widget contributes its value under
        the symbolic name you assigned there.
        """
        sys.stdout.write(f"\n Run method called with arguments: {kwargs}")
