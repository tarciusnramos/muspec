'''
Add:
    xshift
    deal with the legend
    create logs
'''

'''
ERROR SAVING SVG WHEN MULTIPLE FILES ARE LOADED.
FIX IT: 
include the command: s = str(s) 
in the beginning of the function escape_attrib(s)
at Anaconda3/lib/site-packages/matplotlib/backends/backend_svg.py
'''

from PyQt5 import QtCore, QtGui, QtWidgets
import sys, os, logging
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import numpy as np
import Gaussian, Dalton
from tools import *

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

#CONSTANTS
EV2NM = 1240
EV2CM = 8066

INSTALL_DIR = os.path.dirname(__file__) #+ '/'
ICON_DIR = os.path.join(INSTALL_DIR, 'images/')

class XYItem(object):
    """docstring for XYItem"""
    def __init__(self, file):
        super(XYItem, self).__init__()
        self.label = os.path.basename(file)[:50]
        self.read_data(file)
        self.ngroups = 1
        self.xcolumns = np.array([0])
        self.ycolumns = np.array([1])
        self.xyColumnsLabelString = '0,1'
        self.xmult = 1.0
        self.ymult = 1.0
        self.plot_hides = False

    def read_data(self, file):
        comments = ('#', '!', '@')
        try:
            self.data = np.loadtxt(file, comments=comments).T
            self.create_status = True
            self.ncols = len(self.data)
        except:
            self.create_status = False        

    def update_item(self):
        w.validate_xy()

    def draw(self):
        if self.plot_hides:
            return

        self.get_lines()

        for n in range(self.ngroups):
            x = self.get_x(self.spectra_x[n])
            y = self.get_y(self.spectra_y[n])
            self.ax.plot(x, y, gid=self.item_id)
        w.canvas.draw()

    def get_lines(self):
        self.ngroups = len(self.ycolumns)
        self.spectra_x = []
        self.spectra_y = []
        for n in range(self.ngroups):
            spectra_x = self.data[self.xcolumns[n]]
            self.spectra_y += [self.data[self.ycolumns[n]]]
            
            if all(x < 50 for x in spectra_x):
                pass
            elif all(x > 1000 for x in spectra_x):
                logging.info('Converting energies from cm$^{-1}$ to eV')
                spectra_x /= EV2CM
            else:
                logging.info('Converting energies from nm to eV')
                spectra_x = EV2NM / spectra_x

            self.spectra_x += [spectra_x]
            
        self.spectra_x = np.asarray(self.spectra_x)
        self.spectra_y = np.asarray(self.spectra_y)

    def remove_draw(self):
        nmax = len(self.ax.lines)
        for n in range(nmax-1,-1,-1):
            if self.ax.lines[n].get_gid() == self.item_id:
                self.ax.lines[n].remove()

    def get_y(self, y):
        w.update_nomalization(self.item_id)
        return y * self.ymult / w.y_norm

    def get_x(self, x):
        return w.set_energy_unit(x * self.xmult)

    def remove(self):
        self.remove_draw()
        logging.info('Removing item {0}'.format(self.label))
        w.canvas.draw()

class VLItem(Tools):
    """docstring for VLItem"""
    def __init__(self, energies, intensities, label):
        super(VLItem, self).__init__()
        self.energies = energies
        self.intensities = intensities
        self.label = label[:50]
        self.initialize_defaults()
        self.energies_in_ev()
        
        self.initialize_item()
        logging.info('There are {0} states and {1} configurations.'.format(self.nstates, len(self.energies[1])))
        
    def initialize_defaults(self):
        self.lim = [2, 8]
        self.steps = 1000
        self.xmult = 1.0
        self.ymult = 1.0
        self.vline = False
        self.plot_group = False
        self.plot_color_group = False
        self.plot_hides = False
        self.function = 'Lorentzian'
        self.procedure = 'procedure_03'

    def energies_in_ev(self):
        list_energies = []
        for key in self.energies.keys():
            for n in self.energies[key]:
                list_energies.append(n)
        if all(x < 50 for x in list_energies):
            pass
        elif all(x > 1000 for x in list_energies):
            logging.info('Converting energies from cm$^{-1}$ to eV')
            self.energies.update({n: np.array(self.energies[n]) / EV2CM for n in self.energies.keys()})
        else:
            logging.info('Converting energies from nm to eV')
            self.energies.update({n: EV2NM / np.array(self.energies[n]) for n in self.energies.keys()})

    def initialize_item(self):
        states = list(self.energies.keys())
        self.nstates = len(states)
        s = w.groups_to_string([states])
        self.set_groups(s)
        self.fwhm = [ 0.4 for x in range(self.ngroups)]
        self.broad_func = self.broad_bib[self.function]

    def set_groups(self, arg):
        self.groups = []
        for gr in arg.split(';'):
            g_temp = []
            for gs in gr.split(','):
                g_ = gs.split('-')
                g_temp.extend([n + 1 for n in range(int(g_[0])-1, int(g_[-1]))])
            for n in range(len(g_temp)):
                if g_temp[n] > self.nstates:
                    g_temp[n] = self.nstates
            self.groups.append(g_temp)
        self.groups = np.asarray(self.groups)
        self.ngroups = len(self.groups)

    def set_fwhm(self, arg):
        g = arg.split(';')
        g_len = len(g)
        self.fwhm = []

        for ng in range(self.ngroups):
            if ng < g_len:
                self.fwhm += [ float(g[ng]) ]
            else:
                logging.warning('The number of groups and FWHM doesn\'t match. Setting the last parameter as 0.4 eV.')
                for ng in range(self.ngroups - g_len):
                    self.fwhm += [ 0.4 ]

    def update_par(self, *arg):
        for par in arg:
            if par[0] == 'group':
                self.set_groups(par[1])
            elif par[0] == 'fwhm':
                self.set_fwhm(par[1])

    def update_item(self):
        w.validate_groups()
        w.validate_fwhm()
        w.validate_ymult()
        w.validate_xrange()
        self.spectra()

    def spectra(self):
        self.spectra_x = np.linspace(self.lim[0], self.lim[1], self.steps)
        self.procedures_bib[self.procedure]()
        self.get_vlines()

    def get_vlines(self):
        self.vlines_x = []
        self.vlines_y = []
        for ng in self.groups:
            vx, vy = [], []
            for state in ng:
                vx += [*self.energies[state]]
                vy += [*self.intensities[state]]
            self.vlines_x += [np.array(vx)]
            self.vlines_y += [np.array(vy)]
        self.vlines_x = np.asarray(self.vlines_x)
        self.vlines_y = np.asarray(self.vlines_y)

    def draw(self):
        self.remove_draw()
        self.spectra()
        self.draw_spectra()
        self.draw_vlines()

    def remove_draw(self):
        self.remove_lines()
        self.remove_vlines()

    def remove_lines(self):
        nmax = len(self.ax.lines)
        for n in range(nmax-1,-1,-1):
            if self.ax.lines[n].get_gid() == self.item_id:
                self.ax.lines[n].remove()

    def remove_vlines(self):
        nmax = len(self.ax.collections)
        for n in range(nmax-1,-1,-1):
            if self.ax.collections[n].get_gid() == self.item_id:
                self.ax.collections[n].remove()

    def draw_spectra(self):
        if self.plot_hides:
            return
        if self.plot_group:
            linestyle = '-'
            label = ''
        else:
            linestyle = 'None'
            label = '_'
        
        x = self.get_x(self.spectra_x)

        self.colors = []
        text = self.label + '_total'
        y_tot = np.sum(self.get_y(self.spectra_y), axis=0)
        p = self.ax.plot(x, y_tot, gid=self.item_id, label=text, ls='-')
        self.colors += [p[0].get_color()]

        nmax = self.ngroups
        for n in range(nmax):
            y_tmp = self.get_y(self.spectra_y[n])
            if self.ngroups > 1:
                text = label + self.label + '_g' + str(n + 1)
                p = self.ax.plot(x, y_tmp, gid=self.item_id, label=text, ls=linestyle)
                self.colors += [p[0].get_color()]

    def draw_vlines(self):
        if self.plot_hides or not self.vline:
            return
        
        x = self.get_x(self.vlines_x)
        y = self.get_y(self.vlines_y)

        nmax = self.ngroups
        plot_group = self.plot_group or self.plot_color_group
        condition = (nmax > 1) and plot_group

        if condition:
            for n in range(nmax):
                self.ax.vlines(x[n], np.zeros(len(y[n])), y[n], gid=self.item_id, color=self.colors[n + 1], alpha=0.6)
        else:
            xc = np.concatenate(x[:])
            yc = np.concatenate(y[:])
            self.ax.vlines(xc, np.zeros(len(yc)), yc, gid=self.item_id, color=self.colors[0], alpha=0.6)

    def get_y(self, y):
        w.update_nomalization(self.item_id)
        return y * self.ymult / w.y_norm

    def get_x(self, x):
        return w.set_energy_unit(x * self.xmult)

    def remove(self):
        self.remove_lines()
        self.remove_vlines()
        logging.info('Removing item {0}'.format(self.label))
        w.canvas.draw()

class QTextEditLogger(logging.Handler):
    def __init__(self, parent):
        super().__init__()
        self.widget = QtWidgets.QPlainTextEdit(parent)
        self.widget.setReadOnly(True)

    def emit(self, record):
        msg = self.format(record)
        self.widget.appendPlainText(msg)

class QTextAnalyses(QtWidgets.QTextEdit):
    """docstring for QTextAnalyses"""
    def __init__(self, window_name):
        super(QTextAnalyses, self).__init__()
        self.setWindowTitle(window_name)
        self.setWindowIcon(QtGui.QIcon(ICON_DIR + 'muspec.png'))
        self.text = QtWidgets.QTextEdit(self)
        self.text.setMinimumSize(480, 320)
        self.text.setReadOnly(True)
        
        # Layout
        self.layout_analyses = QtWidgets.QGridLayout()
        self.layout_analyses.addWidget(self.text, 1, 0)
        self.layout_analyses.setContentsMargins(5, 5, 5, 5)
        self.layout_analyses.setSpacing(5)

        self.setLayout(self.layout_analyses)   

        self.adjustSize()
        

    def emit(self, message):
        self.text.clear()
        self.text.insertPlainText(message)
        
class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()

        self.setWindowTitle("MUSPEC")
        self.setWindowIcon(QtGui.QIcon(ICON_DIR + 'muspec.png'))
        self.last_dir = os.getcwd()

        self.create_grid()
        self.create_menu()
        self.setCentralWidget(self.frame)
        self.create_actions()
        self.create_toolbar()

        self.item_id = 0
        self.item = {}

        self.row = self.filesView.currentRow()
        self.row_id = None

        self.y_norm = 1.0

        self.fill = {'FileItem':self.fill_item,
                     'VLItem':self.fill_item,
                     'XYItem':self.fill_xy}
        self.infoview = self.infoViewGroupBox

        self.fill_normalization_cbox()

    def create_menu(self):
        mainMenu = self.menuBar()
        
        ## File menu
        fileMenu = mainMenu.addMenu('File')
        #Load QM button
        load_qm_button = QtWidgets.QAction(QtGui.QIcon(ICON_DIR + 'qm.png'), 'Load QM', self)
        #load_qm_button.setShortcut('Ctrl+L')
        load_qm_button.setStatusTip('Load QM files')
        load_qm_button.triggered.connect(self.load_files)
        fileMenu.addAction(load_qm_button)
        #Load VL button
        load_vl_button = QtWidgets.QAction(QtGui.QIcon(ICON_DIR + 'vl.png'), 'Load VL', self)
        #load_vl_button.setShortcut('Ctrl+L')
        load_vl_button.setStatusTip('Load VL files')
        load_vl_button.triggered.connect(self.load_files_vl)
        fileMenu.addAction(load_vl_button)
        #Load XY button
        load_xy_button = QtWidgets.QAction(QtGui.QIcon(ICON_DIR + 'xy.png'), 'Load XY', self)
        #load_xy_button.setShortcut('Ctrl+L')
        load_xy_button.setStatusTip('Load XY files')
        load_xy_button.triggered.connect(self.load_files_xy)
        fileMenu.addAction(load_xy_button)
        #
        fileMenu.addSeparator()
        #Remove item button
        remove_item_button = QtWidgets.QAction(QtGui.QIcon(ICON_DIR + 'remove.png'), 'Remove item', self)
        #load_xy_button.setShortcut('Ctrl+L')
        remove_item_button.setStatusTip('Remove item')
        remove_item_button.triggered.connect(self.unload)
        fileMenu.addAction(remove_item_button)
        #
        fileMenu.addSeparator()
        #
        open_project_button = QtWidgets.QAction(QtGui.QIcon(ICON_DIR + 'open.png'), 'Open project', self)
        open_project_button.setShortcut('Ctrl+O')
        open_project_button.setStatusTip('Open project')
        open_project_button.triggered.connect(self.open_project)
        fileMenu.addAction(open_project_button)
        #
        save_project_button = QtWidgets.QAction(QtGui.QIcon(ICON_DIR + 'save.png'), 'Save project', self)
        save_project_button.setShortcut('Ctrl+S')
        save_project_button.setStatusTip('Save project')
        save_project_button.triggered.connect(self.save_project)
        fileMenu.addAction(save_project_button)
        #
        fileMenu.addSeparator()
        #Exit button
        exitButton = QtWidgets.QAction(QtGui.QIcon(ICON_DIR + 'exit.png'), 'Exit', self)
        exitButton.setShortcut('Ctrl+Q')
        exitButton.setStatusTip('Exit application')
        exitButton.triggered.connect(self.close)
        fileMenu.addAction(exitButton)

        ## Export menu
        exportMenu = mainMenu.addMenu('Export')
        #Export lines button
        export_spectra_Button = QtWidgets.QAction(QtGui.QIcon(ICON_DIR + 'e_spectra.png'), 'Export spectra', self)
        #load_qm_button.setShortcut('Ctrl+L')
        export_spectra_Button.setStatusTip('Analyse the peaks and widhts of the spectra')
        export_spectra_Button.triggered.connect(self.export_lines)
        #Export Vlines button
        export_vlines_Button = QtWidgets.QAction(QtGui.QIcon(ICON_DIR + 'e_vl.png'), 'Export VLines', self)
        #load_qm_button.setShortcut('Ctrl+L')
        export_vlines_Button.setStatusTip('Analyse the peaks and widhts of the spectra')
        export_vlines_Button.triggered.connect(self.export_vlines)
        exportMenu.addAction(export_spectra_Button)
        exportMenu.addAction(export_vlines_Button)

        ## Tools menu
        toolsMenu = mainMenu.addMenu('Tools')
        #Peaks Analyses button
        pw_Button = QtWidgets.QAction(QtGui.QIcon(ICON_DIR + 'pa.png'), 'Peaks Analyses', self)
        #load_qm_button.setShortcut('Ctrl+L')
        pw_Button.setStatusTip('Analyse the peaks and widhts of the spectra')
        pw_Button.triggered.connect(self.analyses_pw)
        toolsMenu.addAction(pw_Button)
        #Average Analyses button
        averaged_Button = QtWidgets.QAction(QtGui.QIcon(ICON_DIR + 'av.png'), 'Averaged Analyses', self)
        #load_qm_button.setShortcut('Ctrl+L')
        averaged_Button.setStatusTip('Analyse the averaged and standard deviation from transitions')
        averaged_Button.triggered.connect(self.analyses_av)
        toolsMenu.addAction(averaged_Button)
        
        ## About menu
        aboutMenu = mainMenu.addMenu('About')
        #Load QM button
        info_Button = QtWidgets.QAction(QtGui.QIcon(ICON_DIR + 'muspec.png'), 'Info', self)
        #load_qm_button.setShortcut('Ctrl+L')
        info_Button.setStatusTip('Analyse the peaks and widhts of the spectra')
        info_Button.triggered.connect(self.window_info)
        aboutMenu.addAction(info_Button)


    def create_grid(self):
        self.frame = QtWidgets.QFrame(self)
        self.grid = QtWidgets.QGridLayout(self.frame)

        #Creating the list of load files
        self.filesViewGroupBox = QtWidgets.QGroupBox("Items:")
        layout = QtWidgets.QHBoxLayout()
        self.filesView = QtWidgets.QTableWidget(0, 2)
        self.filesView.setHorizontalHeaderLabels(("Item Name", "ID"))
        self.filesView.verticalHeader().hide()
        self.filesView.hideColumn(1)
        self.filesView.setColumnWidth(0, 220)
        self.filesView.setShowGrid(False)
        layout.addWidget(self.filesView)
        self.filesViewGroupBox.setLayout(layout)
        self.grid.addWidget(self.filesViewGroupBox, 0, 0)
        
        #Creating the general box
        self.generalViewGroupBox = QtWidgets.QGroupBox("General:")
        layout_general = QtWidgets.QGridLayout(self.generalViewGroupBox)
        ##Normalize
        self.setNormalization = QtWidgets.QComboBox(self)
        self.setNormalizationLabel = QtWidgets.QLabel("Set normalization:")
        layout_general.addWidget(self.setNormalizationLabel, 0, 0)
        layout_general.addWidget(self.setNormalization, 0, 1)
        ##Energy units
        self.setEnergyUnit = QtWidgets.QComboBox(self)
        self.setEnergyUnitLabel = QtWidgets.QLabel("Energy unit:")
        self.fill_energy_unit()
        layout_general.addWidget(self.setEnergyUnitLabel, 1, 0)
        layout_general.addWidget(self.setEnergyUnit, 1, 1)
        ##buttom
        self.update_bt = QtWidgets.QPushButton('Draw', self)
        layout_general.addWidget(self.update_bt)
        ##Set grid
        self.generalViewGroupBox.setLayout(layout_general)
        self.grid.addWidget(self.generalViewGroupBox, 2, 0)
        #Creating figure
        self.create_figure()
        self.grid.addWidget(self.canvas, 0, 1, 3, 1)
        #Creating log
        logTextBox = QTextEditLogger(self)
        # You can format what is printed to text box
        logTextBox.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
        logging.getLogger().addHandler(logTextBox)
        # You can control the logging level
        logging.getLogger().setLevel(logging.INFO)
        # Add the new logging box widget to the layout
        self.grid.addWidget(logTextBox.widget, 3, 0, 1, 2)

        #Creating the info box
        self.create_info_vl()
        self.create_info_xy()

        self.analyses_window = QTextAnalyses('Peaks Analyses')
        self.analyses_averaged = QTextAnalyses('Averaged Analyses')
        
    def create_info_vl(self):
        self.infoViewGroupBox = QtWidgets.QGroupBox("Item options:")
        layout = QtWidgets.QGridLayout(self.infoViewGroupBox)
        nl = 0
        ##Broadening function
        self.broadening = QtWidgets.QComboBox(self)
        self.broadeningLabel = QtWidgets.QLabel("Function:")
        self.broadening.addItem('Lorentzian')
        self.broadening.addItem('Gaussian')
        layout.addWidget(self.broadeningLabel, nl, 0)
        layout.addWidget(self.broadening, nl, 1)
        nl += 1
        ##Procedures
        self.procedures = QtWidgets.QComboBox(self)
        self.proceduresLabel = QtWidgets.QLabel("Procedure:")
        self.fill_procedures()
        layout.addWidget(self.proceduresLabel, nl, 0)
        layout.addWidget(self.procedures, nl, 1)
        nl += 1
        ##Groups
        self.groupsLineEdit = QtWidgets.QLineEdit()
        self.groupsLabel = QtWidgets.QLabel("Groups:")
        self.groupsLabel.setBuddy(self.groupsLineEdit)
        layout.addWidget(self.groupsLabel, nl, 0)
        layout.addWidget(self.groupsLineEdit, nl, 1)
        nl += 1
        ##FWHM
        self.fwhmLineEdit = QtWidgets.QLineEdit()
        self.fwhmLabel = QtWidgets.QLabel("FWHM:")
        self.fwhmLabel.setBuddy(self.fwhmLineEdit)
        layout.addWidget(self.fwhmLabel, nl, 0)
        layout.addWidget(self.fwhmLineEdit, nl, 1)
        nl += 1
        ##XY MULTIPLICATION
        self.yLineEdit = QtWidgets.QLineEdit()
        self.yLabel = QtWidgets.QLabel("X;Y multiplier:")
        self.yLabel.setBuddy(self.yLineEdit)
        layout.addWidget(self.yLabel, nl, 0)
        layout.addWidget(self.yLineEdit, nl, 1)
        nl += 1
        ##XRANGE 
        self.xrangeLineEdit = QtWidgets.QLineEdit()
        self.xrangeLabel = QtWidgets.QLabel("X_init;X_end:")
        self.xrangeLabel.setBuddy(self.xrangeLineEdit)
        layout.addWidget(self.xrangeLabel, nl, 0)
        layout.addWidget(self.xrangeLineEdit, nl, 1)
        nl += 1
        ##Ocultar ID
        self.plot_hide = QtWidgets.QCheckBox('Hide', self)
        layout.addWidget(self.plot_hide)
        ##Vlines
        self.vlines = QtWidgets.QCheckBox('Vertical lines', self)
        layout.addWidget(self.vlines)
        ##Plot groups
        self.plot_groups = QtWidgets.QCheckBox('Plot grous', self)
        layout.addWidget(self.plot_groups)
        ##Plot color groups
        self.plot_color_groups = QtWidgets.QCheckBox('Plot vlines groups', self)
        layout.addWidget(self.plot_color_groups)
        self.grid.addWidget(self.infoViewGroupBox, 1, 0)
        self.infoViewGroupBox.hide()

    def create_info_xy(self):
        self.info_xy_GroupBox = QtWidgets.QGroupBox("XY Item:")
        layout = QtWidgets.QGridLayout(self.info_xy_GroupBox)
        nl = 0
        ##Columns
        self.xyColumnsLineEdit = QtWidgets.QLineEdit()
        self.xyColumnsLabel = QtWidgets.QLabel("XY:")
        self.xyColumnsLabel.setBuddy(self.xyColumnsLineEdit)
        layout.addWidget(self.xyColumnsLabel, nl, 0)
        layout.addWidget(self.xyColumnsLineEdit, nl, 1)
        nl += 1
        ##Ymult
        self.yLineXYEdit = QtWidgets.QLineEdit()
        self.yLabelXY = QtWidgets.QLabel("X;Y multiplier:")
        self.yLabelXY.setBuddy(self.yLineXYEdit)
        layout.addWidget(self.yLabelXY, nl, 0)
        layout.addWidget(self.yLineXYEdit, nl, 1)
        nl += 1
        ##Ocultar ID
        self.plot_hide_xy = QtWidgets.QCheckBox('Hide', self)
        layout.addWidget(self.plot_hide_xy)

        self.grid.addWidget(self.info_xy_GroupBox, 1, 0)
        self.info_xy_GroupBox.hide()

    def create_figure(self):
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        self.axes = self.fig.add_subplot(111)
        self.axes.set_ylabel('Intensity')

    def create_toolbar(self):
        self.fileToolBar = self.addToolBar("File")
        self.fileToolBar.addAction(self.loadAct)
        self.fileToolBar.addAction(self.load_vl)
        self.fileToolBar.addAction(self.load_xy)
        self.fileToolBar.addAction(self.remoAct)

        self.testToolBar = self.addToolBar("Canvas")

        self.pltToolBar = NavigationToolbar(self.canvas, self, coordinates=True)
        self.testToolBar.addWidget(self.pltToolBar)

    def create_toolbar_analyses(self):
    	pass

    def create_actions(self):
        #load file
        self.loadAct = QtWidgets.QAction(QtGui.QIcon(ICON_DIR + 'qm.png'),"Load", self,
                statusTip="Import files",
                triggered=self.load_files)

        #load vl
        self.load_vl = QtWidgets.QAction(QtGui.QIcon(ICON_DIR + 'vl.png'),"Load VL", self,
                statusTip="Import files",
                triggered=self.load_files_vl)

        #load xy
        self.load_xy = QtWidgets.QAction(QtGui.QIcon(ICON_DIR + 'xy.png'),"Load XY", self,
                statusTip="Import files",
                triggered=self.load_files_xy)
        
        #unload file
        self.remoAct = QtWidgets.QAction(QtGui.QIcon(ICON_DIR + 'remove.png'),"Remove files", self,
                statusTip="Remove file.",
                triggered=self.unload)

        #lines edit
        self.groupsLineEdit.editingFinished.connect(self.update_group)
        self.fwhmLineEdit.editingFinished.connect(self.update_fwhm)
        self.yLineEdit.editingFinished.connect(self.validate_ymult)
        self.yLineXYEdit.editingFinished.connect(self.validate_ymult)

        #button
        self.update_bt.clicked.connect(self.update_item)

        #change item table
        self.filesView.itemClicked.connect(self.show_info)

        #change broadening functions
        self.broadening.currentIndexChanged.connect(self.update_broadening)

        #change procedures functions
        self.procedures.currentIndexChanged.connect(self.update_procedures)

        #checkboxes
        self.vlines.stateChanged.connect(self.update_vlines)
        self.plot_groups.stateChanged.connect(self.update_plot_groups)
        self.plot_color_groups.stateChanged.connect(self.update_plot_color_groups)
        self.plot_hide.stateChanged.connect(self.update_plot_hide)
        self.plot_hide_xy.stateChanged.connect(self.update_plot_hide)

    def new_item(self, item):
        row = self.filesView.rowCount()
        self.row_id = self.item_id
        self.item[self.item_id] = item
        self.item[self.item_id].item_id = self.item_id
        self.item[self.item_id].ax = self.axes

        fileNameItem = QtWidgets.QTableWidgetItem(self.item[self.item_id].label)
        fileID = QtWidgets.QTableWidgetItem(str(self.item_id))
        self.filesView.insertRow(row)
        self.filesView.setItem(row, 0, fileNameItem)
        self.filesView.setItem(row, 1, fileID)
        self.item_id += 1        

    def which_program(self, file):
        f = open(file, 'r')
        for n in range(1000):
            line = f.readline()
            if 'gaussian' in line.lower():
                f.close()
                return Gaussian.Initialize(file)
            elif 'dalton' in line.lower():
                f.close()
                return Dalton.Initialize(file) 
        logging.error('The program was not identified.')

    def load_files_qm(self, files):
        files_init  = [self.which_program(x) for x in files if self.which_program(x) is not None]
        files_label = [os.path.basename(x) for x in files]

        energies = {}
        intensities = {}

        label = ''
        nmax = len(files_init)
        for n in range(nmax):
            if len(files_init[n].calc) == 0:
                logging.error('The file {0} has no identified spectrum.'.format(files_init[n].name))
                continue
            files_init[n].get_all_data()
            check = [len(x) for x in files_init[n].intensities.values()]
            if len(check) == 0:
                logging.error('The file {0} has no energies and will be removed.'.format(files_init[n].name))
                continue
            ch_min = min(check)
            check = all(x == ch_min for x in check)
            if not check:
                logging.warning('The number of energies are not the same for all transitions.')

            for key in files_init[n].energies.keys():
                nval = len(files_init[n].energies[key])
                for nval in range(nval):
                    Tools.set_key(self, energies, key, files_init[n].energies[key][nval])
                    Tools.set_key(self, intensities, key, files_init[n].intensities[key][nval])
                files_init[n].energies[key] = np.array(files_init[n].energies[key])
                files_init[n].intensities[key] = np.array(files_init[n].intensities[key])

            label += files_label[n] + ','

        if len(energies) == 0:
            logging.error('There is no file loaded.')
            return
        
        folder = os.path.dirname(files[0])
        logging.info('Opening files from folder {0}'.format(folder))
        logging.info(','.join(files_label))
        item = VLItem(energies, intensities, label)
        self.new_item(item)
        item.draw()
        self.canvas.draw()

    def load_files(self):
        files = QtWidgets.QFileDialog.getOpenFileNames(self, 'Open file', self.last_dir)
        files = files[:][0]
        if len(files) == 0:
            return
        self.last_dir = os.path.dirname(files[0])
        self.load_files_qm(files)
        
    def load_files_xy(self):
        file = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', self.last_dir)
        file_name = file[0]
        if len(file_name) == 0:
            return
        item = XYItem(file_name)
        if item.create_status:
            self.new_item(item)
            self.last_dir = os.path.dirname(file_name)
        else:
            logging.error('The file {0} cannot be read using numpy.loadtxt().'.format(file_name))
            return
        
        logging.info('Opening file from folder {0}'.format(self.last_dir))
        logging.info(item.label)
        logging.info('There are {0} columns in the file.'.format(len(item.data)))
        item.draw()
        
    def load_files_vl(self):
        file = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', self.last_dir)
        file_name = file[0]
        if len(file_name) == 0:
            return

        energies = {}
        intensities = {}
        label = os.path.basename(file_name)

        comments = ('#', '!', '@')
        try:
            data = np.loadtxt(file_name, comments=comments).T
        except:
            logging.error('The file {0} cannot be read using numpy.loadtxt().'.format(file_name))
            return
        nmax = len(data)
        if nmax < 1:
            logging.error('Less than two columns. The file {0} will not be loaded.'.format(file_name))
            return
        mmax = len(data[0])
        if mmax < 1:
            logging.error('Less than one line in file {0}, it will not be loaded'.format(file_name))
            return

        energies = {}
        intensities = {}

        for m in range(mmax):
            if nmax == 2:
                state = m + 1
                n = 0
            elif nmax == 3:
                state = int(data[0][m])
                n = 1
            elif nmax == 4:
                state = int(data[0][m])
                n = 1
            energy = data[n][m]
            intensity = data[-1][m]
            Tools.set_key(self, energies, state, energy)
            Tools.set_key(self, intensities, state, intensity)

        logging.info('The last column was set for the intensities.')
        check = [len(x) for x in intensities.values()]
        ch_min = min(check)
        check = all(x == ch_min for x in check)
        if not check:
            logging.warning('The number of energies are not the same for all transitions.')

        logging.info('Opening files from folder {0}'.format(self.last_dir))
        logging.info(label)
        item = VLItem(energies, intensities, label)
        self.new_item(item)
        item.draw()
        self.canvas.draw()

    def unload(self):
        if self.row == -1:
            return
        self.infoview.hide()
        self.item[self.row_id].remove()
        self.filesView.removeRow(self.row)
        del self.item[self.row_id]
        self.row = -1
        self.fill_normalization_cbox()
        self.canvas.draw()
                
    def fill_energy_unit(self):
        self.setEnergyUnit.addItem('eV')
        self.setEnergyUnit.addItem('cm-1')
        self.setEnergyUnit.addItem('nm')

    def fill_normalization_cbox(self):
        self.setNormalization.clear()
        self.setNormalization.addItem('False', (None, None))
        g = [ self.item[n].ngroups for n in self.item.keys() ]
        if len(g) > 0:
            for n in range(min(g)):
                text = 'g_{0}'.format(n + 1)
                self.setNormalization.addItem(text, ('all', n))
        for k, value in self.item.items():
            for n in range(value.ngroups):
                text = '{0}_{1}'.format(value.label, n + 1)
                self.setNormalization.addItem(text, (k, n))

    def fill_procedures(self):
        procedures = {
        r'Procedure 1: Convolution over the average energy and intensity using a given $\Gamma$.':'procedure_01',
        r'Procedure 2: Convolution over the average energy and intensity using the standard deviation as $\Gamma$.':'procedure_02',
        r'Procedure 3: Convolution over all configurations using a given $\Gamma$.':'procedure_03',
        r'Procedure 4: Convolution over all configurations using the standard deviation as $\Gamma$.':'procedure_04',
        r'Procedure 5: Convolution over all configurations fitting $\Gamma$.':'procedure_05',
        }
        for k, item in procedures.items():
            self.procedures.addItem(item)

    def update_nomalization(self, row_id):
        item_id, g_id = self.setNormalization.currentData()
        if item_id is None:
            self.y_norm = 1.0
        elif item_id != 'all':
            self.y_norm = np.max(self.item[item_id].spectra_y[g_id]) * self.item[item_id].ymult
        elif item_id == 'all':
            self.y_norm = np.max(self.item[row_id].spectra_y[g_id]) * self.item[row_id].ymult

    def update_vlines(self):
        if self.row == -1:
            return
        if self.vlines.isChecked():
            state = self.item[self.row_id].vline
            if state:
                return
            self.item[self.row_id].vline = True
        else:
            self.item[self.row_id].vline = False    

    def update_plot_groups(self):
        if self.row == -1:
            return
        if self.plot_groups.isChecked():
            state = self.item[self.row_id].plot_group
            if state:
                return
            self.item[self.row_id].plot_group = True
        else:
            self.item[self.row_id].plot_group = False

    def update_plot_color_groups(self):
        if self.row == -1:
            return
        if self.plot_color_groups.isChecked():
            state = self.item[self.row_id].plot_color_group
            if state:
                return
            self.item[self.row_id].plot_color_group = True
        else:
            self.item[self.row_id].plot_color_group = False

    def update_plot_hide(self):
        name = self.item[self.row_id].__class__.__name__
        if name == 'VLItem':
            checkbox = self.plot_hide
        elif name == 'XYItem':
            checkbox = self.plot_hide_xy

        if self.row == -1:
            return
        if checkbox.isChecked():
            state = self.item[self.row_id].plot_hides
            if state:
                return
            self.item[self.row_id].plot_hides = True
        else:
            self.item[self.row_id].plot_hides = False

    def update_broadening(self):
        self.item[self.row_id].function = self.broadening.currentText()
        self.item[self.row_id].broad_func = self.item[self.row_id].broad_bib[self.item[self.row_id].function]

    def update_procedures(self):
        self.item[self.row_id].procedure = self.procedures.currentText()

    def update_draw(self):
        self.clear_draw()
        for n in self.item.keys():
            self.item[n].draw()
        self.axes.relim()
        self.canvas.draw()
        
    def update_group(self):
        self.validate_groups()
        self.item[self.row_id].update_par(('group',self.groupsLineEdit.text()))
        self.groupsLineEdit.setText(self.groups_to_string(self.item[self.row_id].groups))
        self.fill_normalization_cbox()
        nfwhm = len(self.item[self.row_id].fwhm)
        if nfwhm < self.item[self.row_id].ngroups:
            text = self.fwhmLineEdit.text()
            text += ';0.4' * (self.item[self.row_id].ngroups - nfwhm)
            self.fwhmLineEdit.setText(text)
        elif nfwhm > self.item[self.row_id].ngroups:
            text = self.fwhmLineEdit.text().split(';')
            text = ';'.join(text[0:self.item[self.row_id].ngroups])
            self.fwhmLineEdit.setText(text)

    def update_fwhm(self):
        self.validate_fwhm()
        self.item[self.row_id].update_par(('fwhm',self.fwhmLineEdit.text()))

    def update_item(self):
        if self.row == -1:
            return

        self.item[self.row_id].update_item()
        self.update_draw()
        self.show_info()

    def clear_draw(self):
        for n in self.item.keys():
            self.item[n].remove_draw()
        self.axes.relim()
        self.fig.gca().set_prop_cycle(None)

    def set_energy_unit(self, x):
        if self.setEnergyUnit.currentText() == 'cm-1':
            self.axes.set_xlabel('Energy $(cm^{-1})$')
            return x * EV2CM
        elif self.setEnergyUnit.currentText() == 'nm':
            self.axes.set_xlabel('Wavelength (nm)')
            return EV2NM / x
        elif self.setEnergyUnit.currentText() == 'eV':
            self.axes.set_xlabel('Energy (eV)')
            return x

    def groups_to_string(self, group):
        s = ''
        for list_gr in group:
            gr = np.split(list_gr, [i+1 for i,j in enumerate(np.diff(list_gr)) if j > 1])
            for n in gr:
                if len(n) > 1:
                    s += '{0}-{1},'.format(n[0],n[-1])
                else:
                    s += '{0},'.format(n[0])
            s += ';'
        return s.replace(',;',';')[:-1]

    def fwhm_to_string(self, fwhm):
        s = ''
        for fwhm_ in fwhm:
            s += '{0};'.format(fwhm_)
        return s[:-1]

    def fill_item(self):
        index_function = self.broadening.findText(self.item[self.row_id].function)
        if index_function >= 0:
             self.broadening.setCurrentIndex(index_function)

        index_function = self.procedures.findText(self.item[self.row_id].procedure)
        if index_function >= 0:
             self.procedures.setCurrentIndex(index_function)

        self.infoview = self.infoViewGroupBox
        groups = self.item[self.row_id].groups
        groups_str = self.groups_to_string(groups)
        self.groupsLineEdit.setText(groups_str)

        fwhm = self.item[self.row_id].fwhm
        fwhm_str = self.fwhm_to_string(fwhm)
        self.fwhmLineEdit.setText(fwhm_str)

        self.yLineEdit.setText('{0};{1}'.format(self.item[self.row_id].xmult, self.item[self.row_id].ymult))
        
        if self.item[self.row_id].vline:
            self.vlines.setChecked(True)
        else:
            self.vlines.setChecked(False)

        if self.item[self.row_id].plot_group:
            self.plot_groups.setChecked(True)
        else:
            self.plot_groups.setChecked(False)        

        if self.item[self.row_id].plot_hides:
            self.plot_hide.setChecked(True)
        else:
            self.plot_hide.setChecked(False)

        if self.item[self.row_id].plot_color_group:
            self.plot_color_groups.setChecked(True)
        else:
            self.plot_color_groups.setChecked(False)
        
        self.fill_normalization_cbox()
        self.ymult_field = self.yLineEdit
        self.xrangeLineEdit.setText('{0};{1}'.format(*self.item[self.row_id].lim))
        self.infoview.show()

    def fill_xy(self):
        self.infoview = self.info_xy_GroupBox
        self.xyColumnsLineEdit.setText(self.item[self.row_id].xyColumnsLabelString)
        self.yLineXYEdit.setText('{0};{1}'.format(self.item[self.row_id].xmult, self.item[self.row_id].ymult))
        self.ymult_field = self.yLineXYEdit

        if self.item[self.row_id].plot_hides:
            self.plot_hide_xy.setChecked(True)
        else:
            self.plot_hide_xy.setChecked(False)

        self.fill_normalization_cbox()
        self.infoview.show()

    def show_info(self):
        self.row = self.filesView.currentRow()
        self.row_id = int(self.filesView.item(self.row, 1).text())

        self.infoview.hide()
        name = self.item[self.row_id].__class__.__name__
        self.fill[name]()

    def validate_groups(self):
        text = self.groupsLineEdit.text()
        text = text.replace('-',';').replace(',',';')
        text = text.split(';')
        for string in text:
            if not string.isdigit():
                t = '1-{0}'.format(self.item[self.row_id].nstates)
                logging.warning('The entry {0} is not an integer number. Setting groups to {1}.'.format(string, t))
                self.groupsLineEdit.setText(t)

    def validate_fwhm(self):
        text = self.fwhmLineEdit.text()
        text = text.split(';')
        nmax = len(text)
        for n in range(nmax):
            try:
                float(text[n])
            except ValueError:
                logging.warning('The entry {0} is not a float number. Setting FWHM to 0.4.'.format(text[n]))
                text[n] = '0.4'
        text = ';'.join(text)
        self.fwhmLineEdit.setText(text)

    def validate_xrange(self):
        text = self.xrangeLineEdit.text()
        text = text.split(';')
        nmax = len(text)
        if nmax != 2:
            logging.warning('The number of entries is not equal to 2. Setting the defaults values: 2;8.')
            text = ['2', '8']
        try:
            self.item[self.row_id].lim[0] = float(text[0])
        except:
            logging.warning('The entry {0} cannot be converted to float. Setting as 2.'.format(text[0]))
            self.item[self.row_id].lim[0] = 2

        try:
            self.item[self.row_id].lim[1] = float(text[1])
        except:
            logging.warning('The entry {0} cannot be converted to float. Setting as 8.'.format(text[1]))
            self.item[self.row_id].lim[1] = 8

    def validate_ymult(self):
        text = self.ymult_field.text().split(';')
        if len(text) == 2:
            try:
                self.item[self.row_id].xmult = float(text[0])
            except:
                logging.warning('The X entry {0} is not a float number. Setting xmult to 1.0'.format(text[0]))
                text[0] = 1.0
                self.item[self.row_id].xmult = text[0]
                self.ymult_field.setText('{0};{1}'.format(text[0], text[1]))
            try:
                self.item[self.row_id].ymult = float(text[1])
            except:
                logging.warning('The Y entry {0} is not a float number. Setting ymult to 1.0'.format(text[1]))
                text[1] = 1.0
                self.item[self.row_id].ymult = text[1]
                self.ymult_field.setText('{0};{1}'.format(text[0], text[1]))
        else:
            logging.warning('The XY multipliers do not match the X;Y format. Setting xmult and ymult to 1.0')
            self.item[self.row_id].xmult = 1.0
            self.item[self.row_id].ymult = 1.0
            self.ymult_field.setText('1.0;1.0')

    def validate_xy(self):
        text = self.xyColumnsLineEdit.text().split(';')
        nplots = len(text)
        nx = []
        ny = []

        for n in range(nplots):
            tmp = text[n].split(',')
            if len(tmp) != 2:
                print('Error: Defining the xy columns')
                tmp = [0, 1]
            nx += [tmp[0]]
            ny += [tmp[1]]

        self.item[self.row_id].xcolumns = np.array(nx, int)
        self.item[self.row_id].ycolumns = np.array(ny, int)
        self.item[self.row_id].xyColumnsLabelString = self.xyColumnsLineEdit.text()

    def export_lines(self):
        text = ''
        text_label = ''
        data = []
        for row_id in self.item.keys():
            name = self.item[row_id].__class__.__name__
            if name != 'VLItem':
                continue

            if text != '':
                text += '   '
                text_label += '   '

            text_label += '{0:12s}'.format(' Label')
            text += '{0:>12s}'.format('Energy (eV)')

            text_label += '{0:>15s}'.format(self.item[row_id].label[:14])*self.item[row_id].ngroups
            text_groups = self.groups_to_string(self.item[row_id].groups)
            text_tmp = ['{0:>15s}'.format(x) for x in text_groups.split(';')]
            text += ''.join(text_tmp)
            
            data += [self.item[row_id].get_x(self.item[row_id].spectra_x)]
            data += [self.item[row_id].get_y(self.item[row_id].spectra_y[n]) for n in range(self.item[row_id].ngroups)]
            if self.item[row_id].ngroups > 1:
                text_label += '{0:>15s}'.format(self.item[row_id].label[:14])
                text += '{0:>15s}'.format(text_groups)
                data += [np.sum(self.item[row_id].get_y(self.item[row_id].spectra_y), axis=0)]
            
        data = np.asarray(data)
        header = '{0}\n{1}'.format(text_label, text)
        name = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File')
        if len(name[0]) != 0:
            np.savetxt(name[0], data.T, fmt='%14.6f', header=header)

    def export_vlines(self):
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select Directory')
        if folder == '':
            return
        for row_id in self.item.keys():
            name = self.item[row_id].__class__.__name__
            if name != 'VLItem':
                continue
            s = self.groups_to_string(self.item[row_id].groups).split(';')
            for n in range(len(s)):
                vl_name = '{0}_g_{1}.txt'.format(self.item[row_id].label, s[n])
                text  = ' Vertical lines for item: {0}.\n'.format(vl_name)
                text += ' Energy (eV)      Intensity'

                data = np.stack((self.item[row_id].vlines_x[n], self.item[row_id].vlines_y[n].T))

                logging.info('Saving the file: {0} in folder: {1}.'.format(vl_name, folder))
                file = os.path.join(folder, vl_name)
                np.savetxt(file, data.T, fmt='%14.6f', header=text)

    def analyses_av(self):
        text = ''
        for row_id in self.item.keys():
            name = self.item[row_id].__class__.__name__
            if name != 'VLItem':
                continue

            text += 'Averages for {0}.\n\n'.format(self.item[row_id].label)
            text += 'DE       - Average transition energy.\n'
            text += 'STD_DE   - Standard deviation transition energy.\n'
            text += 'Int      - Average intensity.\n'
            text += 'STD_int  - Standard deviation intensity.\n'
            text += 'DE*      - Average transition energy by intensity.\n'
            text += 'STD_DE*  - Standard deviation transition energy by intensity.\n'
            text += 'Int*     - Average intensity by intensity.\n'
            text += 'STD_int* - Standard deviation intensity by intensity.\n'
            text += '\n'
            
            header = 'State DE STD_DE Int STD_Int DE* STD_DE* Int* STD_Int*'.split()
            text += '{0:>6}{1:>8}{2:>8}{3:>8}{4:>8}{5:>8}{6:>8}{7:>8}{8:>8}\n'.format(*header)


            fmt = '{0:>6}{1:>8.4f}{2:>8.4f}{3:>8.4f}{4:>8.4f}{5:>8.4f}{6:>8.4f}{7:>8.4f}{8:>8.4f}\n'

            for key in self.item[row_id].energies.keys():

                text += fmt.format(
                                    key,
                                    np.average(self.item[row_id].energies[key]),
                                    np.std(self.item[row_id].energies[key]),
                                    np.average(self.item[row_id].intensities[key]),
                                    np.std(self.item[row_id].intensities[key]),
                                    np.average(self.item[row_id].energies[key], weights=self.item[row_id].intensities[key]),
                                    Tools.weighted_std(self, self.item[row_id].energies[key], self.item[row_id].intensities[key]),
                                    np.average(self.item[row_id].intensities[key], weights=self.item[row_id].intensities[key]),
                                    Tools.weighted_std(self, self.item[row_id].intensities[key], self.item[row_id].intensities[key])
                                  )
            text += '\n\n'
            self.window_analyses_average(text)

    def analyses_pw(self):
        text = ''
        for row_id in self.item.keys():
            name = self.item[row_id].__class__.__name__
            if name != 'VLItem':
                continue

            text += 'Peaks analyses for {0}.\n'.format(self.item[row_id].label)
            x = self.item[row_id].get_x(self.item[row_id].spectra_x)
            if self.item[row_id].procedure == 'procedure_05':
                fwhm = self.item[row_id].fwhm_optimal
            else:
                fwhm = self.item[row_id].fwhm

            for n in range(self.item[row_id].ngroups):
                text += 'List of peaks for group {0}.\n'.format(n + 1)
                text += 'States: {0}\n'.format(self.groups_to_string(self.item[row_id].groups).split(';')[n])
                text += 'FWHM = {0:0.2f}\n'.format(fwhm[n])
                y = self.item[row_id].get_y(self.item[row_id].spectra_y[n])
                text += self.list_peaks(x, y) + '\n'

            if self.item[row_id].ngroups > 1:
                text += 'List of peaks of the total spectra.\n'
                y = np.sum(self.item[row_id].get_y(self.item[row_id].spectra_y), axis=0)
                text += self.list_peaks(x, y)
            text += '\n'
        self.window_analyses(text)

    def list_peaks(self, x, y):
        from scipy.signal import find_peaks
        indices, others_info = find_peaks(y, width=0.5, prominence=np.max(y) * 0.05)
        dx = [x[int(others_info['right_ips'][n])] - x[int(others_info['left_ips'][n])] for n in range(len(indices))]
        peaks = ['( {0:.2f} , {1:.2f} ) - fwhm = {2:.2f}'.format(x[indices[n]], y[indices[n]], dx[n]) for n in range(len(indices))]
        return '\n'.join(peaks)

    def window_info(self):
        message  = 'Author: Tárcius N. Ramos\n'
        message += 'Institution: Universidade de Sao Paulo, Instituto de Fisica\n'
        message += 'Funding: FAPESP 2015/14189-3\n'
        message += '\nVersion: 1.0\n'
        message += '\nContact: tarcius_ramos@hotmail.com\n'
        info = QtWidgets.QMessageBox()
        info.setWindowTitle('About')
        info.setText(message)
        info.setWindowIcon(QtGui.QIcon(ICON_DIR + 'muspec.png'))
        info.exec_()

    def window_analyses(self, message):
        self.analyses_window.emit(message)
        self.analyses_window.show()

    def window_analyses_average(self, message):
        self.analyses_averaged.emit(message)
        self.analyses_averaged.show()

    def save_project(self):
        file = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File')
        file = file[0]
        if len(file) == 0:
            return
        f = open(file, 'w')
        for row_id in self.item.keys():
            settings = dict()
            settings['item'] = self.item[row_id].__class__.__name__
            if settings['item'] == 'VLItem':
                settings['label'] = self.item[row_id].label
                settings['function'] = self.item[row_id].function
                settings['procedure'] = self.item[row_id].procedure
                settings['groups'] = self.item[row_id].groups.tolist()
                settings['FWHM'] = self.item[row_id].fwhm
                settings['X_Y_multiplier'] = [self.item[row_id].xmult, self.item[row_id].ymult]
                settings['X_init_X_final'] =self.item[row_id].lim
                settings['hide'] = self.item[row_id].plot_hides
                settings['vertical_lines'] = self.item[row_id].vline
                settings['plot_groups'] = self.item[row_id].plot_group
                settings['plot_vlines_groups'] = self.item[row_id].plot_color_group
                settings['energies'] = self.item[row_id].energies
                settings['intensities'] = self.item[row_id].intensities
            if settings['item'] == 'XYItem':
                settings['label'] = self.item[row_id].label
                settings['ngroups'] = self.item[row_id].ngroups
                settings['xcolumns'] = self.item[row_id].xcolumns.tolist()
                settings['ycolumns'] = self.item[row_id].ycolumns.tolist()
                settings['xyColumnsLabelString'] = self.item[row_id].xyColumnsLabelString
                settings['xmult'] = self.item[row_id].xmult
                settings['ymult'] = self.item[row_id].ymult
                settings['plot_hides'] = self.item[row_id].plot_hides
                settings['data'] = self.item[row_id].data.tolist()
                settings['ncols'] = self.item[row_id].ncols
            f.write(str(settings) + '\n')
        f.close()

    def open_project(self):
        files = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', self.last_dir)
        file = files[0]
        if len(file) == 0:
            return
        self.last_dir = os.path.dirname(files[0])

        f = open(file, 'r')
        lines = f. readlines()
        f.close()

        nitems = len(lines)

        for n in range(nitems):
            try:
                data = eval(lines[n])
            except:
                logging.warning('Error reading item number {0}.'.format(n+1))
                return
                print('error')
            if data['item'] == 'VLItem':
                item = VLItem(data['energies'], data['intensities'], data['label'])
                item.fwhm = data['FWHM']
                item.lim = data['X_init_X_final']
                item.xmult = data['X_Y_multiplier'][0]
                item.ymult = data['X_Y_multiplier'][1]
                item.groups = np.asarray(data['groups'])
                item.ngroups = len(item.groups)
                item.vline = data['vertical_lines']
                item.plot_group = data['plot_groups']
                item.plot_color_group = data['plot_vlines_groups']
                item.plot_hides = data['hide']
                item.function = data['function']
                item.broad_func = item.broad_bib[item.function]
                item.procedure = data['procedure']
                self.new_item(item)
                
            if data['item'] == 'XYItem':
                item = XYItem(data['label'])
                item.label = data['label']
                item.ngroups = data['ngroups']
                item.xcolumns = np.asarray(data['xcolumns'])
                item.ycolumns = np.asarray(data['ycolumns'])
                item.xyColumnsLabelString = data['xyColumnsLabelString']
                item.xmult = data['xmult']
                item.ymult = data['ymult']
                item.plot_hides = data['plot_hides']
                item.data = np.asarray(data['data'])
                item.ncols = data['ncols']
                self.new_item(item)

            item.draw()
            self.canvas.draw()

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    w = MainWindow()
    w.show()
    sys.exit(app.exec_())